#!/usr/bin/env python3
"""
coverity_export.py -- export defects and per-defect event traces from the FREE
Coverity Scan service (scan.coverity.com) for offline triage.

This drives the same UNDOCUMENTED internal JSON endpoints that the "View
Defects" web UI calls, reusing your authenticated browser session cookie.  It
can break without notice; it is meant for periodic personal exports.  See
README_export.md (cookie capture, troubleshooting).

Two subcommands:

  defects   Page the defect table into a single JSON file (the "grid export").
            Endpoint (open-source tier): POST /views/table.json with a JSON body
            {"projectId":..,"viewId":..,"pageNum":N}.  Pages are 1-based; the
            tool walks pageNum until every row is collected.

  traces    For chosen CIDs, fetch each defect's event trace (numbered steps
            with file+line) from GET /sourcebrowser/source.json, keyed by the
            fileInstanceId + defectInstanceId carried in the grid export.

  self-test Offline logic checks for both (no network/cookie needed).

Typical flow:

    export COVERITY_COOKIE='COVJSESSIONID-build=..; XSRF-TOKEN=..; isAuthenticated=true; ..'
    python coverity_export.py defects -p 12062 -v 63757 -o exported_defects.json
    python coverity_export.py traces \
        --cids-from exported_defects.json --cids '504003 504028 503971' \
        -o traces.json
"""

import argparse
import json
import math
import os
import sys
import time
from datetime import datetime, timezone

try:
    import requests
except ImportError:
    sys.exit("This script needs 'requests'. Install with:  pip install requests")

# An expired session returns HTTP 200 with the SAML/LDAP login config (valid
# JSON); these keys identify that so we abort instead of saving junk.
AUTH_CHALLENGE_KEYS = ("availableSamlSsoConfigurations", "ldapConfigured",
                       "samlSsoConfigured", "availableLdapConfigurations")

# Defect paging is stateful and split across two endpoints:
#   POST /views/table.json   {projectId,viewId,pageNum}  -> sets the page (no rows)
#   GET  /reports/table.json ?projectId&viewId            -> rows for the set page
DEFAULT_VIEWS_URL = "https://scan6.scan.coverity.com/views/table.json"
DEFAULT_REPORTS_URL = "https://scan6.scan.coverity.com/reports/table.json"

# ----- defects (table.json) field aliases -----------------------------------
RESULTSET_KEY = "resultSet"
ROW_KEYS = ("results", "aaData", "data", "rows")
TOTAL_KEYS = ("totalCount", "overallCount", "iTotalRecords",
              "iTotalDisplayRecords", "recordsTotal", "total")
LIMIT_KEYS = ("limit", "rowCount", "pageSize", "iDisplayLength")
OFFSET_KEYS = ("offset", "firstRow", "startRowOffset")
ID_KEYS = ("cid", "id", "mergedDefectId")

# ----- traces (source.json) endpoint + per-event field aliases --------------
DEFAULT_TRACE_URL = (
    "https://scan6.scan.coverity.com/sourcebrowser/source.json"
    "?projectId={projectId}&fileInstanceId={fileInstanceId}"
    "&defectInstanceId={defectInstanceId}&mergedDefectId={mergedDefectId}"
    "&fileStart=&fileEnd=")
EVENT_FIELD_LINE = ("lineNumber", "line", "lineNo", "strippedLineNumber")
EVENT_FIELD_FILE = ("filePathname", "strippedFilePathname", "file", "fileName",
                    "displayFile", "path", "pathname")
EVENT_FIELD_DESC = ("eventDescription", "description", "eventText", "text",
                    "message", "caption")
EVENT_FIELD_NUM = ("eventNumber", "orderLabel", "eventTag", "eventNum",
                   "number", "ordinal")
EVENT_FIELD_MAIN = ("main", "isMain", "primary")
DI_FIELDS = ("lastDefectInstanceId", "defectInstanceId")


# =========================================================== shared helpers
def cookie_value(cookie, name):
    for part in cookie.split(";"):
        k, _, v = part.strip().partition("=")
        if k == name:
            return v
    return None


def build_session(cookie):
    session = requests.Session()
    headers = {
        "Cookie": cookie,
        "Accept": "application/json, text/javascript, */*; q=0.01",
        "X-Requested-With": "XMLHttpRequest",
        "User-Agent": "coverity-export/1.1 (+python-requests)",
    }
    # Coverity uses double-submit CSRF: the XSRF-TOKEN cookie value must be
    # echoed in the X-XSRF-TOKEN header or POSTs are rejected (403).
    xsrf = cookie_value(cookie, "XSRF-TOKEN")
    if xsrf:
        headers["X-XSRF-TOKEN"] = xsrf
    session.headers.update(headers)
    return session


def resolve_cookie(args):
    if getattr(args, "cookie", None):
        return args.cookie
    if getattr(args, "cookie_file", None):
        return open(args.cookie_file, encoding="utf-8").read().strip()
    if os.environ.get("COVERITY_COOKIE"):
        return os.environ["COVERITY_COOKIE"]
    sys.exit("No cookie. Use --cookie/--cookie-file or set COVERITY_COOKIE.")


def is_auth_challenge(payload):
    return isinstance(payload, dict) and any(
        k in payload for k in AUTH_CHALLENGE_KEYS)


def http_json(session, url, *, params=None, json_body=None, method="GET",
              timeout=30, retries=4):
    """HTTP GET/POST returning parsed JSON; retries transient errors and aborts
    on an authentication challenge (expired session)."""
    last_err = None
    for attempt in range(1, retries + 1):
        try:
            if method == "POST":
                resp = session.post(url, params=params, json=json_body,
                                    timeout=timeout)
            else:
                resp = session.get(url, params=params, timeout=timeout)
            resp.raise_for_status()
            try:
                payload = resp.json()
            except ValueError:
                snippet = resp.text[:200].replace("\n", " ")
                raise RuntimeError(
                    "Response was not JSON (cookie expired?). "
                    f"First bytes: {snippet!r}")
            if is_auth_challenge(payload):
                raise SystemExit(
                    "Authentication failed: the server returned a login/SSO "
                    "config instead of data. Your COVERITY_COOKIE has expired "
                    "-- re-capture it from a logged-in browser session.")
            return payload
        except (requests.RequestException, RuntimeError) as err:
            last_err = err
            if attempt < retries:
                backoff = 2 ** (attempt - 1)
                print(f"    ! {err}; retry in {backoff}s [{attempt}/{retries}]",
                      file=sys.stderr)
                time.sleep(backoff)
    raise SystemExit(f"Giving up after {retries} attempts: {last_err}")


def _first(d, keys):
    for k in keys:
        if k in d and d[k] not in (None, ""):
            return d[k]
    return None


# ============================================================ defects (grid)
def _container(payload):
    if isinstance(payload, dict):
        inner = payload.get(RESULTSET_KEY)
        return inner if isinstance(inner, dict) else payload
    return {}


def _rows(payload):
    c = _container(payload)
    for k in ROW_KEYS:
        if isinstance(c.get(k), list):
            return c[k]
    return []


def _int_field(payload, keys, positive=False):
    c = _container(payload)
    for k in keys:
        v = c.get(k)
        if isinstance(v, int) and (v > 0 or not positive):
            return v
    return None


def _row_id(row):
    if isinstance(row, dict):
        for k in ID_KEYS:
            if k in row:
                return (k, row[k])
    return None


def _diagnose_empty(payload, raw_path):
    """A response came back but no rows were found: dump its shape and stop."""
    with open(raw_path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, ensure_ascii=False)

    def find_arrays(node, path=""):
        out = []
        if isinstance(node, dict):
            for k, v in node.items():
                p = f"{path}/{k}"
                if isinstance(v, list) and v:
                    out.append((p, len(v), type(v[0]).__name__))
                elif isinstance(v, (dict, list)):
                    out += find_arrays(v, p)
        elif isinstance(node, list):
            for i, v in enumerate(node[:1]):
                out += find_arrays(v, f"{path}[{i}]")
        return out

    tl = sorted(payload.keys()) if isinstance(payload, dict) else repr(type(payload))
    sys.exit(
        "\nGot a response but extracted 0 rows -- the /views/table.json reply "
        "uses a different layout than the defaults.\n"
        f"  top-level keys: {tl}\n"
        f"  non-empty arrays found: {find_arrays(payload)}\n"
        f"  raw first page saved to: {raw_path}\n"
        "Share that file (or these two lines) so the row/column mapping can be "
        "added.")


def collect_defects(session, views_url, reports_url, project_id, view_id, args):
    """Stateful paging: POST /views/table.json selects the page (server-side,
    keyed to the session cookie), then GET /reports/table.json returns that
    page's rows.  Both run on one session so the page state carries over."""
    pid, vid = int(project_id), int(view_id)

    def set_page(page):
        body = {"projectId": pid, "viewId": vid, args.page_param: page}
        http_json(session, views_url, json_body=body, method="POST",
                  timeout=args.timeout, retries=args.retries)

    def get_rows():
        return http_json(session, reports_url,
                         params={"projectId": pid, "viewId": vid},
                         method="GET", timeout=args.timeout, retries=args.retries)

    set_page(args.page_start)
    first = get_rows()
    if not _rows(first):
        _diagnose_empty(first, args.output + ".raw_page1.json")
    total = _int_field(first, TOTAL_KEYS, positive=True)
    page_size = _int_field(first, LIMIT_KEYS, positive=True) or len(_rows(first)) or 200
    npages = math.ceil(total / page_size) if total else None
    if total is not None:
        print(f"Server reports {total} defects; {page_size}/page -> "
              f"{npages} pages.", file=sys.stderr)

    all_rows, seen = [], set()
    k = 0
    while True:
        page = args.page_start + k
        if k == 0:
            data = first
        else:
            set_page(page)
            data = get_rows()
        rows = _rows(data)
        if not rows:
            break
        new = 0
        for row in rows:
            ident = _row_id(row)
            if ident is not None:
                if ident in seen:
                    continue
                seen.add(ident)
            all_rows.append(row)
            new += 1
        print(f"  page {page}: {len(rows)} rows ({new} new, {len(all_rows)} total)",
              file=sys.stderr)

        if total is not None and len(all_rows) >= total:
            break
        if npages is not None and (k + 1) >= npages:
            break
        if new == 0:
            sys.exit("\nPAGINATION PROBLEM: a full page added 0 new rows -- the "
                     "page-setter POST may not be taking effect for the data GET "
                     "(same session?). Re-capture the requests (README_export.md).")
        k += 1
        if args.delay:
            time.sleep(args.delay)
    return all_rows, total


def cmd_defects(args):
    session = build_session(resolve_cookie(args))
    print(f"Exporting defects (project {args.project_id}, view {args.view_id})...",
          file=sys.stderr)
    defects, total = collect_defects(session, args.views_url, args.reports_url,
                                     args.project_id, args.view_id, args)

    document = {
        "project_id": args.project_id, "view_id": args.view_id,
        "retrieved_at": datetime.now(timezone.utc).isoformat(),
        "server_total": total, "count": len(defects), "defects": defects,
    }
    with open(args.output, "w", encoding="utf-8") as fh:
        json.dump(document, fh, indent=2, ensure_ascii=False)
    note = ""
    if total is not None and len(defects) != total:
        note = f"  (WARNING: server reported {total}, collected {len(defects)})"
    print(f"\nWrote {len(defects)} defects to {args.output}{note}",
          file=sys.stderr)


# ================================================================== traces
def read_defect_index(source, id_field):
    data = json.loads(open(source, encoding="utf-8").read())
    grid_project_id = data.get("project_id") if isinstance(data, dict) else None
    if isinstance(data, dict) and isinstance(data.get("defects"), list):
        rows = data["defects"]
    elif (isinstance(data, dict) and isinstance(data.get("resultSet"), dict)
          and isinstance(data["resultSet"].get("results"), list)):
        rows = data["resultSet"]["results"]
    elif isinstance(data, list):
        rows = data
    else:
        rows = []
    cids, file_of, di_of, name_of = [], {}, {}, {}
    for row in rows:
        if not isinstance(row, dict):
            continue
        val = row.get(id_field, row.get("cid"))
        cid = int(val) if isinstance(val, int) or (
            isinstance(val, str) and val.isdigit()) else None
        if cid is None:
            continue
        if cid not in name_of:
            cids.append(cid)
        fiid = row.get("fileInstanceId")
        if fiid not in (None, "", -1, "-1"):
            file_of[cid] = str(fiid)
        for k in DI_FIELDS:
            di = row.get(k)
            if di not in (None, "", -1, "-1"):
                di_of[cid] = str(di)
                break
        name_of[cid] = row.get("displayFile") or row.get("filePath") or ""
    return {"cids": cids, "file_of": file_of, "di_of": di_of, "name_of": name_of,
            "project_id": grid_project_id}


def cids_from_string(text):
    out = []
    for tok in text.replace(",", " ").split():
        if not tok.isdigit():
            sys.exit(f"Invalid CID {tok!r}: expected digits only.")
        if int(tok) not in out:
            out.append(int(tok))
    return out


def is_detail_payload(payload):
    return isinstance(payload, dict) and (
        isinstance(payload.get("defects"), list)
        or isinstance(payload.get("events"), list))


def find_defect(payload, cid):
    for d in payload.get("defects", []):
        ids = d.get("ids", {}) if isinstance(d, dict) else {}
        if cid in (ids.get("cid"), ids.get("mergedDefectId"),
                   d.get("cid"), d.get("mergedDefectId")):
            return d
    return None


def normalize_event(ev, default_file, index, main_ids):
    return {
        "event": _first(ev, EVENT_FIELD_NUM)
        if _first(ev, EVENT_FIELD_NUM) is not None else index,
        "file": _first(ev, EVENT_FIELD_FILE) or default_file,
        "line": _first(ev, EVENT_FIELD_LINE),
        "description": _first(ev, EVENT_FIELD_DESC),
        "main": ev.get("id") in main_ids if ev.get("id") is not None
        else (bool(_first(ev, EVENT_FIELD_MAIN))
              if _first(ev, EVENT_FIELD_MAIN) is not None else None),
        "raw": ev,
    }


def extract_events_for_defect(payload, cid):
    defect = find_defect(payload, cid)
    if defect is None:
        return None
    di = str(defect.get("ids", {}).get("defectInstanceId", "")
             or defect.get("defectInstanceId", ""))
    default_file = payload.get("filePath") or payload.get("displayFile")
    main_ids = set()
    me = defect.get("mainEvent", {})
    if isinstance(me, dict) and me.get("id"):
        main_ids.add(me["id"])
    all_events = payload.get("events")
    if isinstance(all_events, list) and di:
        prefix = di + "-"
        chosen = [e for e in all_events if isinstance(e, dict)
                  and str(e.get("id", "")).startswith(prefix)]
    else:
        chosen = defect.get("events") if isinstance(
            defect.get("events"), list) else []
    return {"defect": defect,
            "events": [normalize_event(e, default_file, i + 1, main_ids)
                       for i, e in enumerate(chosen)]}


def fill_trace_url(template, project_id, fiid, di, cid):
    return (template.replace("{projectId}", str(project_id))
            .replace("{fileInstanceId}", str(fiid))
            .replace("{defectInstanceId}", str(di))
            .replace("{mergedDefectId}", str(cid))
            .replace("{cid}", str(cid)))


def load_valid_cache(raw_path, fiid, cid):
    if not os.path.exists(raw_path):
        return None
    try:
        cached = json.load(open(raw_path, encoding="utf-8"))
    except ValueError:
        return None
    if (is_detail_payload(cached)
            and str(cached.get("fileInstanceId", "")) == str(fiid)
            and find_defect(cached, cid) is not None):
        return cached
    return None


def cmd_traces(args):
    for tok in ("{fileInstanceId}", "{defectInstanceId}"):
        if tok not in args.url_template:
            sys.exit(f"--url-template must contain {tok}.")
    index = read_defect_index(args.cids_from, args.id_field)
    project_id = args.project_id or index.get("project_id")
    if not project_id:
        sys.exit("No projectId: pass --project-id, or use a --cids-from grid "
                 "export produced by the 'defects' subcommand (it records "
                 "project_id). The projectId is shown in the View Defects URL.")
    cids = cids_from_string(args.cids) if args.cids else index["cids"]
    if args.limit:
        cids = cids[:args.limit]
    miss = [c for c in cids
            if c not in index["file_of"] or c not in index["di_of"]]
    if miss:
        sys.exit(f"Missing fileInstanceId/defectInstanceId for CIDs {miss}; "
                 "pass the grid export from the 'defects' subcommand.")
    print(f"{len(cids)} CIDs to fetch.", file=sys.stderr)

    os.makedirs(args.output_dir, exist_ok=True)
    session = build_session(resolve_cookie(args))
    combined = []
    for n, cid in enumerate(cids, 1):
        fiid, di = index["file_of"][cid], index["di_of"][cid]
        raw_path = os.path.join(args.output_dir, f"cid_{cid}.json")
        payload = None if args.refresh else load_valid_cache(raw_path, fiid, cid)
        status = "cached"
        if payload is None:
            url = fill_trace_url(args.url_template, project_id, fiid, di, cid)
            payload = http_json(session, url, timeout=args.timeout,
                                retries=args.retries)
            got = str(payload.get("fileInstanceId", ""))
            if got and got != str(fiid):
                sys.exit(f"CID {cid}: asked for fileInstanceId {fiid} but got "
                         f"{got}. Check the {{fileInstanceId}} placeholder/host.")
            with open(raw_path, "w", encoding="utf-8") as fh:
                json.dump(payload, fh, ensure_ascii=False)
            status = "fetched"
            if args.delay:
                time.sleep(args.delay)
        res = extract_events_for_defect(payload, cid)
        if res is None:
            print(f"  [{n}/{len(cids)}] CID {cid}: NOT in payload "
                  f"({index['name_of'].get(cid,'?')}).", file=sys.stderr)
            combined.append({"cid": cid, "event_count": 0, "events": [],
                             "error": "cid-not-in-payload"})
            continue
        combined.append({"cid": cid, "file": payload.get("filePath"),
                         "event_count": len(res["events"]),
                         "events": res["events"]})
        short = (payload.get("filePath") or "").split("/lammps/")[-1]
        print(f"  [{n}/{len(cids)}] CID {cid}: {len(res['events'])} events "
              f"({status})  {short}", file=sys.stderr)

    document = {"retrieved_at": datetime.now(timezone.utc).isoformat(),
                "cid_count": len(combined), "traces": combined}
    with open(args.output, "w", encoding="utf-8") as fh:
        json.dump(document, fh, indent=2, ensure_ascii=False)
    empties = sum(1 for t in combined if t["event_count"] == 0)
    print(f"\nWrote {len(combined)} traces to {args.output} "
          f"(raw per-CID JSON in {args.output_dir}/).", file=sys.stderr)
    if empties:
        print(f"NOTE: {empties} CIDs produced 0 events (see 'error' fields).",
              file=sys.stderr)


# ================================================================ self-test
def cmd_self_test(_args):
    # --- defects: POST + pageNum paging/dedup ---
    class _R:
        def __init__(self, payload): self._p = payload
        def raise_for_status(self): pass
        def json(self): return self._p

    class FakeSession:                       # the stateful two-endpoint paging
        def __init__(self, total, page):
            self.total, self.page, self.cur = total, page, 1

        def post(self, url, params=None, json=None, timeout=None):
            self.cur = int(json["pageNum"])  # page-setter: stores page, returns no rows
            return _R({"table": {"pageNum": self.cur, "limit": self.page}})

        def get(self, url, params=None, timeout=None):
            off = (self.cur - 1) * self.page
            rows = [{"cid": i} for i in range(off, min(off + self.page, self.total))]
            return _R({"resultSet": {"offset": off, "limit": self.page,
                                     "totalCount": self.total, "results": rows}})

    class A:
        page_param, page_start, delay, timeout, retries = "pageNum", 1, 0, 5, 1
    defects, total = collect_defects(FakeSession(1813, 200), "views", "reports",
                                     16404, 70697, A())
    assert total == 1813 and len(defects) == 1813, "defects paging failed"
    assert [d["cid"] for d in defects] == list(range(1813)), "order/dedup failed"

    # --- traces: per-defect extraction + auth detection ---
    payload = {
        "filePath": "/x/src/GRANSURF/pair_surf_granular.cpp",
        "fileInstanceId": "193210750", "mainEventId": "777-1",
        "defects": [
            {"ids": {"cid": 504003, "defectInstanceId": 777}, "mainEvent": {"id": "777-1"}},
            {"ids": {"cid": 504036, "defectInstanceId": 888}, "mainEvent": {"id": "888-0"}}],
        "events": [
            {"id": "777-0", "lineNumber": 1615, "description": "decl pt"},
            {"id": "777-1", "lineNumber": 1700, "description": "use pt"},
            {"id": "888-0", "lineNumber": 1589, "description": "neg index"}],
    }
    r = extract_events_for_defect(payload, 504003)
    assert r and [e["line"] for e in r["events"]] == [1615, 1700], "trace filter failed"
    assert r["events"][1]["main"] is True, "main flag failed"
    assert extract_events_for_defect(payload, 999) is None, "absent cid"
    assert is_auth_challenge({"ldapConfigured": True}) and not is_auth_challenge(payload)
    assert fill_trace_url(DEFAULT_TRACE_URL, 16404, "F", "D", 503987).endswith(
        "fileInstanceId=F&defectInstanceId=D&mergedDefectId=503987&fileStart=&fileEnd=")
    assert cookie_value("a=1; XSRF-TOKEN=tok; b=2", "XSRF-TOKEN") == "tok", "xsrf parse"

    import tempfile
    dump = {"defects": [{"cid": 504003, "fileInstanceId": "193210750",
                         "lastDefectInstanceId": "777", "displayFile": "a.cpp"}]}
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as tf:
        json.dump(dump, tf)
        path = tf.name
    idx = read_defect_index(path, "cid")
    os.unlink(path)
    assert idx["file_of"][504003] == "193210750" and idx["di_of"][504003] == "777"
    print("Self-test passed: defects POST paging/dedup + traces extraction/auth/"
          "xsrf/index OK.")


# ===================================================================== main
def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[1],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    def auth(sp):
        sp.add_argument("--cookie", help="Raw Cookie header string.")
        sp.add_argument("--cookie-file", help="File with the Cookie string.")
        sp.add_argument("--delay", type=float, default=0.5,
                        help="Seconds between requests (default: 0.5).")
        sp.add_argument("--timeout", type=float, default=30,
                        help="Per-request timeout (default: 30).")
        sp.add_argument("--retries", type=int, default=4,
                        help="Retry attempts per request (default: 4).")

    d = sub.add_parser("defects", help="Export the defect table (grid).")
    d.add_argument("-p", "--project-id", default=12062,
                   help="projectId from the request body (e.g. 12062).")
    d.add_argument("-v", "--view-id", default=63757,
                   help="viewId from the request body (the view that lists all "
                        "defects, e.g. 63757).")
    d.add_argument("--reports-url", default=DEFAULT_REPORTS_URL,
                   help=f"Row-data endpoint, GET (default: {DEFAULT_REPORTS_URL}).")
    d.add_argument("--views-url", default=DEFAULT_VIEWS_URL,
                   help=f"Page-setter endpoint, POST (default: {DEFAULT_VIEWS_URL}).")
    d.add_argument("--page-param", default="pageNum",
                   help="Body field that selects the page (default: pageNum).")
    d.add_argument("--page-start", type=int, default=1,
                   help="First page number (default: 1; try 0 if rows are "
                        "missing).")
    d.add_argument("-o", "--output", default="exported_defects.json")
    auth(d)
    d.set_defaults(func=cmd_defects)

    t = sub.add_parser("traces", help="Export per-defect event traces.")
    t.add_argument("--cids-from", required=True,
                   help="Grid export from 'defects' (rows carry cid, "
                        "fileInstanceId, lastDefectInstanceId).")
    t.add_argument("--cids", help="Optional comma/space-separated subset of CIDs.")
    t.add_argument("--url-template", default=DEFAULT_TRACE_URL,
                   help="source.json URL; needs {fileInstanceId} and "
                        "{defectInstanceId} (also {projectId}{mergedDefectId}).")
    t.add_argument("--project-id", default=None,
                   help="projectId for the source.json URL. Defaults to the "
                        "'project_id' recorded in the --cids-from grid export; "
                        "pass explicitly to override or if the grid lacks it.")
    t.add_argument("--id-field", default="cid")
    t.add_argument("--output-dir", default="defect_traces")
    t.add_argument("-o", "--output", default="traces.json")
    t.add_argument("--limit", type=int, help="Process only the first N CIDs.")
    t.add_argument("--refresh", action="store_true",
                   help="Re-fetch even if a valid raw file already exists.")
    auth(t)
    t.set_defaults(func=cmd_traces)

    st = sub.add_parser("self-test", help="Offline logic checks; no network.")
    st.set_defaults(func=cmd_self_test)

    args = p.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
