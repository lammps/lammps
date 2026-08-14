# Exporting Coverity Scan defects and traces (`coverity_export.py`)

`coverity_export.py` pulls LAMMPS's static-analysis results off the **free /
open-source** [Coverity Scan](https://scan.coverity.com/) service (project
`LAMMPS`, `projectId=16404`) into local JSON for offline triage. It drives the
same *undocumented* internal endpoints the "View Defects" web UI calls, reusing
your authenticated browser session, so it can break if Coverity changes those
endpoints. It is for periodic personal exports, not automation.

(For the separate Coverity *modeling* file that suppresses false positives, see
`README_modeling.md`.)

## Requirements

```bash
pip install requests
```

## Two steps

| Subcommand | Request | Output |
|---|---|---|
| `defects` | per page: `POST /views/table.json` (sets the page) then `GET /reports/table.json` (returns the rows) | `exported_defects.json` — every defect row (the "grid") |
| `traces`  | `GET /sourcebrowser/source.json?fileInstanceId=..&defectInstanceId=..` | `traces.json` + `defect_traces/cid_*.json` — per-defect event traces (file + line) |

`traces` reads the ids it needs (`fileInstanceId`, `lastDefectInstanceId`) out of
the `defects` output, so always run `defects` first.

## Getting the authentication cookie (do this first)

The endpoints require your logged-in session cookie. It is **not** an API token;
it is the browser session, and it **expires** (the tool aborts with
"Authentication failed … re-capture it" when it does — just repeat these steps).

1. Log in to <https://scan.coverity.com/> in your browser and open the LAMMPS
   project's **View Defects** page.
2. Open the browser **Developer Tools → Network** tab (F12), and reload the
   defects view so requests appear.
3. Find a request to the defect data — `table.json` (the grid) or `source.json`
   (one defect). Right-click it → **Copy → Copy as cURL**.
4. From the copied cURL, take the value of the `-b '…'` / `-H 'cookie: …'`
   header. On the open-source tier it looks like
   `COVJSESSIONID-build=…; XSRF-TOKEN=…; isAuthenticated=true; …`. Export it:

   ```bash
   export COVERITY_COOKIE='COVJSESSIONID-build=...; XSRF-TOKEN=...; isAuthenticated=true; ...'
   ```

   Or put it in a file and pass `--cookie-file FILE`, or inline with `--cookie`.

> **The XSRF token is handled for you.** The POST endpoint uses double-submit
> CSRF protection: the `XSRF-TOKEN` cookie value must be echoed in an
> `X-XSRF-TOKEN` header. The script extracts it from your cookie automatically,
> so just include the whole cookie string (with `XSRF-TOKEN=…` in it).

> **Do not commit the cookie.** It is a live credential. Keep it in
> `COVERITY_COOKIE` or a file outside the repo. The exported JSON is also best
> kept out of the repo.

## Step 1 — export the defect grid (all pages)

You need the project id and the **view id** of the view that lists every defect
(from the `views/table.json` request body, e.g.
`{"projectId":16404,"viewId":70697,"pageNum":1}`). Then:

```bash
python coverity_export.py defects -p 16404 -v 70697 -o exported_defects.json
```

A single command collects the **whole** set. The output records `server_total`
and `count`; a mismatch is flagged. (If your view is 0-based and the first 200
rows are missing, add `--page-start 0`.)

> **How paging works here (it is stateful).** On the open-source tier the page
> is *not* chosen by the data request's URL. For each page the tool makes two
> requests on one session: `POST /views/table.json {projectId,viewId,pageNum:N}`
> to **set** the page server-side, then `GET /reports/table.json?projectId&viewId`
> to **read** that page's rows. That is why scraping the data GET alone (with an
> `offset`/`pageNum` in its URL) only ever returned the first page — the page
> selector is the side-channel POST, and the two must share one session.

## Step 2 — export per-defect traces

Pick the CIDs you want the event trace for (the numbered data-flow steps with
file + line that pinpoint each defect). The ids are taken from the grid export,
and the `source.json` URL is built in, so this needs only the cookie:

```bash
python coverity_export.py traces \
    --cids-from exported_defects.json \
    --cids '504003 504028 503971' \
    -o traces.json
```

Omit `--cids` to fetch traces for every CID in the grid export. Each
`defect_traces/cid_<CID>.json` is the raw `source.json` payload (one source file
plus that defect's events); `traces.json` is the flattened per-CID event list.
Re-runs reuse cached files; pass `--refresh` to force a re-fetch (and to discard
files left over from a failed/auth-expired run).

### Verifying a traces run

A good run has, for each `cid_<CID>.json`, a `filePath` for the expected source
file and a `defects[]` entry whose `ids.cid == CID`. If instead a file contains
SSO/login keys (`availableSamlSsoConfigurations`, `ldapConfigured`) the cookie
expired; if every file is the *same* unrelated source file, the `source.json`
URL is not selecting by `fileInstanceId`.

## Offline check

```bash
python coverity_export.py self-test
```

exercises the POST paging/dedup and trace-extraction/auth-detection/XSRF logic
with no network or cookie.
