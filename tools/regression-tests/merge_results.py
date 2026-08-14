#!/usr/bin/env python3
'''
Merge JUnit XML files from distributed regression test runs (see run_tests.py)
and generate machine readable and human readable summaries.

Only the Python standard library is required.

Input:
    + one or more JUnit XML files, e.g. output-0.xml output-1.xml ...
      (files that are empty or cannot be parsed are skipped with a warning)

Output (all optional, selected via command line flags):
    + --output-xml merged.xml : all test cases merged into a single JUnit XML file
    + --json run.json         : test results and metadata as a JSON file; this is
                                the input format for --compare and is meant to be
                                archived for tracking results over time
    + --summary summary.md    : summary in Markdown format (use "-" for stdout),
                                suitable for $GITHUB_STEP_SUMMARY
    + --cost-file costs.json  : measured run times per input script, to be passed to
                                the next "run_tests.py --analyze" run with --cost-file
                                so that it can distribute the input scripts evenly
                                over the workers; an existing file is updated in
                                place, keeping the times of tests that did not run

With --compare previous-run.json, the Markdown summary also reports the changes
relative to a previous run: new failures, fixed tests, new tests, removed tests.

The Markdown summary also groups the failed runs by where their output starts to
deviate from the reference log file: runs that already differ at the start of the run
point at a difference in the model or in the setup, while runs that only differ towards
the end of the run look like chaotic amplification of rounding differences.

Since the cost of a test depends on how it is run, each test configuration (serial,
4 MPI tasks, ...) needs its own cost file.

The summary also lists the input scripts that a developer has to fix before their result
means anything: input scripts whose initial velocities depend on the number of MPI
processes, and production sized runs that should be shortened.  REPORTING.md in this
folder documents all the fields of the JSON output and how to classify them.

Example:
    python3 merge_results.py --title "Full Regression Test" \
            --output-xml regression.xml --json run.json --summary summary.md \
            --cost-file costs.json --compare previous-run.json output-*.xml
'''

from argparse import ArgumentParser
import datetime
import json
import os
import sys
import xml.etree.ElementTree as ET

# statuses that count as broken for the comparison with a previous run
BAD = ('failed', 'error')

# A classical MD trajectory is chaotic, so the tiniest rounding difference from another
# machine, compiler or number of MPI processes grows exponentially and shows up in the
# thermo output after roughly this many steps.  A deviation that appears later than this
# is expected and says nothing about the correctness of the code, while one that appears
# within the first steps is a difference in what is computed.
CHAOS_STEPS = 1000
SUSPICIOUS_STEPS = 200

'''
    parse a single JUnit XML file as written by run_tests.py (a single <testsuite>
    root element, or a <testsuites> wrapper) and return a tuple of
    (properties dict, tests dict keyed by "classname/name")
'''
def parse_junit_xml(filename):
    tree = ET.parse(filename)
    root = tree.getroot()
    if root.tag == 'testsuites':
        suites = root.findall('testsuite')
    elif root.tag == 'testsuite':
        suites = [root]
    else:
        raise ET.ParseError(f"unexpected root element <{root.tag}>")

    properties = {}
    tests = {}
    for suite in suites:
        for prop in suite.findall('./properties/property'):
            properties.setdefault(prop.get('name', ''), prop.get('value', ''))
        for case in suite.findall('testcase'):
            key = case.get('classname', '') + '/' + case.get('name', '')
            entry = {'status': 'passed', 'time': float(case.get('time', 0.0)), 'message': ''}
            for tag in ('failure', 'error', 'skipped'):
                elem = case.find(tag)
                if elem is not None:
                    entry['status'] = 'failed' if tag == 'failure' else tag
                    entry['message'] = elem.get('message', '')
                    if elem.text:
                        entry['details'] = elem.text.strip()
                    break
            # how far into the run the output starts to deviate from the reference
            if case.get('diverged'):
                entry['diverged'] = int(case.get('diverged'))
                entry['diverged_row'] = int(case.get('diverged-row', 0))
                if case.get('diverged-at') is not None:
                    entry['diverged_at'] = int(case.get('diverged-at'))
            if case.get('attention'):
                entry['attention'] = case.get('attention')
            tests[key] = entry
    return properties, tests

'''
    write the merged test data back out as a single JUnit XML file
'''
def write_merged_xml(filename, title, properties, tests, counts):
    testsuite = ET.Element('testsuite')
    testsuite.set('name', title)
    testsuite.set('timestamp', datetime.datetime.now().isoformat(timespec='seconds'))
    testsuite.set('tests', str(counts['tests']))
    testsuite.set('failures', str(counts['failed']))
    testsuite.set('errors', str(counts['error']))
    testsuite.set('skipped', str(counts['skipped']))
    testsuite.set('time', f"{counts['time']:.3f}")

    if properties:
        props = ET.SubElement(testsuite, 'properties')
        for key, value in properties.items():
            prop = ET.SubElement(props, 'property')
            prop.set('name', str(key))
            prop.set('value', str(value))

    for key in sorted(tests):
        entry = tests[key]
        classname, _, name = key.rpartition('/')
        case = ET.SubElement(testsuite, 'testcase')
        case.set('name', name)
        case.set('classname', classname)
        case.set('time', f"{entry['time']:.3f}")
        if 'diverged' in entry:
            case.set('diverged', str(entry['diverged']))
            case.set('diverged-row', str(entry['diverged_row']))
            if 'diverged_at' in entry:
                case.set('diverged-at', str(entry['diverged_at']))
        if 'attention' in entry:
            case.set('attention', entry['attention'])
        if entry['status'] != 'passed':
            tag = 'failure' if entry['status'] == 'failed' else entry['status']
            elem = ET.SubElement(case, tag)
            elem.set('message', entry['message'])
            if entry.get('details'):
                elem.text = entry['details']

    tree = ET.ElementTree(testsuite)
    ET.indent(tree)
    tree.write(filename, encoding='UTF-8', xml_declaration=True)

'''
    update the file with the measured run times of the individual input scripts

    The file maps the test names ("<folder>/<input script>", relative to the top-level
    examples folder) to the measured wall-clock time of the run in seconds.  It is read
    back by "run_tests.py --analyze --cost-file ..." to distribute the input scripts
    over the workers such that they all take about the same time.

    Tests that did not run in this test run have a recorded time of zero; their previous
    time is kept, since it is a better estimate than assuming that they are free.
'''
def write_cost_file(filename, title, tests):
    costs = {}
    if os.path.isfile(filename):
        try:
            with open(filename, 'r') as f:
                costs = json.load(f).get('costs', {})
        except (OSError, ValueError) as err:
            print(f"WARNING: ignoring unreadable cost file {filename}: {err}", file=sys.stderr)

    updated = 0
    for key, entry in tests.items():
        time = round(float(entry['time']), 3)
        if time > 0.0 or key not in costs:
            costs[key] = time
            updated += 1

    data = {
        'metadata': {
            'title': title,
            'generated': datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
            'updated': updated,
        },
        'costs': dict(sorted(costs.items())),
    }
    with open(filename, 'w') as f:
        json.dump(data, f, indent=1)
        f.write('\n')

'''
    format a capped Markdown list of test names with messages
'''
def markdown_list(keys, tests, maxlen=50):
    lines = []
    for key in sorted(keys)[:maxlen]:
        message = tests[key]['message'] if key in tests else ''
        if message:
            lines.append(f"- `{key}`: {message}")
        else:
            lines.append(f"- `{key}`")
    if len(keys) > maxlen:
        lines.append(f"- ... and {len(keys) - maxlen} more")
    return '\n'.join(lines) + '\n'

if __name__ == "__main__":
    parser = ArgumentParser(description="Merge and summarize JUnit XML regression test results")
    parser.add_argument("xmlfiles", nargs='+', help="JUnit XML files to merge")
    parser.add_argument("--title", default="Regression Test", help="Title used in the merged output")
    parser.add_argument("--output-xml", dest="output_xml", default="", help="Merged JUnit XML output file")
    parser.add_argument("--json", dest="json_file", default="", help="JSON output file with results and metadata")
    parser.add_argument("--summary", dest="summary", default="", help="Markdown summary output file ('-' for stdout)")
    parser.add_argument("--compare", dest="compare", default="", help="JSON file from a previous run to compare against")
    parser.add_argument("--cost-file", dest="cost_file", default="",
                        help="JSON file with the measured run times, updated in place")
    parser.add_argument("--max-details", dest="max_details", type=int, default=0,
                        help="Number of characters of the detailed report of each failed "
                             "test to keep in the JSON output (0: none, the reports are "
                             "always kept in the merged JUnit XML file)")
    parser.add_argument("--commit", dest="commit", default="", help="Commit hash for git revision that was tested")
    parser.add_argument("--branch", dest="branch", default="", help="Git branch of revision that was tested")
    args = parser.parse_args()

    # merge the shards
    properties = {}
    tests = {}
    for xmlfile in args.xmlfiles:
        try:
            shard_properties, shard_tests = parse_junit_xml(xmlfile)
        except ET.ParseError as err:
            print(f"WARNING: skipping {xmlfile}: {err}", file=sys.stderr)
            continue
        for key, value in shard_properties.items():
            properties.setdefault(key, value)
        duplicates = set(shard_tests) & set(tests)
        if duplicates:
            print(f"WARNING: {len(duplicates)} duplicate test(s) in {xmlfile}, keeping the first result",
                  file=sys.stderr)
            for key in duplicates:
                del shard_tests[key]
        tests.update(shard_tests)

    counts = {'tests': len(tests), 'passed': 0, 'failed': 0, 'error': 0, 'skipped': 0, 'time': 0.0}
    for entry in tests.values():
        counts[entry['status']] += 1
        counts['time'] += entry['time']

    if args.output_xml:
        write_merged_xml(args.output_xml, args.title, properties, tests, counts)

    if args.cost_file:
        write_cost_file(args.cost_file, args.title, tests)

    if args.json_file:
        # the detailed reports of all failed tests make the archived JSON files several
        # times larger, so they are only included when asked for and then truncated
        entries = {}
        for key, entry in tests.items():
            entries[key] = { name: value for name, value in entry.items() if name != 'details' }
            if args.max_details > 0 and 'details' in entry:
                entries[key]['details'] = entry['details'][:args.max_details]
        data = {
            'metadata': {
                'title': args.title,
                'generated': datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
                'properties': properties,
                'counts': counts,
            },
            'tests': entries,
        }
        if args.commit:
            data['metadata']['commit'] = args.commit
        if args.branch:
            data['metadata']['branch'] = args.branch
        with open(args.json_file, 'w') as f:
            json.dump(data, f, indent=2)
            f.write('\n')

    # generate the Markdown summary
    if args.summary:
        md = f"## {args.title}\n\n"
        md += "| Tests | Passed | Failed | Errors | Skipped | Walltime |\n"
        md += "|------:|-------:|-------:|-------:|--------:|---------:|\n"
        md += (f"| {counts['tests']} | {counts['passed']} | {counts['failed']} | {counts['error']} |"
               f" {counts['skipped']} | {counts['time']:.0f} s |\n\n")

        slowest = sorted(tests, key=lambda key: -tests[key]['time'])[:10]
        if slowest and tests[slowest[0]]['time'] > 0.0:
            md += f"The {len(slowest)} most expensive runs, which set the lower bound for the"
            md += " walltime of the whole test run:\n\n"
            md += ''.join(f"- `{key}`: {tests[key]['time']:.1f} s\n" for key in slowest) + '\n'

        failed = [key for key, entry in tests.items() if entry['status'] == 'failed']
        errors = [key for key, entry in tests.items() if entry['status'] == 'error']

        # inputs that a developer has to fix before their result means anything
        attention = [key for key, entry in tests.items() if entry.get('attention')]
        if attention:
            md += f"### Input scripts that need to be fixed ({len(attention)})\n\n"
            md += ("These have to be corrected in the repository and their reference log"
                   " files recreated; until then their result says little.\n\n")
            md += ''.join(f"- `{key}`: {tests[key]['attention']}\n" for key in sorted(attention)[:50])
            if len(attention) > 50:
                md += f"- ... and {len(attention) - 50} more\n"
            md += '\n'

        # Sort the failed runs by where their output starts to deviate from the reference
        # log file: that is what says whether a failure is worth investigating.  Inputs
        # that are known not to be able to match their reference log file are left out,
        # since their early deviation is already explained and would bury the rest.
        early = [key for key in failed
                 if not tests[key].get('attention')
                 and (tests[key].get('diverged_row') == 0
                      or tests[key].get('diverged_at', CHAOS_STEPS + 1) <= SUSPICIOUS_STEPS)]
        late = [key for key in failed if tests[key].get('diverged_at', 0) > CHAOS_STEPS]
        if early or late:
            md += "### Where the failed runs start to deviate\n\n"
            if early:
                md += (f"**Deviating within the first {SUSPICIOUS_STEPS} steps ({len(early)}):**"
                       " too early for a rounding difference to have grown that far, so these"
                       " are a difference in what is computed until proven otherwise.\n\n")
                md += markdown_list(early, tests) + '\n'
            if late:
                md += (f"**Deviating only after {CHAOS_STEPS} steps ({len(late)}):** a classical"
                       " MD trajectory is chaotic, so a rounding difference from another"
                       " machine, compiler or number of MPI processes grows into a visible"
                       " deviation on its own. Expect these unless something else is wrong;"
                       " a shorter run or a wider tolerance is the usual answer.\n\n")
                md += markdown_list(late, tests) + '\n'

        if failed:
            md += f"### Failed numerical checks ({len(failed)})\n\n"
            md += markdown_list(failed, tests) + '\n'
        if errors:
            md += f"### Runs with errors ({len(errors)})\n\n"
            md += markdown_list(errors, tests) + '\n'

        if args.compare:
            with open(args.compare, 'r') as f:
                previous = json.load(f)
            tests_prev = previous.get('tests', {})
            generated_prev = previous.get('metadata', {}).get('generated', 'unknown date')

            new_failures = [k for k in tests if k in tests_prev
                            and (tests[k]['status'] in BAD) and (tests_prev[k]['status'] not in BAD)]
            fixed = [k for k in tests if k in tests_prev
                     and (tests[k]['status'] == 'passed') and (tests_prev[k]['status'] in BAD)]
            new_tests = [k for k in tests if k not in tests_prev]
            removed = [k for k in tests_prev if k not in tests]

            md += f"### Changes relative to the previous run (from {generated_prev})\n\n"
            if not (new_failures or fixed or new_tests or removed):
                md += "No changes.\n\n"
            if new_failures:
                md += f"**New failures ({len(new_failures)}):**\n\n"
                md += markdown_list(new_failures, tests) + '\n'
            if fixed:
                md += f"**Fixed since the previous run ({len(fixed)}):**\n\n"
                md += markdown_list(fixed, {}) + '\n'
            if new_tests:
                md += f"**New tests ({len(new_tests)}):**\n\n"
                md += markdown_list(new_tests, tests) + '\n'
            if removed:
                md += f"**Removed tests ({len(removed)}):**\n\n"
                md += markdown_list(removed, {}) + '\n'

        if args.summary == '-':
            sys.stdout.write(md)
        else:
            with open(args.summary, 'w') as f:
                f.write(md)
