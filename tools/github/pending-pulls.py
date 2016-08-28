#!/usr/bin/env python
import pycurl
import json
import sys
from io import BytesIO
from datetime import datetime

#############################################################
# settings:
# github repository from which pull requests are culled
base = 'akohlmey/lammps'
# github userid of pull request assignee
user = 'sjplimp'
# upstream branch into which pull requests are merged
upstream = 'integration'
# whether pull request descriptions should be printed, too.
verbose = True
#############################################################

# for writing separate python 2/3 code blocks to work
# around python's unicode vs. bytes vs. strings insanity.
if (sys.version_info > (3, 0)):
    pyver = 3
else:
    pyver = 2

# set up a query to the github web API to return information about all
# open pull requests for a given repository which are assigned to the
# given userid and store the response in a buffer.
buf = BytesIO()
c = pycurl.Curl()
c.setopt(c.URL,'https://api.github.com/repos/'+base+'/pulls?state=open,assignee='+user)
c.setopt(c.WRITEFUNCTION,buf.write)
c.perform()

# github responds with a JSON format text listing all pull requests.
# we convert the response into a dictionary with individual entries
# being dictionaries as well. The JSON data is in UTF-8 encoding.
result = json.loads(buf.getvalue().decode(encoding="utf-8"));

print('\n%d pending pull requests for repository %s assigned to %s\n' % (len(result),base,user))

# loop over entrees
for pull in result:
    # get pull request id and repo and branch info
    num = pull['number']
    repo = pull['head']['repo']['clone_url']
    branch = pull['head']['ref']
    # print some general info
    print(str('Pending pull request #%d' % num))
    print('Submitted by: %s' % pull['head']['user']['login'])
    print('Last update: %s' % datetime.strptime(pull['updated_at'],'%Y-%m-%dT%H:%M:%SZ').ctime())
    if pyver == 2:
        print('Title: %s' % pull['title'].encode('utf-8'))
    else:
        print('Title: %s' % pull['title'])
    print('URL: https://github.com/'+base+'/pull/'+str(num))
    # print command line suggestions
    print('\nCommand line instructions, step 1:\n')
    print('git checkout -b merge-pull-%d %s' % (num, upstream))
    print('git pull %s %s' % (repo, branch))
    print('\nCommand line instructions, step 2:\n')
    print('git checkout '+upstream)
    print('git merge --no-ff merge-pull-%d' % num)
    # if requested, display pull request description
    if verbose:
        if pyver == 2:
            print('\nDescription:\n%s' % pull['body'].encode('utf-8'))
        else:
            print('\nDescription:\n%s' % pull['body'])
        print('------------------\n')
