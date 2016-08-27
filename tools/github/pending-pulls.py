#!/usr/bin/env python
import pycurl
import json
from io import BytesIO

# settings
base = 'akohlmey/lammps'
user = 'sjplimp'
upstream = 'integration'
verbose = True

buf = BytesIO()
c = pycurl.Curl()
c.setopt(c.URL, 'https://api.github.com/repos/'+base+'/pulls?state=open')
c.setopt(c.WRITEFUNCTION, buf.write)
c.perform()

result = json.loads(buf.getvalue().decode());

print('Open pull requests for repository: '+base+'\n')

for pull in result:
    if pull['assignee'] and pull['assignee']['login'] == user:
        num = pull['number']
        print('Pending pull request #%d' % num)
        print('Assigned to %s' % pull['assignee']['login'])
        print('Submitted by: %s' % pull['head']['repo']['owner']['login'])
        print('Title: '+pull['title'])
        print('URL: https://github.com/'+base+'/pull/'+str(num))
        # get pull repository and branch
        repo = pull['head']['repo']['clone_url']
        branch = pull['head']['ref']
        # instructions
        print('\nCommand line instructions, step 1:\n')
        print('git checkout -b merge-pull-%d %s' % (num,upstream))
        print('git pull %s %s' % (repo, branch))
        print('\nCommand line instructions, step 2:\n')
        print('git checkout '+upstream)
        print('git merge --no-ff merge-pull-%d' % num)
        if verbose:
          print('\nDescription:\n%s' % pull['body'])
        print('------------------\n')
