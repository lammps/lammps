#!/usr/bin/env python
import pycurl
import json
from io import BytesIO

#############################################################
# settings:
# github repository from pull requests are culled
base = 'akohlmey/lammps'
# github userid of pull request assignee
user = 'sjplimp'
# upstream branch into which pull requests are merged
upstream = 'integration'
# whether pull request descriptions should be printed, too.
verbose = True
#############################################################

buf = BytesIO()
c = pycurl.Curl()
c.setopt(c.URL,'https://api.github.com/repos/'+base+'/pulls?state=open,assignee='+user)
c.setopt(c.WRITEFUNCTION,buf.write)
c.perform()

result = json.loads(buf.getvalue().decode(encoding="utf-8"));

print('\nPending pull requests for repository '+base+' assigned to '+user+'\n')

for pull in result:
    if pull['assignee']:
        num = pull['number']
        print('Pending pull request #%d' % num)
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
