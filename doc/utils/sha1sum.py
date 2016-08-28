#!/bin/env python
# simple utility which reimplements sha1sum using Python
import hashlib
import sys
s = hashlib.sha1()
s.update(sys.stdin.read().encode())
print(s.hexdigest())
