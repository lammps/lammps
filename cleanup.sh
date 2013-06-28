#!/bin/sh
# do some cleaning up to fix permissions and remove cruft

find ./ -name \*.orig -print | xargs rm -v
find ./ -name \*~ -print | xargs rm -v

find ./ -name \*.cpp -print | xargs chmod -x
find ./ -name \*.h -print | xargs chmod -x

