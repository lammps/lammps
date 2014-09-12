#!/bin/sh
# do some cleaning up to fix permissions and remove cruft

for f in `find ./ -name \*.orig -or -name \*~ -or -name \*.rej -print`
do rm -v $f
done

for f in `find ./ -name \*.cpp -or -name \*.c -or -name \*.h -print`
do chmod -x $f
done

