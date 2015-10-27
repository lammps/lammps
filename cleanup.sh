#!/bin/sh
# do some cleaning up to fix permissions and remove cruft

for f in `find ./ -name \*.orig -print -or -name \*~ -print -or -name \*.rej -print`
do rm -v $f
done

for f in `find ./ -name \*.cpp -print -or -name \*.c -print -or -name \*.h -print`
do chmod -x $f
done

