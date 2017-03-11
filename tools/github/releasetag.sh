#!/bin/bash

ver=$(sed -e 's,^.*"\(.*\)".*$,\1,' < version.h)
tag=$(date +"patch_%Y%m%d" --date="${ver}")
echo git tag ${tag}
echo git push --tags origin

