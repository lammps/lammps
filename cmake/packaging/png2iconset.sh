#!/bin/sh

if [ $# != 2 ]
then
   echo "usage: $0 <pngfile> <iconset name>"
   exit 1
fi

png="$1"
ico="$2"

if [ ! -f ${png} ]
then
   echo "PNG Image $1 not found"
fi

rm -rf ${ico}.iconset
mkdir ${ico}.iconset
sips -z   16   16 ${png} --out ${ico}.iconset/icon_16x16.png
sips -z   32   32 ${png} --out ${ico}.iconset/icon_16x16@2x.png
sips -z   32   32 ${png} --out ${ico}.iconset/icon_32x32.png
sips -z   64   64 ${png} --out ${ico}.iconset/icon_32x32@2x.png
sips -z  128  128 ${png} --out ${ico}.iconset/icon_128x128.png
sips -z  256  256 ${png} --out ${ico}.iconset/icon_128x128@2x.png
sips -z  256  256 ${png} --out ${ico}.iconset/icon_256x256.png
sips -z  512  512 ${png} --out ${ico}.iconset/icon_256x256@2x.png
sips -z  512  512 ${png} --out ${ico}.iconset/icon_512x512.png
sips -z 1024 1024 ${png} --out ${ico}.iconset/icon_512x512@2x.png
iconutil -c icns ${ico}.iconset
rm -rf ${ico}.iconset
