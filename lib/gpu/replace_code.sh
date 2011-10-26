#!/bin/tcsh

set files = `echo *.h *.cu`
rm -rf /tmp/cpp5678
mkdir /tmp/cpp5678
mkdir /tmp/cpp5678/lgpu

foreach file ( $files )
#	/bin/cp $file /tmp/cpp5678/$file:t:t
	# ------ Sed Replace
	sed -i.bak 's/(numtyp)1\.0\/rsq/ucl_recip(rsq)/g' $file
end

