#!/bin/tcsh

set files = `echo *.h *.cpp *.cu`
rm -rf /tmp/cpp5678
mkdir /tmp/cpp5678
mkdir /tmp/cpp5678/lgpu

foreach file ( $files )
#	/bin/cp $file /tmp/cpp5678/$file:t:t
	# ------ Sed Replace
	sed -i 's/atom->dev_engv/ans->dev_engv/g' $file
end

