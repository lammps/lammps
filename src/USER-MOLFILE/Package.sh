# Update package files in LAMMPS
# copy package file to src if it doesn't exists or is different

for file in molfile_interface.cpp molfile_interface.h molfile_plugin.h \
    dump_molfile.cpp dump_molfile.h reader_molfile.h reader_molfile.cpp \
    vmdplugin.h ; do \
  if (test ! -e ../$file) then
    echo "  creating src/$file"
    cp $file ..
  elif ! cmp -s $file ../$file ; then
    echo "  updating src/$file"
    cp $file ..
  fi
done

