SHELL=/bin/sh
CXX=g++
CXXFLAGS=-Wall -fno-rtti -fno-exceptions -O2 -g -I../STUBS -I. -I..
LIBS=-Wl,--start-group ../liblmp_serial.a ../../lib/meam/libmeam.a ../../lib/poems/libpoems.a \
	../../lib/colvars/libcolvars.a ../STUBS/libmpi.a -Wl,--end-group -lgfortran -ldl

default: molfile_tester.x

clean:
	rm -vf *.o *.x *~

molfile_tester.x: molfile_tester.o molfile_interface.o dump_molfile.o
	$(CXX) -o $@ $^ $(LIBS)

molfile_tester.o: molfile_tester.cpp molfile_interface.h dump_molfile.h
	$(CXX) $(CXXFLAGS) -o $@ -c $<

molfile_interface.o: molfile_interface.cpp molfile_interface.h
	$(CXX) $(CXXFLAGS) -o $@ -c $<

dump_molfile.o: dump_molfile.cpp dump_molfile.h molfile_interface.h
	$(CXX) $(CXXFLAGS) -o $@ -c $<


