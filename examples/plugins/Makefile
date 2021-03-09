CXX=mpicxx
CXXFLAGS=-I../../src -Wall -O3 -fPIC
LD=$(CXX) -shared -rdynamic

morse2plugin.so: morse2plugin.o pair_morse2.o
	$(LD) -o $@ $^

.cpp.o:
	$(CXX) -o $@ $(CXXFLAGS) -c $<

pair_morse2.o: pair_morse2.cpp pair_morse2.h
morse2plugin.o: morse2plugin.cpp pair_morse2.h

clean:
	rm -f *~ *.so *.o log.lammps


