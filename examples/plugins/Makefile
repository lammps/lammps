CXX=mpicxx
CXXFLAGS=-I../../src -Wall -Wextra -O3 -fPIC -I../../src/USER-OMP -fopenmp
LD=$(CXX) -shared -rdynamic -fopenmp

default: morse2plugin.so nve2plugin.so helloplugin.so

helloplugin.so: helloplugin.o
	$(LD) -o $@ $^

morse2plugin.so: morse2plugin.o pair_morse2.o pair_morse2_omp.o
	$(LD) -o $@ $^

nve2plugin.so: nve2plugin.o fix_nve2.o
	$(LD) -o $@ $^

.cpp.o:
	$(CXX) -o $@ $(CXXFLAGS) -c $<

helloplugin.o: helloplugin.cpp

pair_morse2.o: pair_morse2.cpp pair_morse2.h
pair_morse2_omp.o: pair_morse2_omp.cpp pair_morse2_omp.h pair_morse2.h
morse2plugin.o: morse2plugin.cpp pair_morse2.h pair_morse2_omp.h

fix_nve2.o: fix_nve2.cpp fix_nve2.h
nve2plugin.o: nve2plugin.cpp fix_nve2.h

clean:
	rm -f *~ *.so *.o log.lammps

