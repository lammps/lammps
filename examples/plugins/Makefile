CXX=mpicxx
CXXFLAGS=-I../../src -Wall -O3 -fPIC -I../../src/USER-OMP -fopenmp
LD=$(CXX) -shared -rdynamic -fopenmp

morse2plugin.so: morse2plugin.o pair_morse2.o pair_morse2_omp.o
	$(LD) -o $@ $^

.cpp.o:
	$(CXX) -o $@ $(CXXFLAGS) -c $<

pair_morse2.o: pair_morse2.cpp pair_morse2.h
pair_morse2_omp.o: pair_morse2_omp.cpp pair_morse2_omp.h pair_morse2.h
morse2plugin.o: morse2plugin.cpp pair_morse2.h pair_morse2_omp.h

clean:
	rm -f *~ *.so *.o log.lammps


