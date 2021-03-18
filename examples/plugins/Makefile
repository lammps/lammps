CXX=mpicxx
CXXFLAGS=-I../../src -Wall -Wextra -O3 -fPIC -I../../src/USER-OMP -fopenmp
LD=$(CXX) -shared -rdynamic -fopenmp
DSOEXT=.so

include Makefile.common
