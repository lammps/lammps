SHELL=/bin/sh
#########################
# adjust as needed
CXX=g++
CXXFLAGS=-Wall -O2
# set to .exe for windows
EXT=
#########################

all: abf_integrate$(EXT)

clean:
	-rm *~ *.o abf_integrate$(EXT) *.exe

abf_integrate$(EXT): abf_integrate.o abf_data.o
	$(CXX) -o $@ $(CXXFLAGS) $^

%.o: %.cpp
	$(CXX) -o $@ -c $(CXXFLAGS) $<

# dependencies
abf_integrate.o: abf_integrate.cpp abf_data.h
abf_data.o: abf_data.cpp abf_data.h

