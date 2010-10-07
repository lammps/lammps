#
#

CXXFLAGS=-O2 -Wall -W

txt2html: txt2html.cpp
	$(CXX) -o $@ $(CXXFLAGS) $^

install: txt2html
	cp -pv $^ $(HOME)/bin/

clean:
	rm -f txt2html txt2html.o *~
