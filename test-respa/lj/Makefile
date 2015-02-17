
ALLIN=$(wildcard in.*)
ALLOUT=$(ALLIN:in.%=log.%)
ALLOMP=$(ALLIN:in.%=log.%-omp1) $(ALLIN:in.%=log.%-omp4)
LMP=../../src/lmp_openmpi-omp

all : $(ALLOUT) # $(ALLOMP)


log.% : in.%
	$(LMP) -in $< -log $@

clean:
	rm -f log.* *~

