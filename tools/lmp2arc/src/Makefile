CC=gcc
OBJS = lmp2arc.o ReadCarFile.o ProcessPosFile01.o ProcessPosFile05.o unwrap_molecules.o WriteArcFrame.o

lmp2arc.exe: $(OBJS)
	$(CC) -o lmp2arc.exe $(OBJS)

lmp2arc.o: lmp2arc.c lmp2.h
ReadCarFile.o: ReadCarFile.c lmp2.h
ProcessPosFile01.o: ProcessPosFile01.c lmp2.h
ProcessPosFile05.o: ProcessPosFile05.c lmp2.h
WriteArcFrame.o: WriteArcFrame.c lmp2.h
unwrap_molecules.o: unwrap_molecules.c lmp2.h

.c.o:
	$(CC) -O2 -c $<

clean: 
	rm lmp2arc.exe *.o
