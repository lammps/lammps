
TARGET	 = msi2lmp.exe

SRCS     = msi2lmp.c \
           ReadCarFile.c \
           ReadMdfFile.c \
	   MakeLists.c \
           ReadFrcFile.c \
           InitializeItems.c \
	   SearchAndFill.c \
           GetParameters.c \
	   CheckLists.c \
           WriteDataFile.c

OBJS     = $(SRCS:.c=.o)

HEADERS  = msi2lmp.h Forcefield.h

CC       = gcc
CFLAGS   = -O -Wall -W -g 
LDFLAGS  =
FRCFILE  = cvff.frc
FRCFILE2 = cff91.frc
README   = README
MKFILE   = Makefile

$(TARGET) : $(OBJS)
	$(CC) $(LDFLAGS) $(CFLAGS) -o $(TARGET) $(OBJS) -lm

.c.o:
	$(CC) $(CFLAGS) -c $<
clean:
	rm -f $(OBJS) $(TARGET)

# dependencies
CheckLists.o: CheckLists.c msi2lmp.h
GetParameters.o: GetParameters.c msi2lmp.h Forcefield.h
InitializeItems.o: InitializeItems.c msi2lmp.h Forcefield.h
MakeLists.o: MakeLists.c msi2lmp.h
msi2lmp.o: msi2lmp.c msi2lmp.h
ReadCarFile.o: ReadCarFile.c msi2lmp.h
ReadFrcFile.o: ReadFrcFile.c msi2lmp.h Forcefield.h
ReadMdfFile.o: ReadMdfFile.c msi2lmp.h
SearchAndFill.o: SearchAndFill.c msi2lmp.h Forcefield.h
WriteDataFile.o: WriteDataFile.c msi2lmp.h Forcefield.h
