/* **************************************************************** */
/* toascii.c - binary to ascii xmovie format converter              */
/*                                                                  */
/* Mike Uttormark - 8/11/92                                         */
/* Sandia Nat'l Labs 1421                                           */
/* On leave from University of Wisconsin--Madison                   */
/* **************************************************************** */
/* format of binary file:                                           */
/*                                                                  */
/* magic number   - 4-byte integer                                  */
/* version        - 4-byte integer                                  */
/* <records>                                                        */
/*                                                                  */
/* generic record format:                                           */
/* record id      - 4-byte integer                                  */
/* length (bytes) - 4-byte integer                                  */
/* data                                                             */
/* **************************************************************** */

#define MagicNumber 	12344321
#define Version		1001

#define TimeID		100
#define BoundID		101
#define PositionID	102
#define	BondID		103
#define CopyBondID	104
#define CopyAtomID	105

#define LINELEN		256

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef MISSINGDEFS
#include <unistd.h>
int fprintf(FILE *, char *, ...);
int fread(void *, size_t, size_t, FILE *);
int fseek(FILE *, long int, int);
#endif

#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS	0
#endif
#ifndef EXIT_FAILURE
#define EXIT_FAILURE	1
#endif

#undef TRUE
#undef FALSE
#define TRUE	1
#define FALSE	0

/* **************************************************************** */
/* local typedefs */

typedef long	INT4;
typedef float	REAL4;
typedef void	(*PFV)();

typedef struct {
	INT4	index;
	INT4	type;
	REAL4	coord[3];
	}	POSITION;

typedef struct {
	INT4	type;
	INT4	index1;
	INT4	index2;
	}	BOND;

typedef struct {
	REAL4	low[3];
	REAL4	high[3];
	}	BOUND;

typedef struct {
	REAL4	time;
	}	TIME;

typedef INT4	RECORD;
typedef INT4	LENGTH;
typedef INT4	MAGIC;
typedef INT4	VERSION;

typedef struct {
	RECORD	record_id;
	PFV	writer;
	}	WRITER;

typedef struct {
	RECORD	record_id;
	PFV	reader;
	}	READER;

/* **************************************************************** */
/* function proto-types */

int	main(int argc, char **argv);
int	GetRecord(FILE *file);
void	PutRecord(FILE *file);

void	PutTime(FILE *f);
void	PutBound(FILE *f);
void	PutPosition(FILE *f);
void	PutBond(FILE *f);
void	PutCopyBond(FILE *f);
void	PutCopyAtom(FILE *f);

WRITER	*GetWriter(RECORD record);

void	PrintError(char *s);
void	*Realloc(void *ptr, size_t amt);

READER	*GetReader(RECORD record);

void	ReadHeader(FILE *f);
int	ReadRecordHeader(FILE *f);
void	ReadTime(FILE *f);
void	ReadBound(FILE *f);
void	ReadPosition(FILE *f);
void	ReadBond(FILE *f);
void	ReadDummy(FILE *f);

/* **************************************************************** */
/* local data */

static MAGIC	magic = MagicNumber;
static VERSION	version = Version;
static RECORD	record;
static LENGTH	length;
static POSITION	*positions = (POSITION *) NULL;
static INT4	npositions = 0;
static INT4	maxpositions = 0;
static BOUND	bounds;
static TIME	time;
static BOND	*bonds = (BOND *) NULL;
static INT4	nbonds = 0;
static INT4	maxbonds = 0;

static WRITER	write_table[] = {
	{ TimeID,	PutTime },
	{ BoundID,	PutBound },
	{ PositionID,	PutPosition },
	{ BondID,	PutBond },
	{ CopyBondID,	PutCopyBond },
	{ CopyAtomID,	PutCopyAtom },
	};

static READER	read_table[] = {
	{ TimeID,	ReadTime },
	{ BoundID,	ReadBound },
	{ PositionID,	ReadPosition },
	{ BondID,	ReadBond },
	{ CopyBondID,	ReadDummy },
	{ CopyAtomID,	ReadDummy },
	};

/* **************************************************************** */

/* ARGSUSED */
int main(int argc, char **argv)
{
	ReadHeader(stdin);

	while(GetRecord(stdin))
		PutRecord(stdout);

	exit(EXIT_SUCCESS);
}

/* **************************************************************** */

int GetRecord(FILE *f)
{
	register READER	*r;

	while(1) {
		if (!ReadRecordHeader(f)) return(FALSE);
		r = GetReader(record);
		if (r == (READER *) NULL) {		/* dont recognize */
			fseek(f, length, SEEK_CUR);	/* skip */
			continue;
		}
		(*r->reader)(f);
		break;
	}
	return(TRUE);
}

/* **************************************************************** */

READER *GetReader(RECORD record)
{
	register READER	*r;
	register int	i;

	i = sizeof(read_table)/sizeof(read_table[0]);
	for(r = read_table; i ; i--, r++)
		if (r->record_id == record) return(r);

	return((READER *) NULL);
}

/* **************************************************************** */

void PutTime(FILE *f)
{
	fprintf(f, "ITEM: TIME\n"
		"%g\n", time.time);

}

/* **************************************************************** */

void PutBound(FILE *f)
{
	register int	i;

	fprintf(f, "ITEM: BOUNDS\n");

	for(i = 0; i < 3; i++) 
		fprintf(f, "%g %g\n", bounds.low[i], bounds.high[i]);
}

/* **************************************************************** */

void PutPosition(FILE *f)
{
	register POSITION	*p;
	register int		i;

	fprintf(f,"ITEM: POSITIONS %i\n", npositions);

	for(p = positions, i = npositions; i; i--, p++)
		fprintf(f, "%i %i %g %g %g\n",p->index, p->type,
			p->coord[0], p->coord[1], p->coord[2]);
}

/* **************************************************************** */

void PutBond(FILE *f)
{
	register BOND	*b;
	register int	i;

	fprintf(f,"ITEM: BONDS %i\n", nbonds);

	for(b = bonds, i = nbonds; i; i--, b++)
		fprintf(f, "%i %i %i\n", b->type, b->index1, b->index2);
}

/* **************************************************************** */

void PutCopyBond(FILE *f)
{
	fprintf(f, "ITEM: COPYBONDS\n");
}

/* **************************************************************** */

void PutCopyAtom(FILE *f)
{
	fprintf(f, "ITEM: COPYATOMS\n");
}

/* **************************************************************** */

void PrintError(char *s)
{
	fprintf(stderr,"%s", s);
	exit(EXIT_FAILURE);
}

/* **************************************************************** */

void *Realloc(void *ptr, size_t amt)
{
	ptr = (ptr == NULL) ? malloc(amt) : realloc(ptr, amt);

	if (ptr != NULL) return(ptr);

	PrintError("Error: unable to allocate space.\n");
}

/* **************************************************************** */

void PutRecord(FILE *f)
{
	register WRITER	*w;

	w = GetWriter(record);
	if (w == (WRITER *) NULL)
		PrintError("Internal error: no writer.\n");

	(*w->writer)(f);
}

/* **************************************************************** */

WRITER *GetWriter(RECORD r)
{
	register int	i;
	register WRITER	*w;

	i = sizeof(write_table)/sizeof(write_table[0]);
	for(w = write_table; i; i--, w++)
		if (w->record_id == r) return(w);

	return((WRITER *) NULL);
}

/* **************************************************************** */


void ReadHeader(FILE *f)
{
	if (1 > fread(&magic, sizeof(magic), 1, f)) goto error;
	if (1 > fread(&version, sizeof(version), 1, f)) goto error;

	if (magic == MagicNumber) return;

	error:

	PrintError(	"Error: magic number not found.\n"
			"File is probably not binary.\n");
}

/* **************************************************************** */

void ReadTime(FILE *f)
{
	if (1 == fread(&time.time, length, 1, f)) return;

	PrintError("Error: unable to read time.\n");
}

/* **************************************************************** */

void ReadBound(FILE *f)
{
	register int	i;

	for(i = 0; i < 3; i++) {
		if (1 > fread(&bounds.low[i], sizeof(bounds.low[0]), 1, f))
			goto error;
		if (1 > fread(&bounds.high[i], sizeof(bounds.high[0]), 1, f))
			goto error;
	}
	return;

	error:
		PrintError("Error: unable to read bounds.\n");
}

/* **************************************************************** */

void ReadPosition(FILE *f)
{
	register int		i;
	register POSITION	*p;

	npositions = length / 
		(sizeof(p->index) + sizeof(p->type) + 3*sizeof(p->coord[0]));

	if (npositions > maxpositions) {
		maxpositions = npositions + 16;
		positions = (POSITION *)
			Realloc(positions, maxpositions * sizeof(*positions));
	}

	for(i = npositions, p = positions; i; i--, p++) {
		if (1 > fread(&p->index, sizeof(p->index), 1, f)) goto error;
		if (1 > fread(&p->type, sizeof(p->type), 1, f)) goto error;
		if (1 > fread(&p->coord[0], sizeof(p->coord[0]), 1, f)) goto error;
		if (1 > fread(&p->coord[1], sizeof(p->coord[0]), 1, f)) goto error;
		if (1 > fread(&p->coord[2], sizeof(p->coord[0]), 1, f)) goto error;
	}
	return;

	error:
		PrintError("Error: unable to get atom positions.\n");
}

/* **************************************************************** */

void ReadBond(FILE *f)
{
	register int	i;
	register BOND	*b;

	nbonds = length / 
		(sizeof(b->type) + sizeof(b->index1) + sizeof(b->index2));

	if (nbonds > maxbonds) {
		maxbonds = nbonds + 16;
		bonds = (BOND *)
			Realloc(bonds, maxbonds * sizeof(*bonds));
	}

	for(i = nbonds, b = bonds; i; i--, b++) {
		if (1 > fread(&b->type,   sizeof(b->type),   1, f)) goto error;
		if (1 > fread(&b->index1, sizeof(b->index1), 1, f)) goto error;
		if (1 > fread(&b->index2, sizeof(b->index2), 1, f)) goto error;
	}
	return;

	error:
		PrintError("Error: unable to get bonds.\n");
}

/* **************************************************************** */

int ReadRecordHeader(FILE *f)
{
	if (1 > fread(&record, sizeof(record), 1, f)) return(FALSE);
	if (1 > fread(&length, sizeof(length), 1, f)) return(FALSE);
}

/* **************************************************************** */

void ReadDummy(FILE *f)
{}
