/* **************************************************************** */
/* tobinary.c - binary xmovie format converter                      */
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
#define BondID		103
#define CopyBondID	104
#define CopyAtomID	105

#define LINELEN		256

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef MISSINGDEFS
long strtol(char *s, char **t, int base);
double strtod(char *s, char **t);
int fprintf(FILE *, char *, ...);
int fwrite(void *, size_t, size_t, FILE *);
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

#define WHITESPACE	" \t\n\v\r"
#define SEPARATORS	WHITESPACE ",;"

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
	char	*item;
	RECORD	record_id;
	PFV	reader;
	}	PARSER;

typedef struct {
	RECORD	record_id;
	PFV	writer;
	}	WRITER;

/* **************************************************************** */
/* function proto-types */

int	main(int argc, char **argv);
int	GetRecord(FILE *file);
void	PutRecord(FILE *file);

void	ReadTime(FILE *f);
void	ReadBound(FILE *f);
void	ReadPosition(FILE *f);
void	ReadBond(FILE *);
void	ReadDummy(FILE *F);

int	GetLine(char *s, FILE *f);
void	UnGetLine(char *s);
int	LineIsBlank(char *s);
int	IsItem(char *s, char **t);
PARSER	*GetParser(char *s);

int	GetInt4(char *s, char **t, INT4 *i);
int	GetReal4(char *s, char **t, REAL4 *r);

void	PrintError(char *s);
void	*Realloc(void *ptr, size_t amt);

WRITER	*GetWriter(RECORD record);

void	WriteHeader(FILE *f);
void	WriteRecordHeader(FILE *f);
void	WriteTime(FILE *f);
void	WriteBound(FILE *f);
void	WritePosition(FILE *f);
void	WriteBond(FILE *f);
void	WriteDummy(FILE *f);

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

static PARSER	parse_table[] = {
	{ "TIME",	TimeID,		ReadTime },
	{ "BOUNDS",	BoundID,	ReadBound },
	{ "POSITIONS",	PositionID,	ReadPosition },
	{ "BONDS",	BondID,		ReadBond },
	{ "COPYBONDS",	CopyBondID,	ReadDummy },
	{ "COPYATOMS",	CopyAtomID,	ReadDummy },
	};

static WRITER	write_table[] = {
	{ TimeID,	WriteTime },
	{ BoundID,	WriteBound },
	{ PositionID,	WritePosition },
	{ BondID,	WriteBond },
	{ CopyBondID,	WriteDummy },
	{ CopyAtomID,	WriteDummy },
	};

/* **************************************************************** */

/* ARGSUSED */
int main(int argc, char **argv)
{
	WriteHeader(stdout);

	while(GetRecord(stdin))
		PutRecord(stdout);

	exit(EXIT_SUCCESS);
}

/* **************************************************************** */

int GetRecord(FILE *f)
{
	char	line[LINELEN];
	char	*t;
	PARSER	*p;

	while(1) {
		if (!GetLine(line, f)) return(FALSE);
		if (!IsItem(line, &t)) continue;
		p = GetParser(t);
		if (p == (PARSER *) NULL) continue;
		record = p->record_id;
		(*p->reader)(f);
		break;
	}
	return(TRUE);
}

/* **************************************************************** */

static int	have_line = 0;
static char	hold_line[LINELEN];

int GetLine(char *s, FILE *f)
{
	int	result;

	if (have_line) {
		have_line = FALSE;
		strcpy(s, hold_line);
		return(TRUE);
	}

	do {
		result = (NULL != fgets(s, LINELEN, f));
	}
	while ( result && LineIsBlank(s) );

	return(result);
}

void UnGetLine(char *s)
{
	have_line = TRUE;
	strcpy(hold_line, s);
}

int LineIsBlank(char *s)
{
	register char	*t;

	t = s + strspn(s, WHITESPACE);
	return(*t == '\0');
}

/* **************************************************************** */

int IsItem(char *s, char **t)
{
	register char	*item;

	item = strtok(s, WHITESPACE );
	if (strcmp(item, "ITEM:")) return(FALSE);

	*t = strtok(NULL, WHITESPACE);
	return(TRUE);
}

/* **************************************************************** */

PARSER *GetParser(char *s)
{
	register PARSER	*p;
	register int	i;

	i = sizeof(parse_table)/sizeof(parse_table[0]);
	for(p = parse_table; i ; i--, p++)
		if (!strcmp(p->item, s)) return(p);

	return((PARSER *) NULL);
}

/* **************************************************************** */

void ReadTime(FILE *f)
{
	char	line[80];
	char	*t;

	if (!GetLine(line, f)) 
		PrintError("Error: unable to get line for time.\n");

	if (!GetReal4(line, &t, &time.time))
		PrintError("Error: unable to convert time.\n");

}

/* **************************************************************** */

void ReadBound(FILE *f)
{
	char		line[80];
	char		*t;
	register int	i;

	for(i = 0; i < 3; i++) {
		if (!GetLine(line, f))
			PrintError("Error: unable to get line for bounds.\n");
		if (!GetReal4(line, &t, &bounds.low[i]))
			PrintError("Error: unable to get low bound.\n");
		if (!GetReal4(t, &t, &bounds.high[i]))
			PrintError("Error: unable to get high bound.\n");
	}
}

/* **************************************************************** */

void ReadPosition(FILE *f)
{
	char		line[LINELEN];
	char		*t;	
	POSITION	p;
	register int	i;

	npositions = 0;

	while(1) {
		if (!GetLine(line, f)) return;
		if (!GetInt4(line, &t, &p.index)) {
			UnGetLine(line);
			return;
		}
		if (!GetInt4(t, &t, &p.type))
			PrintError("Error: unable to get atoms type.\n");
		for(i = 0; i < 3; i++)
			if (!GetReal4(t, &t, &p.coord[i]))
				PrintError("Error: unable to atom position.\n");

		if (npositions >= maxpositions) {
			maxpositions += 128;
			positions = (POSITION *) Realloc(positions, 
				maxpositions * sizeof(*positions));
		}
		positions[npositions++] = p;
	}
}

/* **************************************************************** */

void ReadBond(FILE *f)
{
	char		line[LINELEN];
	char		*t;	
	BOND		b;

	nbonds = 0;

	while(1) {
		if (!GetLine(line, f)) return;
		if (!GetInt4(line, &t, &b.type)) {
			UnGetLine(line);
			return;
		}
		if (!GetInt4(t, &t, &b.index1))
			PrintError("Error: unable to get bond index 1.\n");

		if (!GetInt4(t, &t, &b.index2))
			PrintError("Error: unable to get bond index 2.\n");

		if (nbonds >= maxbonds) {
			maxbonds += 128;
			bonds = (BOND *) Realloc(bonds,
				maxbonds * sizeof(*bonds));
		}
		bonds[nbonds++] = b;
	}
}
/* **************************************************************** */

void ReadDummy(FILE *f)
{}

/* **************************************************************** */

int GetInt4(char *s, char **t, INT4 *i)
{
	s += strspn(s, SEPARATORS);

	*i = strtol(s, t, 10);

	return(*t > s);
}

int GetReal4(char *s, char **t, REAL4 *r)
{
	s += strspn(s, SEPARATORS);

	*r = strtod(s, t);

	return(*t > s);
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
	WRITER	*w;

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


void WriteHeader(FILE *f)
{
	fwrite(&magic, sizeof(magic), 1, f);
	fwrite(&version, sizeof(version), 1, f);
}

/* **************************************************************** */

void WriteTime(FILE *f)
{
	length = sizeof(time.time);

	WriteRecordHeader(f);

	fwrite(&time.time, length, 1, f);
}

/* **************************************************************** */

void WriteBound(FILE *f)
{
	register int	i;

	length = 3 * (sizeof(bounds.low[0]) + sizeof(bounds.high[0]));

	WriteRecordHeader(f);

	for(i = 0; i < 3; i++) {
		fwrite(&bounds.low[i], sizeof(bounds.low[0]), 1, f);
		fwrite(&bounds.high[i], sizeof(bounds.high[0]), 1, f);
	}
}

/* **************************************************************** */

void WritePosition(FILE *f)
{
	register int		i;
	register POSITION	*p;

	length = npositions * 
		(sizeof(p->index) + sizeof(p->type) + 3*sizeof(p->coord[0]));

	WriteRecordHeader(f);

	for(i = npositions, p = positions; i; i--, p++) {
		fwrite(&p->index, sizeof(p->index), 1, f);
		fwrite(&p->type, sizeof(p->type), 1, f);
		fwrite(&p->coord[0], sizeof(p->coord[0]), 1, f);
		fwrite(&p->coord[1], sizeof(p->coord[0]), 1, f);
		fwrite(&p->coord[2], sizeof(p->coord[0]), 1, f);
	}
}

/* **************************************************************** */

void WriteBond(FILE *f)
{
	register int		i;
	register BOND		*b;

	length = nbonds * 
		(sizeof(b->type) + sizeof(b->index1) + sizeof(b->index2));

	WriteRecordHeader(f);

	for(i = nbonds, b = bonds; i; i--, b++) {
		fwrite(&b->type,   sizeof(b->type),   1, f);
		fwrite(&b->index1, sizeof(b->index1), 1, f);
		fwrite(&b->index2, sizeof(b->index2), 1, f);
	}
}

/* **************************************************************** */

void WriteDummy(FILE *f)
{
	length = 0;

	WriteRecordHeader(f);
}

/* **************************************************************** */

void WriteRecordHeader(FILE *f)
{
	fwrite(&record, sizeof(record), 1, f);
	fwrite(&length, sizeof(length), 1, f);
}
