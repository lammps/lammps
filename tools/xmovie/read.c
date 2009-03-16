/* **************************************************************** */
/* functions to read data from files */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

/* commented this out for SGI gcc compiler
#ifdef INCL_FLOAT
double	strtod(char *s, char **t);
long	strtol(char *s, char **t, int base);
#endif
*/

#include <X11/Intrinsic.h>

extern FILE    *popen(const char *, const char *);
extern int     pclose(FILE *);

#include "xmovie.h"

#define LINELEN		256
#define WHITESPACE	" \f\n\r\t\v"
#define SEPARATORS	",;:" WHITESPACE
#define NPRINT		500
#define NREAD		50
#define NBREAD		2000

#define ERROR			-1
#define END_OF_FILE		1
#define BEING_NICE		2
#define READ_OK			3
#define LOOK_FOR_ITEM		10
#define LOOK_FOR_BOUNDS		11
#define LOOK_FOR_POSITIONS	12
#define LOOK_FOR_TIME		13
#define LOOK_FOR_BONDS		14

#define BAD_STATE		-2

#define ALL_DONE		TRUE
#define MORE_TO_READ		FALSE

/* **************************************************************** */
/* local typedefs */

typedef struct {
	char	*string;
	int	state;
	}	PARSETABLE;

#define	NPARSE(x)	(sizeof(x)/sizeof(PARSETABLE))

typedef struct {
	int	nfiles;
	char	**fnames;
	FILE	*file;
	int	in_read;
	int	is_compressed;
	int	is_binary;
	}	READDATA;

typedef int	(*PFI)();

typedef struct {
	RECORD	record_id;
	PFI	reader;
	}	READER;

/* **************************************************************** */
/* local prototypes */

PRIVATE int	DoSomeReading(void);
PRIVATE int	ReadSomeBinary(void);
PRIVATE int	parser(char *token, PARSETABLE *table, int cnt);
PRIVATE int	get_bounds(FILE *file, char *line);
PRIVATE int	get_positions(FILE *file, char *line);
PRIVATE int	get_time(FILE *file, char *line);
PRIVATE int	get_bonds(FILE *file, char *line);
PRIVATE int	copy_bonds(void);
PRIVATE void	add_point(int index, int type, float x, float y, float z);
PRIVATE void	add_bond(int type, int index1, int index2);
PRIVATE void	NewData(void);
PRIVATE	void	ShowStatus(void);
PRIVATE char	*RemoveNL(char *s);
PRIVATE int	getline(FILE *file, char *line);
PRIVATE int	get_float(char **s, float *f);
PRIVATE int	get_int(char **s, int *i);
PRIVATE void	ungetline(void);
PRIVATE int	LineIsBlank(char *s);
PRIVATE void	SortData(DATA *dptr);
PRIVATE int	compare_atypes(POSITION *a, POSITION *b);
PRIVATE int	compare_btypes(BOND *a, BOND *b);
PRIVATE int	IsCompressed(char *s);
PRIVATE int	IsBinary(char *s);
PRIVATE int	TestBinary(FILE *f, int rewind);
PRIVATE FILE	*GetCompressed(char *s);
PRIVATE FILE	*OpenFile(char *s, int compressed, int binary);
PRIVATE READER	*GetReader(RECORD record_id);
PRIVATE int	ReadTime(FILE *f, LENGTH length);
PRIVATE int	ReadBound(FILE *f, LENGTH length);
PRIVATE int	ReadPosition(FILE *f, LENGTH length);
PRIVATE int	ReadBond(FILE *f, LENGTH length);
PRIVATE int	CopyBond(FILE *f, LENGTH length);
PRIVATE int	ReadRecordHeader(FILE *f, RECORD *record, LENGTH *length);

/* **************************************************************** */
/* local data */

/* include last 3 out-of-date item names in PARSETABLE for compatability */

static PARSETABLE	parse_table[] = {
	{ "TIMESTEP",	LOOK_FOR_TIME },
	{ "BOX BOUNDS",	LOOK_FOR_BOUNDS },
	{ "ATOMS",	LOOK_FOR_POSITIONS },
	{ "BONDS",	LOOK_FOR_BONDS },
	{ "TIME",	LOOK_FOR_TIME },
	{ "BOUNDS",	LOOK_FOR_BOUNDS },
	{ "POSITIONS",	LOOK_FOR_POSITIONS },
	};

static DATA	*current_data = (DATA *) NULL;
static int	seen_time = 0;
static int	npositions = 0;
static int	nlines = 0;
static int	nbonds = 0;
static int	read_this_time = 0;

static READDATA	ReadData = {
	0,	/* nfiles */
	NULL,	/* filenames */
	NULL,	/* file pointer */
	0,	/* not in read */
	0,	/* not compressed */
	0,	/* not binary */
	};

static READER	read_table[] = {
	{ TimeID,	ReadTime },
	{ BoundID,	ReadBound },
	{ PositionID,	ReadPosition },
	{ BondID,	ReadBond },
	{ CopyBondID,	CopyBond },
	};

/* **************************************************************** */

void InitRead(int nfiles, char **fnames)
{
	static char	*Stdin = "stdin";

	if (nfiles > 0) {
		ReadData.nfiles = nfiles;
		ReadData.fnames = fnames;
		return;
	}

	if (isatty(fileno(stdin))) {
		Usage();
		exit(0);
	}

	ReadData.nfiles	= 1;
	ReadData.fnames	= &Stdin;
	ReadData.file	= stdin;
	ReadData.in_read = 1;
	ReadData.is_binary = TestBinary(stdin, TRUE);

	if (ReadData.is_binary) 
		fseek(stdin, 2*sizeof(INT4), SEEK_CUR);	/* skip header */
}

Boolean ReadProc(XtPointer client_data)
{
	char	err_msg[512];
	char	*file_type;

	if (ReadData.in_read) goto do_read;

	next_file:

		if (ReadData.nfiles <= 0) {
			ShowStatus();
			SortData(current_data);
			return(ALL_DONE);
		}

		ReadData.is_compressed = FALSE;
		ReadData.is_binary     = FALSE;

		ReadData.is_compressed = IsCompressed(*ReadData.fnames);

		if (ReadData.is_compressed) {
			file_type = "compressed";
			goto open_file;
		}

		ReadData.is_binary = IsBinary(*ReadData.fnames);
		if (ReadData.is_binary) {
			file_type = "binary";
			goto open_file;
		}
			

		file_type = "ascii";
		goto open_file;

	open_file:

		ReadData.file = OpenFile(*ReadData.fnames,
			ReadData.is_compressed, ReadData.is_binary);
		if (ReadData.file == (FILE *) NULL) {
			sprintf(err_msg, 
				"Frames: %i  Positions: %i  Bonds: %i\n"
				"Error: unable to open %s file \"%s\".\n"
				"Read aborted.",
				Common.ndata, npositions, nbonds,
				file_type, *ReadData.fnames);
				SetReadString(err_msg);
				SortData(current_data);
				return(ALL_DONE);
		}

		ReadData.in_read = 1;
		nlines = 0;
		ShowStatus();

	do_read:

		switch(DoSomeReading()) {
		case END_OF_FILE:
			if (ReadData.is_compressed) 
				pclose(ReadData.file);
			else
				if (strcmp("stdin",*ReadData.fnames))
					fclose(ReadData.file);
			ReadData.in_read = 0;
			ReadData.nfiles--;
			ReadData.fnames++;
			goto next_file;
		case ERROR:
			SortData(current_data);
			return(ALL_DONE);
		case BEING_NICE:
			break;
		default:
			fprintf(stderr,"Don't know how I got here: ReadProc.\n");
			exit(0);
		}

	return(MORE_TO_READ);
}

PRIVATE int IsCompressed(char *s)
{
	register char	*t;

	t = s + strlen(s) - 3;
	if (t < s) return(FALSE);
	return(!strcmp(t,".gz"));
}

PRIVATE int IsBinary(char *s)
{
	FILE	*f;
	int	ok;

	if ((FILE *) NULL == (f = fopen(s, "r"))) return(FALSE);

	ok = TestBinary(f, FALSE);
	fclose(f);

	return(ok);
}

PRIVATE int TestBinary(FILE *f, int rewind)
{
	union {
		INT4	magic;
		char	c[sizeof(INT4)];
	} u;
	int	i;

	if (1 != fread(&u.magic, sizeof(u.magic), 1, f)) return(FALSE);

	if (rewind)
		for(i = sizeof(u); i; i--)
			ungetc(u.c[i-1], f);
		
	return(u.magic == MagicNumber);
}

PRIVATE FILE *GetCompressed(char *name)
{
	char	cmd[256];

	sprintf(cmd, "gunzip -c %s", name);
	return(popen(cmd,"r"));
}

PRIVATE FILE *OpenFile(char *name, int compressed, int binary)
{
	FILE	*f;

	if (compressed) return(GetCompressed(name));

	f = fopen(name, "r");
	if (f == (FILE *) NULL) return(f);

	if (binary) fseek(f, 2*sizeof(INT4), SEEK_CUR);

	return(f);
}

PRIVATE int DoSomeReading(void)
{
	static int	state = LOOK_FOR_ITEM;
	static char	line[LINELEN];
	static char	token[LINELEN];
	static char	*t;
	int		ntries;
	int		newstate;
	char		err_msg[512];

	if (ReadData.is_binary) return(ReadSomeBinary());

	read_this_time = 0;
	
	begin:

	switch(state) {
	case LOOK_FOR_ITEM:
		ntries = 0;
		while(1) {
			ntries++;
			if (ntries > 100) return(BEING_NICE);

			if (!getline(ReadData.file, line))
				return(END_OF_FILE);
			t = strtok(line, WHITESPACE);
			if (!t) continue;
			if (strcmp(t, "ITEM:")) continue;

			t = strtok(NULL, WHITESPACE);
			if (!t) continue;
			sprintf(token,"%s",t);

			t = strtok(NULL, WHITESPACE);
			while (t != NULL) {
			  sprintf(token,"%s %s",token,t);
			  t = strtok(NULL, WHITESPACE);
			}

			newstate = parser(token,
					  parse_table,NPARSE(parse_table));

			if (newstate == BAD_STATE) continue;
			state = newstate;
			goto begin;
		}
	case LOOK_FOR_TIME:
		if (!get_time(ReadData.file, line)) goto print_error;
		state = LOOK_FOR_ITEM;
		goto begin;
	case LOOK_FOR_BOUNDS:
		if (!get_bounds(ReadData.file, line)) goto print_error;
		state = LOOK_FOR_ITEM;
		goto begin;
	case LOOK_FOR_POSITIONS:
	case LOOK_FOR_BONDS:
		switch((state == LOOK_FOR_POSITIONS) ?
			get_positions(ReadData.file, line) : 
			get_bonds(ReadData.file, line)) {
		case END_OF_FILE:
			state = LOOK_FOR_ITEM;
			return(END_OF_FILE);
		case BEING_NICE:
			return(BEING_NICE);
		case LOOK_FOR_ITEM:
			state = LOOK_FOR_ITEM;
			return(BEING_NICE);
		case ERROR:
		default:
			goto print_error;
		}

	case BAD_STATE:
	default:
		fprintf(stderr,"Error: read state corrupted: %i\n", state);
		exit(0);
	}

	print_error:

	sprintf(err_msg,"Frames: %i  Positions: %i  Bonds: %i\n"
			"Error occurred at line %i\n"
			"Processing ITEM: %s\n"
			"Line = \"%s\"\n"
			"Read aborted.",
			Common.ndata, npositions, nbonds, nlines, token, 
			RemoveNL(line));

	SetReadString(err_msg);

	return(ERROR);

	/* exit(0); */

}

/* **************************************************************** */

PRIVATE int parser(char *token, PARSETABLE *table, int cnt)
{
	PARSETABLE	*t;
	int		i;

	for(i = cnt, t = table; i; i--, t++)
	  if (strstr(token,t->string) == token)
	    return t->state;
	return(BAD_STATE);
}

/* **************************************************************** */

PRIVATE int get_bounds(FILE *file, char *line)
{
	int		i, n, changed;
	char		*t;
	BOUND		*b, *cb;

	if (!current_data) NewData();

	n = (Common.two_d) ? 2 : 3;

	for(i = 0; i < n; i++){
		if (!getline(file, line)) return(0);

		t = line;

		if (!get_float(&t, &(current_data->bounds[i].low)))
			return(0);

		if (!get_float(&t, &(current_data->bounds[i].high)))
			return(0);
	}

	if (Common.two_d) {
		current_data->bounds[2].low = -0.01;
		current_data->bounds[2].high = 0.01;
	}

	changed = 0;
	for(b=current_data->bounds,cb=Common.bounds,i=3; i;i--,b++,cb++){
		if (b->low < cb->low) {
			cb->low = b->low;
			changed = 1;
		}
		if (b->high > cb->high) {
			cb->high = b->high;
			changed = 1;
		}
	}

	if (changed) NewDataSetup();

	return(1);
}

PRIVATE int get_positions(FILE *file, char *line)
{
	char		*t;
	float		x, y, z;
	float           xlow,xhigh,ylow,yhigh,zlow,zhigh,dx,dy,dz;
	int		index, type;

	if (!current_data) NewData();

	if (Common.copy_bond) copy_bonds();

	if (Common.remap || Common.scaleflag) {
	  xlow = current_data->bounds[0].low;
	  xhigh = current_data->bounds[0].high;
	  dx = xhigh - xlow;
	  ylow = current_data->bounds[1].low;
	  yhigh = current_data->bounds[1].high;
	  dy = yhigh - ylow;
	  if (!Common.two_d) {
	    zlow = current_data->bounds[2].low;
	    zhigh = current_data->bounds[2].high;
	    dz = zhigh - zlow;
	  }
	}

	z = 0.0;

	while(1) {

		if (read_this_time++ > NREAD) return(BEING_NICE);

		if (!getline(file, line)) return(END_OF_FILE);

		t = line;

		if (!get_int(&t, &index)){
			ungetline();
			return(LOOK_FOR_ITEM);
		}

		if (!get_int(&t, &type))  return(ERROR);

		if (!get_float(&t, &x)) return(ERROR);
		if (!get_float(&t, &y)) return(ERROR);

		if (!Common.two_d)
			if (!get_float(&t, &z)) return(ERROR);

		if (type < 1 || type > Common.natomcolors) type = 1;

		if (Common.remap) {
		  if (x-xlow >= 0.0)
		    x = fmod(x-xlow,dx) + xlow;
		  else
		    x = xhigh + fmod(x-xlow,dx);
		  if (y-ylow >= 0.0)
		    y = fmod(y-ylow,dy) + ylow;
		  else
		    y = yhigh + fmod(y-ylow,dy);
		  if (!Common.two_d) {
		    if (z-zlow >= 0.0)
		      z = fmod(z-zlow,dz) + zlow;
		    else
		      z = zhigh + fmod(z-zlow,dz);
		  }
		}

		if (Common.scaleflag) {
		  x = xlow + x*dx;
		  y = ylow + y*dy;
		  if (!Common.two_d) z = zlow + z*dz;
		}

		add_point(index, type, x, y, z);

		npositions++;

		if ((nbonds + npositions) % NPRINT == 0) ShowStatus();

	}

	return(1);
}

PRIVATE int get_bonds(FILE *file, char *line)
{
	char		*t;
	int		type, index1, index2;

	if (!current_data) NewData();

	while(1) {

		if (read_this_time++ > NREAD) return(BEING_NICE);

		if (!getline(file, line)) return(END_OF_FILE);

		t = line;

		if (!get_int(&t, &type)){
			ungetline();
			return(LOOK_FOR_ITEM);
		}

		if (!get_int(&t, &index1)) return(ERROR);
		if (!get_int(&t, &index2)) return(ERROR);

		if (type < 1 || type > Common.nbondcolors) type = 1;

		add_bond(type, index1, index2);

		nbonds++;

		if ((nbonds + npositions) % NPRINT == 0) ShowStatus();
	}

	return(1);
}

PRIVATE int get_time(FILE *file, char *line)
{
	char	*t;

	if (!current_data || seen_time) NewData();
	seen_time = 1;

	if (!getline(file, line)) return(0);

	t = line;

	if (!get_float(&t, &(current_data->time))) return(0);

	return(1);
}
	
PRIVATE int copy_bonds()
{
	DATA	*last_data;
	int	i;
	BOND	*b1, *b2;

	if (Common.ndata < 2) return(0);

	last_data = current_data-1;

	current_data->nbonds = last_data->nbonds;
	current_data->maxbonds = last_data->maxbonds;

	current_data->bonds = 
		(BOND *) XtMalloc(current_data->maxbonds * sizeof(BOND));

	b1 = current_data->bonds;
	b2 = last_data->bonds;
	for(i = current_data->nbonds; i ; i--, b1++, b2++)
		*b1 = *b2;
 
	return(1);
}

/* **************************************************************** */


PRIVATE char *RemoveNL(char *s)
{
	register char	*t;

	if (!s) return(s);

	if ((t = s + strlen(s) - 1) < s) return(s);
	if (*t == '\n') *t = '\0';

	return(s);
}
	
static int got_line = 0;

PRIVATE int getline(FILE *file, char *line)
{
	int	result;

	if (got_line) {
		got_line = 0;
		return(1);
	}

	do {
		nlines++;
		result = (fgets(line, LINELEN, file) != NULL);
	} while( result && LineIsBlank(line) );

	return(result);
}

PRIVATE void ungetline(void)
{
	got_line = 1;
}

PRIVATE int LineIsBlank(char *s)
{
	register char	*t;

	t = s + strspn(s, WHITESPACE);
	return( *t == '\0' );
}

PRIVATE int get_int(char **s, int *i)
{
	register char	*t;

	t = *s + strspn(*s, SEPARATORS);
	*i = strtol(t, s, 10);

	return(t != *s);
}

PRIVATE int get_float(char **s, float *f)
{
	register char	*t;

	t = *s + strspn(*s, SEPARATORS);
	*f = strtod(t, s);

	return(t != *s);
}

/* **************************************************************** */

PRIVATE void NewData(void)
{
	int	i;
	DATA	*last_data,*olddataptr;

	olddataptr = Common.dataptr;
	last_data = current_data;

	if (Common.ndata >= Common.maxdata) {	/* need to realloc */
		Common.maxdata += 256;
		Common.dataptr = (DATA *)
		  XtMalloc(Common.maxdata * sizeof(DATA));
		memcpy(Common.dataptr,olddataptr,Common.ndata*sizeof(DATA));
	}

	current_data = Common.dataptr + Common.ndata;
	Common.ndata++;

	current_data->natypes = 
		(INT4 *) XtMalloc(Common.natomcolors * sizeof(INT4));

	for(i = 0; i < Common.natomcolors; i++)
		current_data->natypes[i] = 0;

	current_data->nbtypes = 
		(INT4 *) XtMalloc(Common.nbondcolors * sizeof(INT4));

	for(i = 0; i < Common.nbondcolors; i++)
		current_data->nbtypes[i] = 0;

	if (Common.ndata > 1) {
		current_data->time = last_data->time;

		for(i = 0; i < 3; i++)
			current_data->bounds[i] = last_data->bounds[i];

		current_data->maxatoms = last_data->natoms;
		current_data->positions = (POSITION *)
			XtMalloc(last_data->natoms * sizeof(POSITION));

		current_data->natoms = 0;

		current_data->maxbonds = last_data->nbonds;
		current_data->bonds = (BOND *)
			XtMalloc(last_data->nbonds * sizeof(BOND));

		current_data->nbonds = 0;

		SortData(last_data);
	}
	else {
		current_data->time = 0.0;

		for(i = 0; i < 3; i++) {
			current_data->bounds[i].low = -1.0;
			current_data->bounds[i].high = 1.0;
		}

		current_data->natoms = 0;
		current_data->positions = NULL;
		current_data->maxatoms = 0;

		current_data->nbonds = 0;
		current_data->bonds = NULL;
		current_data->maxbonds = 0;
	}

	if (olddataptr != Common.dataptr) XtFree((char *) olddataptr);
}

PRIVATE void add_point(int index, int type, float x, float y, float z)
{
	POSITION	*p;

	if (current_data->natoms >= current_data->maxatoms) {
		current_data->maxatoms += 16;
		current_data->positions = (POSITION *)
			XtRealloc((char *) current_data->positions, 
				  current_data->maxatoms * sizeof(POSITION));
	}

	p = &(current_data->positions[current_data->natoms]);

	p->index = index;
	p->type  = type;
	p->x     = x;
	p->y     = y;
	p->z     = z;

	current_data->natoms++;
}


PRIVATE void add_bond(int type, int index1, int index2)
{
	BOND	*b;

	if (current_data->nbonds >= current_data->maxbonds) {
		current_data->maxbonds += 16;
		current_data->bonds = (BOND *)
			XtRealloc((char *) current_data->bonds,
				current_data->maxbonds * sizeof(BOND));
	}

	b = &(current_data->bonds[current_data->nbonds]);

	b->type   = type;
	b->index1 = index1;
	b->index2 = index2;
	b->atom1  = NULL;
	b->atom2  = NULL;

	current_data->nbonds++;
}

PRIVATE void ShowStatus(void)
{
	static char	buffer[256];

	if (ReadData.in_read)
		sprintf(buffer,	"Reading file: %s\n"
				"Frames: %i  Positions: %i  Bonds: %i",
				*ReadData.fnames, Common.ndata, npositions,
				nbonds);
	else
		sprintf(buffer,	"Reading done.\n"
				"Frames: %i  Positions: %i  Bonds %i",
				Common.ndata, npositions, nbonds);

	SetReadString(buffer);
}

/* **************************************************************** */

PRIVATE void SortData(DATA *dptr)
{
	POSITION	*p, **p_by_index, **pp;
	register int	i, j;
	int		maxindex, minindex, nindex, n;
	register BOND	*b, *b2;
	float		bond_length, shortest, dx, dy, dz;

	/* add hpsort in place of standard C qsort, which is very slow
           when atoms or bonds are already in order (N^2 instead of NlnN)
           calling syntax does not change except for first argument 
           becomes arg-1, qsort -> hpsort in 2 places
	*/
	void hpsort(POSITION *, int, int, int (*)(POSITION *, POSITION *));

	if (!dptr) return;

	/* sort by type */

	hpsort((dptr->positions)-1, dptr->natoms, sizeof(POSITION), 
		compare_atypes);

	/* find out how many of each type */

	for(i = dptr->natoms, p = dptr->positions; i ;i--, p++)
		dptr->natypes[p->type-1]++;

	if (dptr->nbonds < 1) return;

	/* find max, min indicies */

	minindex = 1000000000;
	maxindex = -minindex;

	for(i = dptr->natoms, p = dptr->positions; i; i--, p++) {
		if (p->index > maxindex) maxindex = p->index;
		if (p->index < minindex) minindex = p->index;
	}

	nindex = maxindex - minindex + 1;

	if (nindex < 0) {
		dptr->nbonds = 0;
		return;
	}
		  

	if (nindex < 1000000) {

		p_by_index = (POSITION **) XtMalloc(nindex*sizeof(POSITION *));

		for(pp = p_by_index, i = nindex; i; i--, pp++)
			*pp = (POSITION *) NULL;

		for(p = dptr->positions, i = dptr->natoms; i; i--, p++){
			pp = p_by_index + p->index - minindex;
			if (*pp != (POSITION *) NULL) {
				fprintf(stderr,
					"Error: atom %i has more than one position at time %g.\n",
					p->index, dptr->time);
				exit(EXIT_FAILURE);
			}
			*pp = p;
		}

		for(b = dptr->bonds, i = dptr->nbonds; i ; i--, b++) {
			b->atom1 = p_by_index[b->index1 - minindex];
			b->atom2 = p_by_index[b->index2 - minindex];
			if (b->atom1 == (POSITION *) NULL) {
				fprintf(stderr,
					"Error: atom %d has no position at time %g.\n", 
					b->index1, dptr->time);
				exit(EXIT_FAILURE);
			}
			if (b->atom2 == (POSITION *) NULL) {
				fprintf(stderr,
					"Error: atom %d has no position at time %g.\n", 
					b->index2, dptr->time);
				exit(EXIT_FAILURE);
			}
		}

		XtFree((char *) p_by_index);
	}
	else {
		for(b = dptr->bonds, i = dptr->nbonds; i; i--, b++) {
			for(p = dptr->positions, j = dptr->natoms; j; j--, p++)
				if (b->index1 == p->index) {
					b->atom1 = p;
					break;
				}
			if (j < 1) {
				fprintf(stderr, 
					"Error: atom %i has no position at time %g.\n",
					b->index1, dptr->time);
				exit(EXIT_FAILURE);
			}
			for(p = dptr->positions, j = dptr->natoms; j; j--, p++)
				if (b->index2 == p->index) {
					b->atom2 = p;
					break;
				}
			if (j < 1) {
				fprintf(stderr, 
					"Error: atom %i has no position at time %g.\n",
					b->index2, dptr->time);
				exit(EXIT_FAILURE);
			}
		}
	}

	/* throw out bonds longer than half the shortest cell size */
        /* only do if periodic-bond switch is set */

	if (Common.pbc_bond) {

	  n = (Common.two_d) ? 2 : 3;

	  shortest = 1e30;
	  for(i = 0; i < n; i++){
	    bond_length = dptr->bounds[i].high - dptr->bounds[i].low;
	    if (bond_length < shortest) shortest = bond_length;
	  }
	  
	  shortest *= 0.5;
	  shortest *= shortest;
	  
	  b = dptr->bonds;
	  b2 = b;
	  for(i = dptr->nbonds; i; i--, b++) {
	    dx = b->atom1->x - b->atom2->x;
	    dy = b->atom1->y - b->atom2->y;
	    dz = b->atom1->z - b->atom2->z;
	    bond_length = dx*dx + dy*dy + dz*dz;
	    if (bond_length < shortest) {
	      if (b2 != b) *b2 = *b;
	      b2++;
	    }
	  }
	  
	  dptr->nbonds -= b - b2;

	}

	/* sort by type */

	hpsort((dptr->bonds)-1, dptr->nbonds, sizeof(BOND), 
		compare_btypes);

	/* find out how many of each type */

	for(b = dptr->bonds, i = dptr->nbonds; i; i--, b++)
		dptr->nbtypes[b->type-1]++;
}

PRIVATE int compare_atypes(POSITION *a, POSITION *b)
{
	if (a->type  > b->type) return(1);
	if (a->type == b->type) return(0);
	if (a->type  < b->type) return(-1);
}

PRIVATE int compare_btypes(BOND *a, BOND *b)
{
	if (a->type  > b->type) return(1);
	if (a->type == b->type) return(0);
	if (a->type  < b->type) return(-1);
}

/* routine to print out all data for debugging */

#if 0
	
print_them(void)
{
	int	i;

	for(i = 0; i < Common.ndata; i++)
		print_data(Common.dataptr[i]);
}

print_data(DATA *data)
{
	int		i;
	POSITION	*p;

	fprintf(stderr,"Time: %f\n", data->time);
	
	for(i = 0; i < 3; i++)
		fprintf(stderr,"Bounds %i: %f %f\n", i, data->bounds[i].low,
				data->bounds[i].high);

	for(i = 0; i < data->natoms; i++) {
		p = &(data->positions[i]);
		fprintf(stderr,"Atom: %i %i %f %f %f\n",
			p->index, p->type, p->x, p->y, p->z);
	}
}

#endif

/* **************************************************************** */

PRIVATE int ReadSomeBinary(void)
{
	RECORD	record_id;
	LENGTH	length;
	READER	*r;
	char	err_msg[256];
	int	state;
	int	old_positions, old_bonds;

	old_positions = npositions;
	old_bonds = nbonds;

	while(1) {
		if (!ReadRecordHeader(ReadData.file, &record_id, &length))
			return(END_OF_FILE);

		r = GetReader(record_id);
		if (r == (READER *) NULL) {	/* dont recognize */
			fseek(ReadData.file, length, SEEK_CUR);
			continue;
		}

		switch(state = (*r->reader)(ReadData.file, length)) {
			case BEING_NICE:
				if (old_positions + old_bonds + NBREAD > 
					nbonds + npositions)
					continue;
				return(state);
			case END_OF_FILE:
				return(state);
			case READ_OK:
				continue;
			case ERROR:
			default:	
				sprintf(err_msg,"Frames: %i  Positions: %i  Bonds: %i\n"
					"Error occurred during binary read.\n"
					"Read aborted.",
					Common.ndata, npositions, nbonds); 
				SetReadString(err_msg);
				return(ERROR);
		}
	}
}

PRIVATE READER *GetReader(RECORD record_id)
{
	register int	i;
	register READER	*r;

	i = sizeof(read_table)/sizeof(read_table[0]);
	for(r = read_table; i; i--, r++)
		if (r->record_id == record_id) return(r);

	return((READER *) NULL);
}

PRIVATE int ReadRecordHeader(FILE *f, RECORD *record, LENGTH *length)
{
	if (1 != fread(record, sizeof(*record), 1, f)) return(FALSE);
	if (1 != fread(length, sizeof(*length), 1, f)) return(FALSE);

	return(TRUE);
}

/*ARGSUSED*/
PRIVATE int ReadBound(FILE *file, LENGTH length)
{
	int		i, changed;
	BOUND		*b, *cb;

	if (!current_data) NewData();

	

	for(b = current_data->bounds, i = 3; i; i--, b++){
		if (1 != fread(&(b->low), sizeof(b->low), 1, file))
			return(ERROR);
		if (1 != fread(&(b->high), sizeof(b->high), 1, file))
			return(ERROR);
	}

	changed = 0;
	for(b=current_data->bounds,cb=Common.bounds,i=3; i;i--,b++,cb++){
		if (b->low < cb->low) {
			cb->low = b->low;
			changed = 1;
		}
		if (b->high > cb->high) {
			cb->high = b->high;
			changed = 1;
		}
	}

	if (changed) NewDataSetup();

	return(READ_OK);
}

PRIVATE int ReadPosition(FILE *file, LENGTH length)
{
	POSITION	*p;
	int		n;
	int		i;

	if (!current_data) NewData();

	n = length / (sizeof(p->index) + sizeof(p->type) + 3*sizeof(p->x));

	if (current_data->maxatoms < n) {
		free(current_data->positions);
		current_data->maxatoms = n;
		current_data->positions = 
			(POSITION *) XtMalloc(n * sizeof(POSITION));
	}

	current_data->natoms = 0;
	p = current_data->positions;
	i = 0;

	for(; i < n; i++, p++) {
		if (1 != fread(&(p->index), sizeof(p->index), 1, file))
			return(ERROR);
		if (1 != fread(&(p->type),  sizeof(p->type),  1, file))
			return(ERROR);
		if (1 != fread(&(p->x),     sizeof(p->x),     1, file))
			return(ERROR);
		if (1 != fread(&(p->y),     sizeof(p->y),     1, file))
			return(ERROR);
		if (1 != fread(&(p->z),     sizeof(p->z),     1, file))
			return(ERROR);

		if (p->type < 1 || p->type > Common.natomcolors)
			p->type = 1;

		npositions++;
		current_data->natoms++;

		if ((npositions + nbonds)% NPRINT == 0) ShowStatus();
	}

	return(BEING_NICE);
}

PRIVATE int ReadBond(FILE *file, LENGTH length)
{
	BOND		*b;
	int		n;
	int		i;

	if (!current_data) NewData();

	n = length / (sizeof(b->type) + sizeof(b->index1) + sizeof(b->index2));

	if (current_data->maxbonds < n) {
		free(current_data->bonds);
		current_data->maxbonds = n;
		current_data->bonds = 
			(BOND *) XtMalloc(n * sizeof(BOND));
	}

	current_data->nbonds = 0;
	b = current_data->bonds;
	i = 0;

	for(; i < n; i++, b++) {
		if (1 != fread(&(b->type),   sizeof(b->type),    1, file))
			return(ERROR);
		if (1 != fread(&(b->index1), sizeof(b->index1),  1, file))
			return(ERROR);
		if (1 != fread(&(b->index2), sizeof(b->index2),  1, file))
			return(ERROR);

		if (b->type < 1 || b->type > Common.nbondcolors)
			b->type = 1;

		nbonds++;
		current_data->nbonds++;

		if ((npositions + nbonds)% NPRINT == 0) ShowStatus();
	}

	return(BEING_NICE);
}

/*ARGSUSED*/
PRIVATE int CopyBond(FILE *file, LENGTH length)
{
	return( (copy_bonds()) ? READ_OK : ERROR );
}



/*ARGSUSED*/
PRIVATE int ReadTime(FILE *file, LENGTH length)
{
	if (!current_data || seen_time) NewData();
	seen_time = 1;

	if (1 != fread(&(current_data->time), sizeof(current_data->time),
			1, file))
		return(ERROR);
	return(READ_OK);
}
