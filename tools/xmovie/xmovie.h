#ifndef PRIVATE
#ifdef USEPRIVATE
#define PRIVATE static
#else
#define PRIVATE
#endif
#endif

#define VERSION	"Version 9.0"

#define NCOLORS		4
#define MAXCOLORS	24
#define MAXDIAM		100
#define MAXTHICK	100

#define MagicNumber	12344321
#define	TimeID		100
#define BoundID		101
#define PositionID	102
#define BondID		103
#define CopyBondID	104
#define CopyAtomID	105

typedef long	INT4;
typedef float	REAL4;

typedef INT4	RECORD;
typedef INT4	LENGTH;

typedef struct {
	INT4	index;
	INT4	type;
	REAL4	x;
	REAL4	y;
	REAL4	z;
	} POSITION;

typedef struct {
	INT4		type;
	INT4		index1;
	INT4		index2;
	POSITION	*atom1;
	POSITION	*atom2;
	} BOND;

typedef struct	{
	REAL4	low;
	REAL4	high;
	} BOUND;

typedef struct {
	REAL4		time;
	BOUND		bounds[3];
	INT4		natoms;
	INT4		maxatoms;
	POSITION	*positions;
	INT4		*natypes;
	INT4		nbonds;
	INT4		maxbonds;
	BOND		*bonds;
	INT4		*nbtypes;
	} DATA;

typedef struct {
	Bool		*atoms_visible;
	Bool		*bonds_visible;
	Bool		hollow;
	Bool		opaque;
	Bool		two_d;
        Bool            pbc_bond;
        Bool            remap;
        Bool            scaleflag;
        Bool            copy_bond;
	Bool		version;
	int		natomcolors;
	int		nbondcolors;
	Dimension	*diameter;
	int		init;
	int		motion;
	int		saveflag;
	unsigned long	delay;
	int		next_pos;
	int		step;
	int		dstep;
	int		axis;
	int		direction;
	BOUND		bounds[3];
	int		ndata;
	int		maxdata;
	DATA		*dataptr;
	float		position;
	float		thickness;
	float		offset[3];
	float		scale;
	} CommonData;

extern CommonData	Common;
extern Widget		TopLevel;
extern char		*Progname;

Widget	CreateScene(Widget parent, char *name);
Widget	CreateControl(Widget parent, char *name);

void	InitRead(int argc, char **argv);
Boolean	ReadProc(XtPointer client_data);

void	RemoveMotion(void);
void	InstallMotion(void);
void	SceneUpdate(void);
void	SceneSave(void);
void	SetTime(char *s);
int	Usage(void);

void	PositionUpdate(void);
void	SpeedUpdate(void);
void	ThicknessUpdate(void);
void	UpdateRadios(void);
void	SceneSize(Dimension *width, Dimension *height);
void	Setup(void);
void	NewDataSetup(void);

void	ExposeScene(Widget w, XEvent *event, String *strings, 
		Cardinal *nstrings);

void	ExposeAxes(Widget w, XEvent *event, String *strings, 
		Cardinal *nstrings);

int	CoerceStep(int step);

void	SetReadString(char *s);

void	SetAtomColors(Pixel *fg);
void	SetBondColors(Pixel *fg, Dimension *thick);
void	SetBGColor(Pixel bg);

void	Version(void);

#ifdef MISSINGDEFS

/* commented this out for SGI gcc compiler
#ifdef stdin
int fprintf(FILE *file, char *fmt, ... );
int fclose(FILE *file);
int printf(char *fmt, ...);
int fflush(FILE *f);
int pclose(FILE *);
int fread(void *, size_t, size_t, FILE *);
int fseek(FILE *, long int, int);
int ungetc(int c, FILE *f);
#endif
*/

#if 0
void XtTranslateCoords(Widget w, Position x, Position y, Position *rx,
	Position *ry);
#endif

#define EXIT_SUCCESS	0
#define EXIT_FAILURE	1

#endif
