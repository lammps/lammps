/* ****************************************************************
 * xmovie - a simple X based movie program
 *
 * Mike Uttormark - 7/13/92
 * Sandia Nat'l Labs 1421
 * On leave from University of Wisconsin--Madison
*/

#include <stdio.h>
#include <stdlib.h>

#include <X11/StringDefs.h>

#include <X11/Intrinsic.h>
#include <X11/Xaw/Form.h>
#include <X11/Xaw/Cardinals.h>

#include "xmovie.h"
#include "resource.h"

extern FILE    *popen(const char *, const char *);
extern int     pclose(FILE *);

/* **************************************************************** */
/* local function prototypes */

int		main(int argc, char **argv);
PRIVATE void	CheckResources(void);

/* **************************************************************** */
/* common data */

CommonData	Common = {
	NULL,	/* atoms visible */
	NULL,	/* bonds visible */
	FALSE,	/* hollow */
	FALSE,	/* opaque */
	FALSE,	/* 2d mode */
	FALSE,	/* periodic boundary conditions on bond drawing */
	FALSE,	/* remap atoms into box bounds if necessary */
	FALSE,	/* scale atoms to fill box bounds */
	FALSE,	/* bond copying */
	FALSE,	/* version */
	4,	/* number of atom types */
	4,	/* number of bond types */
	NULL,	/* atom diameters */
	0,	/* init */
	0,	/* motion off */
	0,	/* saveflag off */
	250L,	/* delay interval (ms) */
	0,	/* next drawing position */
	0,	/* step number */
	1,	/* dstep */
	2,	/* z axis */
	0,	/* view direction */
	{ { 1e20, 1e-20 }, { 1e20, 1e-20 }, { 1e20, 1e-20 } }, /* bounds */
	0,	/* ndata */
	0,	/* maxdata */
	NULL,
	};

Widget	TopLevel;
char	*Progname;


/* **************************************************************** */
/* local data */

static XrmOptionDescRec options[] = {
	{ "-2D",	"*twoDimensions",	XrmoptionNoArg,	"True"	},
	{ "-2d",	"*twoDimensions",	XrmoptionNoArg,	"True"	},
	{ "-pbc",	"*pbcBond",		XrmoptionNoArg,	"True"	},
	{ "-remap",	"*remap",		XrmoptionNoArg,	"True"	},
	{ "-scale",	"*scale",		XrmoptionNoArg,	"True"	},
	{ "-copy",	"*copyBond",		XrmoptionNoArg,	"True"	},
	{ "-hollow",	"*hollow",		XrmoptionNoArg,	"True"	},
	{ "-opaque",	"*opaque",		XrmoptionNoArg,	"True"	},
	{ "-V",		"*version",		XrmoptionNoArg,	"True"	},
	{ "-atomcolors","*atomColors",		XrmoptionSepArg, NULL	},
	{ "-bondcolors","*bondColors",		XrmoptionSepArg, NULL	}, 
	};

#define Offset(x,y) (((char *) &(x.y)) - ((char *) &(x)))

static XtResource	resources[] = {
	{	"hollow",		"Hollow",
		XtRBool,		sizeof(Bool),
		Offset(Common, hollow),
		XtRImmediate,		(XtPointer) FALSE	},

	{	"opaque",		"Opaque",
		XtRBool,		sizeof(Bool),
		Offset(Common, opaque),
		XtRImmediate,		(XtPointer) FALSE	},

	{ 	"twoDimensions",	"TwoDimensions",
		XtRBool,		sizeof(Bool),
		Offset(Common, two_d),
		XtRImmediate,		(XtPointer) FALSE	},

	{ 	"pbcBond",       	"PbcBond",
		XtRBool,		sizeof(Bool),
		Offset(Common, pbc_bond),
		XtRImmediate,		(XtPointer) FALSE	},

	{ 	"remap",         	"Remap",
		XtRBool,		sizeof(Bool),
		Offset(Common, remap),
		XtRImmediate,		(XtPointer) FALSE	},

	{ 	"scale",         	"Scale",
		XtRBool,		sizeof(Bool),
		Offset(Common, scaleflag),
		XtRImmediate,		(XtPointer) FALSE	},

	{ 	"copyBond",       	"CopyBond",
		XtRBool,		sizeof(Bool),
		Offset(Common, copy_bond),
		XtRImmediate,		(XtPointer) FALSE	},

	{ 	"version",		"Version",
		XtRBool,		sizeof(Bool),
		Offset(Common, version),
		XtRImmediate,		(XtPointer) FALSE 	},

	{	"atomColors",		"AtomColors",
		XtRInt,			sizeof(int),
		Offset(Common,natomcolors),
		XtRImmediate,		(XtPointer) NCOLORS	},

	{	"bondColors",		"BondColors",
		XtRInt,			sizeof(int),
		Offset(Common,nbondcolors),
		XtRImmediate,		(XtPointer) NCOLORS	},

	};

static XtActionsRec	actions[] = {
	{ "ExposeScene",	ExposeScene },
	{ "ExposeAxes",		ExposeAxes },
	};

int main(int argc, char **argv)
{
	XtAppContext	AppContext;     /* the whole process variable */

	/* Get X stuff going */

	Progname = argv[0];

	TopLevel = XtAppInitialize(&AppContext, "XMovie", options, 
			XtNumber(options), &argc, argv, 
			FallbackResources, NULL, ZERO);

	XtVaGetApplicationResources(TopLevel, (XtPointer) &Common,
			resources, XtNumber(resources), NULL);

	CheckResources();

	XtAppAddActions(AppContext, actions, XtNumber(actions));

	/* initialize the reading stuff */

	InitRead(argc-1, argv+1);

	/* Create the scene box */

	(void) CreateScene(TopLevel,"scene");

	/* Create the Control panel */

	(void) CreateControl(TopLevel,"control");

	XtRealizeWidget(TopLevel);

	/* put reading into background */

	XtAppAddWorkProc(AppContext, ReadProc, (XtPointer) NULL);

	Setup();
	Common.init = 1;

	/* Enter Event Loop */

	XtAppMainLoop(AppContext);
}

int Usage(void)
{
	static char	*msg[] = {
	  "xmovie - a simple & fast atom/molecule visualizer",
	  "    written by Mike Uttormark while at Sandia, 1992",
	  "    updated by Steve Plimpton, Sandia National Labs",
	  "    contact info: sjplimp@sandia.gov",
	  "",
	  "Usage: xmovie [-Xoptions] [-options] file1 [file2 ...]",
	  "       xmovie [-Xoptions] [-options] < file",
	  "",
	  "where -Xoptions are any standard XToolkit options",
	  "-options are any of",
	  "    -2d            reads two dimensional data files",
	  "    -hollow        draws hollow atoms",
	  "    -opaque        draws opaque hollow atoms (implies -hollow)",
	  "    -V             prints version number and exits",
	  "    -pbc           does not draw bonds > 1/2 box size",
	  "    -remap         remap atoms back into bounding box if necessary",
	  "    -scale         scale atom positions to fill bounding box",
	  "    -copy          copies bonds from one timestep to next",
	  "    -atomcolors n  sets number of colors for atoms",
	  "    -bondcolors n  sets number of colors for bonds",
	  "",
	  "and 'files' are data files in the format",
	  "",
	  "    ITEM: TIMESTEP",
	  "        time",
	  "    ITEM: BOX BOUNDS",
	  "        xlow xhigh",
	  "        ylow yhigh",
	  "        zlow zhigh",
	  "    ITEM: BONDS",
	  "        type1 atom1a atom1b",
	  "        type2 atom2a atom2b",
	  "        etc.",
	  "    ITEM: ATOMS",
	  "        index1 type1 x1 y1 z1",
	  "        index2 type2 x2 y2 z2",
	  "        etc.",
	  "",
	  "Notes:",
	  "  A TIMESTEP item starts a new frame",
	  "  A new item (or EOF) ends the list of ATOMS or BONDS",
	  "  BOX BOUNDS persist unless reset",
	  "  Any items not matching these patterns are skipped",
	  "  Data files can be gzipped - e.g. xmovie dump.gz",

	  NULL,
	};
		
	char	**s;
	char	*pager;
	FILE	*pipe;

	/* find user's pager */

	pager = getenv("PAGER");
	if (!pager || !*pager) pager = "more";
	pipe = popen(pager, "w");
	if (!pipe) pipe = stderr;	
	
	for(s = msg; *s; s++)
		fprintf(pipe,"%s\n", *s);

	if (pipe != stderr) pclose(pipe);

	return(0);

}

PRIVATE void CheckResources(void)
{
	if (Common.version) Version();

	if (Common.natomcolors < 0 || Common.natomcolors > MAXCOLORS)
		Common.natomcolors = 1;

	if (Common.nbondcolors < 0 || Common.nbondcolors > MAXCOLORS)
		Common.nbondcolors = 1;

	if (Common.opaque) Common.hollow = TRUE;

	Common.atoms_visible =
		 (Bool *) XtMalloc(Common.natomcolors * sizeof(Bool));

	Common.diameter = 
		(Dimension *) XtMalloc(Common.natomcolors * sizeof(Dimension));

	Common.bonds_visible =
		 (Bool *) XtMalloc(Common.nbondcolors * sizeof(Bool));

}


