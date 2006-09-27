/* **************************************************************** */
/* functions to deal with drawing atoms, etc. */

#include <stdio.h>

#include <X11/Intrinsic.h>
#include <X11/StringDefs.h>
#include <X11/Shell.h>

#include <X11/Xaw/Label.h>
#include <X11/Xaw/Cardinals.h>

#include <X11/xpm.h>

#include "xmovie.h"

/* hard offsets so we dont draw on exact edge of screen */

#define HOFFSET	5
#define VOFFSET	5

typedef XArc		POINT;
typedef XSegment	LINE;

/* **************************************************************** */
/* local prototypes */

PRIVATE LINE	*ClipAndScaleLines(BOND *b, int cnt, int *nlines);
PRIVATE POINT	*ClipAndScalePoints(POSITION *p, int cnt, Dimension diam,
	int *npoints);
PRIVATE void	DrawPoints(Display *display, Drawable drawable, GC gc, 
	POINT *points, int npoints);
PRIVATE void	DrawLines(Display *display, Drawable drawable, GC gc,
	LINE *lines, int nlines);
PRIVATE void	SceneDraw(void);
PRIVATE void	SetAllColors(void);
PRIVATE void	MotionUpdate(void);

/* **************************************************************** */
/* local data */

static Widget	Scene;
static Widget	Shell;
static Pixmap	Buffer;

static GC	*gc = (GC *) NULL;
static GC	*bondgc = (GC *) NULL;
static GC	bggc;
static int	first = 1;
static Dimension	Height, Width;

XtIntervalId	TimeOut = (XtIntervalId) NULL;
static String	null_string = "";

/* **************************************************************** */

Widget CreateScene(Widget parent, char *name)
{
        Shell = XtCreatePopupShell("sceneshell", transientShellWidgetClass,
		parent, NULL, ZERO);

	Scene = XtVaCreateManagedWidget(name, labelWidgetClass, Shell,
		XtNlabel, &null_string, NULL);

	gc = (GC *) XtMalloc(Common.natomcolors * sizeof(GC));
	bondgc = (GC *) XtMalloc(Common.nbondcolors * sizeof(GC));

	XtPopup(Shell, XtGrabNone);

	return(Scene);
}

PRIVATE void SceneDraw(void)
{
	DATA		*dptr;
	POINT		*points;
	LINE		*lines;
	POSITION	*p;
	BOND		*b;
	int		i;
	int		npoints, nlines;
	char		s[40];
	static int	last_step = -1;

	if (!XtIsRealized(Scene)) return;

	if (first) {
		SetAllColors();
		first = 0;
	}

	/* Clear everything */

	if (Common.step >= Common.ndata) return;
	if (Common.step < 0) return;

	XFillRectangle(XtDisplay(Scene), Buffer, bggc, 0, 0, 
			Width, Height);

	/* XClearWindow(XtDisplay(Scene), XtWindow(Scene)); */

	/* find the data */

	dptr = Common.dataptr + Common.step;

	/* loop over colors */

	b = dptr->bonds;
	for(i = 0; i < Common.nbondcolors; i++) {
		if (Common.bonds_visible[i]) {
			lines = ClipAndScaleLines(b, dptr->nbtypes[i],&nlines);
			DrawLines(XtDisplay(Scene), Buffer, bondgc[i], 
				lines, nlines);
		}
		b += dptr->nbtypes[i];
	}

	p = dptr->positions;
	for(i = 0; i < Common.natomcolors; i++){
		if (Common.atoms_visible[i]) {
			points = ClipAndScalePoints(p, dptr->natypes[i], 
					Common.diameter[i], &npoints);
			DrawPoints(XtDisplay(Scene), Buffer, gc[i], points, 
				npoints);
		}
		p += dptr->natypes[i];
	}

	XCopyArea(XtDisplay(Scene), Buffer, XtWindow(Scene), bggc, 
			0, 0, Width, Height, 0, 0);

	XFlush(XtDisplay(Scene));

	if (Common.step == last_step) return;
	last_step = Common.step;
			
	sprintf(s,"Time: %g  Frame: %i", dptr->time, Common.step+1);
	SetTime(s);
}

PRIVATE void MotionUpdate(void)
{
	int	next_step;

	SceneUpdate();

	if (Common.saveflag) {
	  char str[24];
	  if (Common.step < 10)
	    sprintf(str,"image00%d.xpm",Common.step);
	  else if (Common.step < 100)
	    sprintf(str,"image0%d.xpm",Common.step);
	  else
	    sprintf(str,"image%d.xpm",Common.step);
	  XpmWriteFileFromPixmap(XtDisplay(Scene),str,Buffer,NULL,NULL);
	}
	
	if (!Common.motion) return;
	
	next_step = CoerceStep(Common.step + Common.dstep);

	if (next_step == Common.step) {
		RemoveMotion();
		return;
	}

	Common.step = next_step;

	TimeOut = (XtIntervalId) NULL;

	InstallMotion();
}

void SceneUpdate(void)
{

	if (!Common.init) return;

	/* Common.step = CoerceStep(Common.step); */

	SceneDraw();
}

void SceneSave(void)
{
  char str[24];
  sprintf(str,"image.%d.xpm",Common.step);
  XpmWriteFileFromPixmap(XtDisplay(Scene),str,Buffer,NULL,NULL);
}

int CoerceStep(int step)
{
	if (step >= Common.ndata) return (Common.ndata-1);
	if (step < 0) return(0);
	return(step);
}
		
void ExposeScene(Widget w, XEvent *event, String *strings, 
	Cardinal *nstrings)
{
	NewDataSetup();
	SceneDraw();
}


void InstallMotion(void)
{
	Common.motion = 1;

	if (TimeOut == (XtIntervalId) NULL)
		TimeOut = XtAppAddTimeOut(XtWidgetToApplicationContext(Scene),
				Common.delay, 
				(XtTimerCallbackProc) MotionUpdate, NULL);
}

void RemoveMotion(void)
{
	if (!Common.motion) return;

	Common.motion = 0;
	Common.step = CoerceStep(Common.step - Common.dstep);

	if (TimeOut != (XtIntervalId) NULL) XtRemoveTimeOut(TimeOut);
	TimeOut = (XtIntervalId) NULL;
}
	
void SceneSize(Dimension *width, Dimension *height)
{
	XtVaGetValues(Scene, XtNwidth, width, XtNheight, height, NULL);
}

PRIVATE POINT *ClipAndScalePoints(POSITION *pos, int cnt, Dimension diam,
	int *npoints)
{
	register int		i;
	BOUND			range;
	register POINT		*p;
	Dimension		width, height;
	static int		max_points = 0;
	static POINT		*points = (POINT *) NULL;

	range.low  = Common.position - Common.thickness * 0.5;
	range.high = Common.position + Common.thickness * 0.5;

	/* use a static buffer so minimize number of allocations */
	/* didnt use Realloc because dont want data copied */

	*npoints = cnt;
	if (*npoints > max_points) {
		XtFree((char *) points);
		points = (POINT *) XtMalloc( *npoints * sizeof(POINT));
		max_points = *npoints;
	}
	p = points;

	if (cnt < 1) return(points);

	SceneSize(&width, &height);
	width  -= 2*HOFFSET;
	height -= 2*VOFFSET;

	/* translate x, y, z to x and y in window. */
	/* note: index and type are also passed, but may not be used */

	switch(3*Common.direction + Common.axis){
	case 0:			/* negative direction, x axis */
		for(i = cnt; i; i--, pos++) {
			if (pos->x > range.high) continue;
			if (pos->x < range.low)  continue;
			p->x = width  - 
				(pos->z * Common.scale + Common.offset[2]);
			p->y = height - 
				(pos->y * Common.scale + Common.offset[1]);
			p->width = diam;
			p->height = diam;
			p++;
			}	
		break;
	case 1:			/* negative direction, y axis */
		for(i = cnt; i; i--, pos++) {
			if (pos->y > range.high) continue;
			if (pos->y < range.low)  continue;
			p->x = pos->x * Common.scale + Common.offset[0];
			p->y = pos->z * Common.scale + Common.offset[2];
			p->width = diam;
			p->height = diam;
			p++;
			}	
		break;
	case 2:			/* negative direction, z axis */
		for(i = cnt; i; i--, pos++) {
			if (pos->z > range.high) continue;
			if (pos->z < range.low)  continue;
			p->x = pos->x * Common.scale + Common.offset[0];
			p->y = height -
				 (pos->y * Common.scale + Common.offset[1]);
			p->width = diam;
			p->height = diam;
			p++;
			}	
		break;
	case 3:			/* positive direction, x axis */
		for(i = cnt; i; i--, pos++) {
			if (pos->x > range.high) continue;
			if (pos->x < range.low)  continue;
			p->x = pos->z * Common.scale + Common.offset[2];
			p->y = height - 
				(pos->y * Common.scale + Common.offset[1]);
			p->width = diam;
			p->height = diam;
			p++;
			}	
		break;

	case 4:			/* positive direction, y axis */
		for(i = cnt; i; i--, pos++) {
			if (pos->y > range.high) continue;
			if (pos->y < range.low)  continue;
			p->x = pos->x * Common.scale + Common.offset[0];
			p->y = height - 
				(pos->z * Common.scale + Common.offset[2]);
			p->width = diam;
			p->height = diam;
			p++;
			}	
		break;


	case 5:			/* postive direction, z axis */
		for(i = cnt; i; i--, pos++) {
			if (pos->z > range.high) continue;
			if (pos->z < range.low)  continue;
			p->x = pos->x * Common.scale + Common.offset[0];
			p->y = pos->y * Common.scale + Common.offset[1];
			p->width = diam;
			p->height = diam;
			p++;
			}	
		break;
	}

	*npoints = p - points;

	/* add the hard offsets so we dont draw on edge of screen */
	/* center drawing based on width, height */

	for(i = *npoints, p = points; i; i--, p++) {
		p->x += HOFFSET - p->width/2;
		p->y += VOFFSET - p->height/2;
	}

	return(points);
}
		
	
PRIVATE void DrawPoints(Display *display, Drawable drawable, GC gc, 
	POINT *points, int npoints)
{
	register int		full_circle;
	register int		i;
	register POINT		*p;

	if (npoints <= 0) return;

	/* this version has POINT typedef'd to XArc, so we need to */
	/* fill in other fields before drawing */

	full_circle = 64*360;

	for(i = npoints, p = points; i ; i--, p++){
		p->angle1 = 0;
		p->angle2 = full_circle;
	}

	if (Common.hollow) {
		if (Common.opaque)
			XFillArcs(display, drawable, bggc, points, npoints);
		XDrawArcs(display, drawable, gc, points, npoints);
	}
	else
		XFillArcs(display, drawable, gc, points, npoints);

}

PRIVATE LINE *ClipAndScaleLines(BOND *bond, int cnt, int *nlines)
{
	register int		i;
	BOUND			range;
	register LINE		*l;
	Dimension		width, height;
	static int		max_lines = 0;
	static LINE		*lines = (LINE *) NULL;

	range.low  = Common.position - Common.thickness * 0.5;
	range.high = Common.position + Common.thickness * 0.5;

	/* use a static buffer so minimize number of allocations */
	/* didnt use Realloc because dont want data copied */

	*nlines = cnt;
	if (*nlines > max_lines) {
		XtFree((char *) lines);
		lines = (LINE *) XtMalloc( *nlines * sizeof(LINE));
		max_lines = *nlines;
	}
	l = lines;

	if (cnt < 1) return(lines);

	SceneSize(&width, &height);
	width  -= 2*HOFFSET;
	height -= 2*VOFFSET;

	/* translate x, y, z to x and y in window. */

	switch(3*Common.direction + Common.axis){
	case 0:			/* negative direction, x axis */
		for(i = cnt; i; i--, bond++) {
			if (bond->atom1->x > range.high) continue;
			if (bond->atom1->x < range.low)  continue;
			if (bond->atom2->x > range.high) continue;
			if (bond->atom2->x < range.low)  continue;
			l->x1 = width  - 
				(bond->atom1->z*Common.scale+Common.offset[2]);
			l->y1 = height - 
				(bond->atom1->y*Common.scale+Common.offset[1]);
			l->x2 = width  - 
				(bond->atom2->z*Common.scale+Common.offset[2]);
			l->y2 = height - 
				(bond->atom2->y*Common.scale+Common.offset[1]);
			l++;
			}	
		break;
	case 1:			/* negative direction, y axis */
		for(i = cnt; i; i--, bond++) {
			if (bond->atom1->y > range.high) continue;
			if (bond->atom1->y < range.low)  continue;
			if (bond->atom2->y > range.high) continue;
			if (bond->atom2->y < range.low)  continue;
			l->x1 = bond->atom1->x*Common.scale + Common.offset[0];
			l->y1 = bond->atom1->z*Common.scale + Common.offset[2];
			l->x2 = bond->atom2->x*Common.scale + Common.offset[0];
			l->y2 = bond->atom2->z*Common.scale + Common.offset[2];
			l++;
			}	
		break;
	case 2:			/* negative direction, z axis */
		for(i = cnt; i; i--, bond++) {
			if (bond->atom1->z > range.high) continue;
			if (bond->atom1->z < range.low)  continue;
			if (bond->atom2->z > range.high) continue;
			if (bond->atom2->z < range.low)  continue;
			l->x1 = bond->atom1->x*Common.scale + Common.offset[0];
			l->y1 = height -
				(bond->atom1->y*Common.scale+Common.offset[1]);
			l->x2 = bond->atom2->x*Common.scale + Common.offset[0];
			l->y2 = height -
				(bond->atom2->y*Common.scale+Common.offset[1]);
			l++;
			}	
		break;
	case 3:			/* positive direction, x axis */
		for(i = cnt; i; i--, bond++) {
			if (bond->atom1->x > range.high) continue;
			if (bond->atom1->x < range.low)  continue;
			if (bond->atom2->x > range.high) continue;
			if (bond->atom2->x < range.low)  continue;
			l->x1 = bond->atom1->z*Common.scale + Common.offset[2];
			l->y1 = height - 
				(bond->atom1->y*Common.scale+Common.offset[1]);
			l->x2 = bond->atom2->z*Common.scale + Common.offset[2];
			l->y2 = height - 
				(bond->atom2->y*Common.scale+Common.offset[1]);
			l++;
			}	
		break;

	case 4:			/* positive direction, y axis */
		for(i = cnt; i; i--, bond++) {
			if (bond->atom1->y > range.high) continue;
			if (bond->atom1->y < range.low)  continue;
			if (bond->atom2->y > range.high) continue;
			if (bond->atom2->y < range.low)  continue;
			l->x1 = bond->atom1->x*Common.scale + Common.offset[0];
			l->y1 = height - 
				(bond->atom1->z*Common.scale+Common.offset[2]);
			l->x2 = bond->atom2->x*Common.scale + Common.offset[0];
			l->y2 = height - 
				(bond->atom2->z*Common.scale+Common.offset[2]);
			l++;
			}	
		break;


	case 5:			/* postive direction, z axis */
		for(i = cnt; i; i--, bond++) {
			if (bond->atom1->z > range.high) continue;
			if (bond->atom1->z < range.low)  continue;
			if (bond->atom2->z > range.high) continue;
			if (bond->atom2->z < range.low)  continue;
			l->x1 = bond->atom1->x*Common.scale + Common.offset[0];
			l->y1 = bond->atom1->y*Common.scale + Common.offset[1];
			l->x2 = bond->atom2->x*Common.scale + Common.offset[0];
			l->y2 = bond->atom2->y*Common.scale + Common.offset[1];
			l++;
			}	
		break;
	}

	*nlines = l - lines;

	/* add the hard offsets so we dont draw on edge of screen */

	for(i = *nlines, l = lines; i; i--, l++) {
		l->x1 += HOFFSET;
		l->y1 += VOFFSET;
		l->x2 += HOFFSET;
		l->y2 += VOFFSET;
	}

	return(lines);
}
		
	
PRIVATE void DrawLines(Display *display, Drawable drawable, GC gc, 
	LINE *lines, int nlines)
{
	if (nlines <= 0) return;

	XDrawSegments(display, drawable, gc, lines, nlines);
}

void Setup(void)
{
	NewDataSetup();
	SpeedUpdate();
	UpdateRadios();
}
	
void NewDataSetup(void)
{
	static int	have_pixmap = 0;

	BOUND		*cb;
	int		i;
	float		longest, f;
	Dimension	width, height;

	SceneSize(&Width, &Height);

	if (have_pixmap) XFreePixmap(XtDisplay(Scene), Buffer);

	Buffer = XCreatePixmap(XtDisplay(Scene), 
		RootWindowOfScreen(XtScreen(Scene)), Width, Height,
		DefaultDepthOfScreen(XtScreen(Scene)));

	have_pixmap = 1;

	/* determine global scaling and offset factors */
	/* offset + scale * coordinate = pixel */

	longest = 0.0;
	for(i = 3, cb = Common.bounds; i ; i--, cb++)
		if ((f = (cb->high) - (cb->low)) > longest) longest = f;

	SceneSize(&Width, &Height);
	width = Width;
	height = Height;

	width  -= HOFFSET*2;
	height -= VOFFSET*2;

	Common.scale = (width < height) ? width/longest : height/longest;

	for(i = 0; i < 3; i++)
		Common.offset[i] = - Common.bounds[i].low * Common.scale;
	
	PositionUpdate();
	ThicknessUpdate();
}

/* **************************************************************** */

void SetAtomColors(Pixel *fg)
{
	int		i;
	XGCValues	xgc;

	XtVaGetValues(Scene, XtNbackground, &xgc.background, NULL);

	xgc.function = GXcopy;

	for(i = 0; i < Common.natomcolors; i++) {
		if (gc[i]) XFreeGC(XtDisplay(Scene), gc[i]);
		xgc.foreground = fg[i];
		gc[i] = XCreateGC(XtDisplay(Scene), XtWindow(Scene),
			GCFunction | GCBackground | GCForeground, &xgc);
	}

	if (bggc) XFreeGC(XtDisplay(Scene), bggc);
	xgc.foreground = xgc.background;
	bggc = XCreateGC(XtDisplay(Scene), XtWindow(Scene),
		GCFunction | GCBackground | GCForeground, &xgc);

}

void SetBGColor(Pixel bg)
{
	XGCValues	xgc;

	XtVaSetValues(Scene, XtNbackground, bg, NULL);

	if (bggc) XFreeGC(XtDisplay(Scene), bggc);

	xgc.function = GXcopy;
	xgc.foreground = xgc.background = bg;

	bggc = XCreateGC(XtDisplay(Scene), XtWindow(Scene),
		GCFunction | GCBackground | GCForeground, &xgc);

}

/* **************************************************************** */

void SetBondColors(Pixel *fg, Dimension *thick)
{
	int		i;
	XGCValues	xgc;

	XtVaGetValues(Scene, XtNbackground, &xgc.background, NULL);
	
	xgc.function = GXcopy;

	for(i = 0; i < Common.nbondcolors; i++) {
		if (bondgc[i]) XFreeGC(XtDisplay(Scene), bondgc[i]);
		xgc.foreground = fg[i];
		xgc.line_width = thick[i];
		bondgc[i] = XCreateGC(XtDisplay(Scene), XtWindow(Scene),
			GCFunction | GCBackground | GCForeground | GCLineWidth,
			&xgc);
	}
}

/* **************************************************************** */

PRIVATE void SetAllColors(void)
{
	int		i;
	XGCValues	xgc;

	xgc.function = GXcopy;
	xgc.line_width = 1;

	XtVaGetValues(Scene, 
			XtNbackground, &xgc.background,
			XtNforeground, &xgc.foreground, 
			NULL);

	for(i = 0; i < Common.natomcolors; i++)
		gc[i] = XCreateGC(XtDisplay(Scene), XtWindow(Scene),
			GCFunction | GCBackground | GCForeground, &xgc);

	for(i = 0; i < Common.nbondcolors; i++)
		bondgc[i] = XCreateGC(XtDisplay(Scene), XtWindow(Scene),
			GCFunction | GCBackground | GCForeground | GCLineWidth,
			&xgc);

	xgc.foreground = xgc.background;
	bggc = XCreateGC(XtDisplay(Scene), XtWindow(Scene),
			GCFunction | GCBackground | GCForeground, &xgc);

}

