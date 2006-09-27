/* **************************************************************** */
/* functions to handle control panel (sliders, labels, etc.) */

#include <stdio.h>
#include <stdlib.h>

#include <X11/Intrinsic.h>
#include <X11/StringDefs.h>
#include <X11/Shell.h>
#include <X11/Xos.h>

#include <X11/Xaw/Cardinals.h>

#include <X11/Xaw/Box.h>
#include <X11/Xaw/Command.h>
#include <X11/Xaw/Dialog.h>
#include <X11/Xaw/Form.h>
#include <X11/Xaw/Label.h>
#include <X11/Xaw/MenuButton.h>
#include <X11/Xaw/Scrollbar.h>
#include <X11/Xaw/SimpleMenu.h>
#include <X11/Xaw/Sme.h>
#include <X11/Xaw/SmeBSB.h>
#include <X11/Xaw/Toggle.h>

#include "xmovie.h"

#define POPUPSHELLCLASS overrideShellWidgetClass

/* **************************************************************** */
/* local typedefs */

typedef void	(*InitFunc)(Widget parent);

typedef struct {
	char		*name;		
	XtCallbackProc	exec_func;
	XtPointer	data;
	InitFunc	init_func;
	}		ButtonStruct;

#define NBTN(a)	(sizeof(a)/sizeof(ButtonStruct))

typedef struct {
	char	*name;
	int	value;
	Widget	w;
	}	ToggleStruct;

#define TGLCNT(x) (sizeof(x)/sizeof(ToggleStruct))

typedef struct {
	char		*name;
	XtCallbackProc	callback;
	ToggleStruct	*buttons;
	int		nbuttons;
	Widget		lw;
	}		RadioStruct;

typedef struct {
	char		*labelname;
	Widget		*label;
	char		*slidername;
	Widget		*slider;
	XtCallbackProc	callback;
	} SliderLabelStruct;

typedef struct {
	Widget		*w;
	XtGrabKind	grab;
	Bool		stop_motion;
	} PopupData;

/* **************************************************************** */
/* local prototypes */

PRIVATE void SetSpeed(float percent);
PRIVATE void NewSpeed(Widget w, XtPointer client_data, XtPointer percent);
PRIVATE void SetPosition(float percent);
PRIVATE void NewPosition(Widget w, XtPointer client_data, XtPointer percent);
PRIVATE void SetThickness(float percent);
PRIVATE void NewThickness(Widget w, XtPointer client_data, XtPointer percent);

PRIVATE void AxisSelect(Widget w, XtPointer client_data, XtPointer call_data);
PRIVATE void DirectionSelect(Widget w, XtPointer client_data, XtPointer call_data);
PRIVATE void MotionSelect(Widget w, XtPointer client_data, XtPointer call_data);

PRIVATE void FindAndDisplaySlider(Widget w, XtCallbackProc func);

PRIVATE void DrawAxes(void);

PRIVATE void MakeRadioGroup(Widget parent, RadioStruct *radio);
PRIVATE void MakeSliderWithLabel(Widget parent, SliderLabelStruct *sl);

PRIVATE void SetRadio(RadioStruct *radio);

PRIVATE void do_start(Widget w, XtPointer client_data, XtPointer call_data);
PRIVATE void do_stop(Widget w, XtPointer client_data, XtPointer call_data);
PRIVATE void do_restart(Widget w, XtPointer client_data, XtPointer call_data);
PRIVATE void do_step(Widget w, XtPointer client_data, XtPointer call_data);
PRIVATE void do_back(Widget w, XtPointer client_data, XtPointer call_data);

PRIVATE void do_save_off(Widget w, XtPointer client_data, XtPointer call_data);
PRIVATE void do_save_on(Widget w, XtPointer client_data, XtPointer call_data);
PRIVATE void do_save_one(Widget w, XtPointer client_data, XtPointer call_data);

PRIVATE void do_popup(Widget w, XtPointer client_data, XtPointer call_data);
PRIVATE void cancel_popup(Widget w, XtPointer client_data, XtPointer call_data);

PRIVATE void init_quit(Widget w);
PRIVATE void quit_ok(Widget w, XtPointer client_data, XtPointer call_data);

PRIVATE void init_atoms(Widget w);
PRIVATE void init_bonds(Widget w);
PRIVATE void init_bg(Widget w);
PRIVATE void init_file(Widget w);

PRIVATE void color_apply(Widget w, XtPointer client_data, XtPointer call_data);
PRIVATE void bond_apply(Widget w, XtPointer client_data, XtPointer call_data);
PRIVATE void bg_apply(Widget w, XtPointer client_data, XtPointer call_data);
PRIVATE void file_apply(Widget w, XtPointer client_data, XtPointer call_data);

PRIVATE Pixel ConvertColorToPixel(Widget w, String color);
PRIVATE XtActionProc NoOp(Widget w, XEvent event, String *params,
		Cardinal nparams);

PRIVATE void visible_toggle(Widget w, XtPointer client_data, XtPointer
		call_data);

PRIVATE void bvisible_toggle(Widget w, XtPointer client_data, XtPointer
		call_data);

/* **************************************************************** */
/* local data */

static Widget	QuitPopup	= (Widget) 0,
		AtomPopup	= (Widget) 0,
		BondPopup	= (Widget) 0,
		BackgroundPopup	= (Widget) 0,
		FileColorPopup	= (Widget) 0,
		TimeLabel	= (Widget) 0,
		SpeedBar	= (Widget) 0,
		SpeedLabel	= (Widget) 0,
		PositionBar	= (Widget) 0,
		PositionLabel	= (Widget) 0,
		ThicknessBar	= (Widget) 0,
		ThicknessLabel	= (Widget) 0,
		Axes		= (Widget) 0,
		ReadInfo	= (Widget) 0;

static Widget 	BGDialog	= (Widget) 0,
         	FileDialog	= (Widget) 0,
		*FGDialog	= (Widget *) NULL,
		*Visible	= (Widget *) NULL,
		*Invisible	= (Widget *) NULL,
		*SizeDialog	= (Widget *) NULL,
		*BFGDialog	= (Widget *) NULL,
		*BVisible	= (Widget *) NULL,
		*BInvisible	= (Widget *) NULL,
		*ThicknessDialog= (Widget *) NULL;

static PopupData
	QuitData	= { &QuitPopup,		XtGrabExclusive, True},
	AtomData 	= { &AtomPopup,		XtGrabNone,	False },
	BondData	= { &BondPopup,		XtGrabNone,	False },
	BackgroundData	= { &BackgroundPopup,	XtGrabNone,	False },
        FileColorData	= { &FileColorPopup,	XtGrabNone,	False };

static ButtonStruct Buttons[] = {
	{ "quit",	do_popup,	(XtPointer) &QuitData,	init_quit },
	{ "start",	do_start,	NULL,			NULL },
	{ "stop",	do_stop,	NULL,			NULL },
	{ "restart",	do_restart,	NULL,			NULL },
	{ "step",	do_step,	NULL,			NULL },
	{ "back",	do_back,	NULL,			NULL },
	};

static ButtonStruct SaveButtons[] = {
	{ "Off while animating", do_save_off, NULL, NULL },
	{ "On  while animating", do_save_on , NULL, NULL },
	{ "This snapshot",       do_save_one, NULL, NULL },
	};

static ButtonStruct ColorButtons[] = {
	{ "atoms",	do_popup,	(XtPointer) &AtomData, init_atoms },
	{ "bonds",	do_popup,	(XtPointer) &BondData, init_bonds },
	{ "background", do_popup,	(XtPointer) &BackgroundData, init_bg },
	{ "filecolor", do_popup,	(XtPointer) &FileColorData, init_file }
	};

static int	in_motion;

static SliderLabelStruct
	SpeedSL = {
		"speed",	&SpeedLabel,
		"speedbar",	&SpeedBar,
		NewSpeed },
	PositionSL = {
		"position",	&PositionLabel,
		"positionbar",	&PositionBar,
		NewPosition },
	ThicknessSL = {
		"thickness",	&ThicknessLabel,
		"thicknessbar",	&ThicknessBar,
		NewThickness }
	;
		
static ToggleStruct 
	AxisButtons[] = {
		{ "x",	0 },
		{ "y",	1 },
		{ "z",	2 } },
	DirectionButtons[] = { 
		{ "vplus", 1 },
		{ "vminus", 0 } },
	MotionButtons[]	= { 
		{ "mplus", 1 },
		{ "mminus", -1 } }
	;

static RadioStruct 
	AxisRadio = {
		"axis", 
		AxisSelect,
		AxisButtons,
		TGLCNT(AxisButtons) },
	DirectionRadio = {
		"viewdir",
		DirectionSelect,
		DirectionButtons,
		TGLCNT(DirectionButtons) },
	MotionRadio = {
		"motiondir",
		MotionSelect,
		MotionButtons,
		TGLCNT(MotionButtons) }
	;
		
static int	normal = 1;
static char	*null_string = "";

static XtActionsRec actions[] = {
        { "NoOp",       (XtActionProc) NoOp},
        };

/* **************************************************************** */

Widget CreateControl(Widget parent, char *name)
{
	Widget		form, command, button, menu;
	ButtonStruct	*bptr;
	int		bcnt;

	/* register the new actions */

	XtAppAddActions(XtWidgetToApplicationContext(parent),
		actions, XtNumber(actions));

	/* create a form widget to handle layout */

	form = XtCreateManagedWidget(name, formWidgetClass, parent,
			NULL, ZERO);

	/* create all the command buttons */

	for(bcnt = NBTN(Buttons), bptr = Buttons;bcnt > 0; bcnt--, bptr++){

		/* create a button */

		command = XtCreateManagedWidget(bptr->name,
				commandWidgetClass, form, NULL, ZERO);

		/* Add the function as a callback */

		XtAddCallback(command, XtNcallback, bptr->exec_func, 
				bptr->data);

		/* call initialization routine */

		if (bptr->init_func) (bptr->init_func)(command);

	}

	/* create the save menu */

	command = XtCreateManagedWidget("save", menuButtonWidgetClass,
		form, NULL, ZERO);

	menu = XtCreatePopupShell("menu", simpleMenuWidgetClass, 
		command, NULL, ZERO);

	for(bcnt=NBTN(SaveButtons), bptr=SaveButtons; bcnt>0; bcnt--,bptr++){
		button = XtCreateManagedWidget(bptr->name, smeBSBObjectClass,
				menu, NULL, ZERO);
		XtAddCallback(button, XtNcallback, bptr->exec_func,bptr->data);
		if (bptr->init_func) (bptr->init_func)(button);
	}

	/* create the color menu */

	command = XtCreateManagedWidget("color", menuButtonWidgetClass,
		form, NULL, ZERO);

	menu = XtCreatePopupShell("menu", simpleMenuWidgetClass, 
		command, NULL, ZERO);

	for(bcnt=NBTN(ColorButtons), bptr=ColorButtons; bcnt>0; bcnt--,bptr++){
		button = XtCreateManagedWidget(bptr->name, smeBSBObjectClass,
				menu, NULL, ZERO);
		XtAddCallback(button, XtNcallback, bptr->exec_func,bptr->data);
		if (bptr->init_func) (bptr->init_func)(button);
	}

	/* create a label widget to display time */

	TimeLabel = XtCreateManagedWidget("time", labelWidgetClass,
			form, NULL, ZERO);

	SetTime("Time: 0");

	/* ************************************************ */
	/* create all the sliders                           */
	/* ************************************************ */


	MakeSliderWithLabel(form, &SpeedSL);
	MakeSliderWithLabel(form, &PositionSL);
	MakeSliderWithLabel(form, &ThicknessSL);

	/* ******************************************************** */
	/* create a radio group for axis, view, motion directions   */
	/* ******************************************************** */

	MakeRadioGroup(form, &AxisRadio);
	MakeRadioGroup(form, &DirectionRadio);
	MakeRadioGroup(form, &MotionRadio);

	/* create a label widget to put axes in */

	/* SJP - 10/02, had to add last NULL argument to each call */

	Axes = XtVaCreateManagedWidget("axes", labelWidgetClass, form,
				       XtNlabel, null_string,NULL);

	ReadInfo = XtVaCreateManagedWidget("read", labelWidgetClass,
					   form, XtNlabel, null_string,NULL);

}

PRIVATE XtActionProc NoOp(Widget w, XEvent event, String *params,
		Cardinal nparams)
{
        return;
}

void SetTime(char *s)
{
	XtVaSetValues(TimeLabel, XtNlabel, s, NULL);
}

PRIVATE void SetSpeed(float percent)
{
	char	string[40];
	int	speed;

	speed = percent * 100;
	(void) sprintf(string,"Speed: %i", speed);

	XtVaSetValues(SpeedLabel, XtNlabel, string, NULL);
}
	
PRIVATE void NewSpeed(Widget w, XtPointer client_data, XtPointer percent_ptr)
{
	float	percent;

	percent = * (float *) percent_ptr;

	/* delay between frames in milliseconds */

	Common.delay = 30 + 10.0/(0.01 + percent);

	SetSpeed(percent);
}

PRIVATE void SetPosition(float position)
{
	char	string[40];

	static char	*axis[] = { "X", "Y", "Z" };

	(void) sprintf(string,"%s Position: %6.3f", axis[Common.axis],
			position);

	XtVaSetValues(PositionLabel, XtNlabel, string, NULL);
}
	
PRIVATE void NewPosition(Widget w, XtPointer client_data, XtPointer percent_ptr)
{
	float	percent;
	int	i;
	float	old_pos;

	percent = * (float *) percent_ptr;

	i = Common.axis;

	old_pos = Common.position;

	Common.position = 
		(Common.bounds[i].high - Common.bounds[i].low) * percent +
		Common.bounds[i].low;
				
	SetPosition(Common.position);

	if (Common.position == old_pos) return;

	if (normal) SceneUpdate();
}

PRIVATE void SetThickness(float thickness)
{
	char	string[40];

	(void) sprintf(string,"Thickness: %6.3f", thickness);

	XtVaSetValues(ThicknessLabel, XtNlabel, string, NULL);
}
	
PRIVATE void NewThickness(Widget w, XtPointer client_data, XtPointer percent_ptr)
{
	float	percent;
	int	i;
	float	old_thick;

	percent = * (float *) percent_ptr;
	i = Common.axis;

	old_thick = Common.thickness;
	Common.thickness = 
		(Common.bounds[i].high - Common.bounds[i].low) * percent;

	if (old_thick == Common.thickness) return;

	SetThickness(Common.thickness);

	if (normal) SceneUpdate();
}

/* **************************************************************** */

PRIVATE void AxisSelect(Widget w, XtPointer client_data, XtPointer call_data)
{
	int	new_axis;

	new_axis  = (int) client_data;

	if (new_axis == Common.axis) return;

	Common.axis = new_axis;

	DrawAxes();

	normal = 0;

	FindAndDisplaySlider(PositionBar, NewPosition);
	FindAndDisplaySlider(ThicknessBar, NewThickness);

	normal = 1;

	SceneUpdate();
}

PRIVATE void DirectionSelect(Widget w, XtPointer client_data, XtPointer call_data)
{
	int	new_dir;

	new_dir = (int) client_data;

	if (new_dir == Common.direction) return;

	Common.direction = new_dir;

	DrawAxes();
	SceneUpdate();
}

PRIVATE void MotionSelect(Widget w, XtPointer client_data, XtPointer call_data)
{
	int	new_dir;

	new_dir = (int) client_data;

	Common.dstep = new_dir;
}


PRIVATE void FindAndDisplaySlider(Widget w, XtCallbackProc func)
{
	float	percent;

	XtVaGetValues(w, XtNtopOfThumb, &percent, NULL);

	(func)(NULL, NULL, &percent);
}

void PositionUpdate(void)
{
	FindAndDisplaySlider(PositionBar, NewPosition);
}

void ThicknessUpdate(void)
{
	FindAndDisplaySlider(ThicknessBar, NewThickness);
}

void SpeedUpdate(void)
{
	FindAndDisplaySlider(SpeedBar, NewSpeed);
}

/* **************************************************************** */

typedef struct {
	int		x;
	int		y;
	char		*name;
	} LABEL;

typedef struct {
	XSegment	lines[6];
	LABEL		labels[2];
	} AXIS;

static AXIS	axis_data[] = {
		{ { { 120, 120, 120, 30 },	{ 115, 35, 120, 30 },
		    { 120, 30, 125, 35 },	{ 120, 120, 30, 120 },
		    { 35, 125, 30, 120 },	{ 30, 120, 35, 115 } },
		  { { 15, 127, "Z" }, 		{ 115, 25, "Y" } },
		},
		{ { { 20, 20, 20, 110 },	{ 15, 105, 20, 110 },
		    { 20, 110, 25, 105 },	{ 20, 20, 110, 20 },
		    { 105, 15, 110, 20 },	{ 110, 20, 105, 25 } },
		  { { 115, 27, "X" }, 		{ 15, 127,  "Z" } },
		},
		{ { { 20, 120, 20, 30 },	{ 15, 35, 20, 30 },
		    { 20, 30, 25, 35 },		{ 20, 120, 110, 120 },
		    { 105, 125, 110, 120 },	{ 110, 120, 105, 115 } },
		  { { 115, 127, "X" }, 		{ 15, 27,  "Y" } },
		},
		{ { { 20, 120, 20, 30 },	{ 15, 35, 20, 30 },
		    { 20, 30, 25, 35 },		{ 20, 120, 110, 120 },
		    { 105, 125, 110, 120 },	{ 110, 120, 105, 115 } },
		  { { 115, 127, "Z" }, 		{ 15, 27,  "Y" } },
		},
		{ { { 20, 120, 20, 30 },	{ 15, 35, 20, 30 },
		    { 20, 30, 25, 35 },		{ 20, 120, 110, 120 },
		    { 105, 125, 110, 120 },	{ 110, 120, 105, 115 } },
		  { { 115, 127, "X" }, 		{ 15, 27,  "Z" } },
		},
		{ { { 20, 20, 20, 110 },	{ 15, 105, 20, 110 },
		    { 20, 110, 25, 105 },	{ 20, 20, 110, 20 },
		    { 105, 15, 110, 20 },	{ 110, 20, 105, 25 } },
		  { { 115, 27, "X" }, 		{ 15, 127,  "Y" } },
		},
	};

PRIVATE void DrawAxes(void)
{
	static int		first = 1;
	static GC		gc;
	static Dimension	height, width;
	XGCValues		xgc;
	AXIS			*a;
	int			i;
	XFontStruct		*font;

	if (	Common.axis < 0 ||
		Common.direction < 0 ||
		!Axes ||
		!XtIsRealized(Axes) ) return;

	if (first) {
		xgc.function = GXcopy;

		XtVaGetValues(Axes,
				XtNbackground,	&xgc.background,
				XtNforeground,	&xgc.foreground, 
				XtNfont,	&font,
				XtNheight,	&height,
				XtNwidth,	&width,
				NULL );

		xgc.font = font->fid;

		gc = XCreateGC(XtDisplay(Axes), XtWindow(Axes),
			GCFunction | GCBackground | GCForeground | GCFont,
			&xgc);

		first = 0;
	}

	XClearWindow(XtDisplay(Axes), XtWindow(Axes));

	a = axis_data + Common.axis + 3*Common.direction;

	XDrawSegments(XtDisplay(Axes), XtWindow(Axes), gc, a->lines, 6);

	for(i = 0; i < 2; i++)
		XDrawString(XtDisplay(Axes), XtWindow(Axes), gc, 
			a->labels[i].x, a->labels[i].y,
			a->labels[i].name,
			strlen(a->labels[i].name));
		
}
void ExposeAxes(Widget w, XEvent *event, String *strings, Cardinal *nstrings)
{
	DrawAxes();
}


/* **************************************************************** */

PRIVATE void MakeRadioGroup(Widget parent, RadioStruct *radio)
{
	ToggleStruct	*t;
	int		i;

	if (radio->name)
		radio->lw = XtCreateManagedWidget(radio->name, 
			labelWidgetClass, parent , NULL, ZERO );

	for(t = radio->buttons, i = radio->nbuttons; i; i--, t++) {
		t->w = XtCreateManagedWidget(t->name, toggleWidgetClass,
			parent, NULL, ZERO);
		XtAddCallback(t->w, XtNcallback, radio->callback,
			(XtPointer) t->value);
	}

}

/* **************************************************************** */

PRIVATE void MakeSliderWithLabel(Widget parent, SliderLabelStruct *sl)
{
	/* create the slider */

	*(sl->slider) = XtCreateManagedWidget(sl->slidername, 
		scrollbarWidgetClass, parent, NULL, ZERO);

	/* add the callback */

	XtAddCallback(*(sl->slider), XtNjumpProc, sl->callback, NULL);

	/* create the label */

	*(sl->label) = XtCreateManagedWidget(sl->labelname,
		labelWidgetClass, parent, NULL, ZERO);

	/* Initialize the display */

	FindAndDisplaySlider(*(sl->slider), sl->callback);
}


void UpdateRadios(void)
{
	SetRadio(&AxisRadio);
	SetRadio(&DirectionRadio);
	SetRadio(&MotionRadio);
}

PRIVATE void SetRadio(RadioStruct *radio)
{
	ToggleStruct	*t;
	int		i;
	char		set_name[80], *rg;

	t = radio->buttons;

	/* find out which toggle is currently set */
	strcpy(set_name, (char *) XawToggleGetCurrent(t->w));

	for(i = radio->nbuttons; i; i--, t++)
		if (!strcmp(set_name, t->name)) goto found_set;

	t = radio->buttons;

	XtVaGetValues(t->w, XtNradioGroup, &rg, NULL);
	fprintf(stderr,"Error: radioData for radioGroup %s inconsistent.\n",
		rg);
        XtDestroyApplicationContext(XtWidgetToApplicationContext(TopLevel));
        exit(0);

	found_set:

	/* call the callback with the proper data */

	(radio->callback)((Widget) t->w, (XtPointer) t->value, (XtPointer) 0);
}

void SetReadString(char *s)
{
	if (ReadInfo == (Widget) 0) return;

	XtVaSetValues(ReadInfo, XtNlabel, s, NULL);
}

/* **************************************************************** */

PRIVATE void init_quit(Widget parent)
{
	Widget	dialog;

	/* create a popup shell */

	QuitPopup = XtCreatePopupShell("popup", POPUPSHELLCLASS,
			parent, NULL, ZERO);

	/* put a dialog box in it */

        dialog = XtCreateManagedWidget("dialog", dialogWidgetClass,
                        QuitPopup, NULL, ZERO);

        XawDialogAddButton(dialog, "ok", quit_ok, (XtPointer) dialog);
        XawDialogAddButton(dialog, "cancel", cancel_popup,
				(XtPointer) &QuitData);

}

/* **************************************************************** */

PRIVATE void do_popup(Widget w, XtPointer client_data, XtPointer call_data)
{
	Position	x, y;
	PopupData	*pdata;

	pdata = (PopupData *) client_data;

	/* stop motion for now */

	if (pdata->stop_motion) {
		in_motion = Common.motion;
		RemoveMotion();
	}

	/* translate to corner of calling widgets window */

	XtTranslateCoords(w, 0, 0, &x, &y);

        /* Set us to pop up there */

	XtVaSetValues(*(pdata->w), XtNx, x, XtNy, y, NULL);

	/* do the popping */

	XtPopup(*(pdata->w), pdata->grab);
}


PRIVATE void quit_ok(Widget w, XtPointer client_data, XtPointer call_data)
{
        XtDestroyApplicationContext(XtWidgetToApplicationContext(w));
        exit(0);
}

PRIVATE void cancel_popup(Widget w, XtPointer client_data, XtPointer call_data)
{
	PopupData	*pdata;

	pdata = (PopupData *) client_data;

	if (pdata->stop_motion && in_motion) InstallMotion();

        XtPopdown(*(pdata->w));
}

/* **************************************************************** */

PRIVATE void do_start(Widget w, XtPointer client_data, XtPointer call_data)
{
	InstallMotion();
}

PRIVATE void do_stop(Widget w, XtPointer client_data, XtPointer call_data)
{
	RemoveMotion();
}

PRIVATE void do_restart(Widget w, XtPointer client_data, XtPointer call_data)
{
	RemoveMotion();

	if (Common.dstep > 0)
		Common.step = 0;
	else
		Common.step = Common.ndata-1;

	InstallMotion();
}

PRIVATE void do_step(Widget w, XtPointer client_data, XtPointer call_data)
{
	RemoveMotion();

	Common.step = CoerceStep(Common.step + 1);

	SceneUpdate();
}

PRIVATE void do_back(Widget w, XtPointer client_data, XtPointer call_data)
{
	RemoveMotion();

	Common.step = CoerceStep(Common.step - 1);

	SceneUpdate();
}

PRIVATE void do_save_off(Widget w, XtPointer client_data, XtPointer call_data)
{
  Common.saveflag = 0;
}

PRIVATE void do_save_on(Widget w, XtPointer client_data, XtPointer call_data)
{
  Common.saveflag = 1;
}

PRIVATE void do_save_one(Widget w, XtPointer client_data, XtPointer call_data)
{
  SceneSave();
}

/* **************************************************************** */

PRIVATE void init_atoms(Widget parent)
{
	Widget	outerform, form;
	Widget	apply, cancel, last;
	char	formname[80], fglabel[80], sizelabel[80];
	char	fgname[80], sizename[80], vname[80], iname[80];
	int	i, j;
	size_t	size;

	/* create a popup shell for atom colors */

	AtomPopup = XtCreatePopupShell("atomcolors",
			transientShellWidgetClass, TopLevel, NULL, ZERO);

	/* put a form widget in it for layout */

	outerform = XtCreateManagedWidget("form", formWidgetClass, AtomPopup,
		NULL, ZERO);
	
	/* add the buttons */

	apply = XtCreateManagedWidget("apply", commandWidgetClass, 
		outerform, NULL, ZERO);
	XtAddCallback(apply, XtNcallback, color_apply, 
		(XtPointer) NULL);

	cancel = XtVaCreateManagedWidget("cancel", commandWidgetClass,
		outerform, XtNfromHoriz, apply, NULL);
	XtAddCallback(cancel, XtNcallback, cancel_popup, 
		(XtPointer) &AtomData);

	/* put the dialog boxes in it for the colors */

	/* ************************************************ */
	/* allocate all the widgets                         */
	/* ************************************************ */

	size = Common.natomcolors * sizeof(Widget);

	FGDialog 	= (Widget *) XtMalloc(size);
	SizeDialog	= (Widget *) XtMalloc(size);
	Visible		= (Widget *) XtMalloc(size);
	Invisible	= (Widget *) XtMalloc(size);

	/* ************************************************ */
	/* build the rest of the panel                      */
	/* ************************************************ */

	last = apply;

	for(i = 0; i < Common.natomcolors; i++) {
		j = i+1;

		sprintf(formname, 	"colorform%i",		j);
		sprintf(fgname,		"color%i",		j);
		sprintf(fglabel,	"Type %i Color",	j);
		sprintf(sizename,	"size%i",		j);
		sprintf(sizelabel,	"Type %i Diameter",	j);
		sprintf(vname,		"visible%i",		j);
		sprintf(iname,		"invisible%i",		j);

		form = XtVaCreateManagedWidget(formname, formWidgetClass,
			outerform, XtNfromVert, last, NULL);

		FGDialog[i] = XtVaCreateManagedWidget(fgname, 
			dialogWidgetClass, form, XtNlabel, fglabel, 
			XtNvalue, "white", NULL);

		SizeDialog[i] = XtVaCreateManagedWidget(sizename,
			dialogWidgetClass, form, XtNlabel, sizelabel,
			XtNvalue, "10", XtNfromHoriz, FGDialog[i], NULL);

		Visible[i] = XtVaCreateManagedWidget(vname, 
			toggleWidgetClass, form, XtNlabel, "Visible", 
			XtNstate, True, XtNfromVert, FGDialog[i], NULL);

		Invisible[i] = XtVaCreateManagedWidget(iname,	
			toggleWidgetClass, form, XtNlabel, "Invisible",
			XtNstate, False, XtNradioGroup, Visible[i],
			XtNfromHoriz, Visible[i], XtNfromVert, FGDialog[i],
			NULL);

		XtAddCallback(Visible[i], XtNcallback, visible_toggle, 
			(XtPointer) i);

		XtAddCallback(Invisible[i], XtNcallback, visible_toggle,
			(XtPointer) i);

		Common.atoms_visible[i] = True;
		Common.diameter[i] = 10;
		last = form;
	}
}

/* **************************************************************** */

PRIVATE void init_bg(Widget parent)
{
	Widget	outerform;
	Widget	apply, cancel;

	/* create a popup shell for background color */

	BackgroundPopup = XtCreatePopupShell("backcolors",
			transientShellWidgetClass, TopLevel, NULL, ZERO);

	/* put a form widget in it for layout */

	outerform = XtCreateManagedWidget("form", formWidgetClass,
		BackgroundPopup, NULL, ZERO);
	
	/* add the buttons */

	apply = XtCreateManagedWidget("apply", commandWidgetClass, 
		outerform, NULL, ZERO);
	XtAddCallback(apply, XtNcallback, bg_apply, (XtPointer) NULL);

	cancel = XtVaCreateManagedWidget("cancel", commandWidgetClass,
		outerform, XtNfromHoriz, apply, NULL);
	XtAddCallback(cancel, XtNcallback, cancel_popup, 
		(XtPointer) &BackgroundData);

	/* put the dialog box in it for the color */

	BGDialog = XtVaCreateManagedWidget("bg", dialogWidgetClass,
		outerform, XtNfromVert, apply, NULL);
}

/* **************************************************************** */

PRIVATE void init_file(Widget parent)
{
	Widget	outerform;
	Widget	apply, cancel;

	/* create a popup shell for color file */

	FileColorPopup = XtCreatePopupShell("filecolors",
			transientShellWidgetClass, TopLevel, NULL, ZERO);

	/* put a form widget in it for layout */

	outerform = XtCreateManagedWidget("form", formWidgetClass,
		FileColorPopup, NULL, ZERO);
	
	/* add the buttons */

	apply = XtCreateManagedWidget("apply", commandWidgetClass, 
		outerform, NULL, ZERO);
	XtAddCallback(apply, XtNcallback, file_apply, (XtPointer) NULL);

	cancel = XtVaCreateManagedWidget("cancel", commandWidgetClass,
		outerform, XtNfromHoriz, apply, NULL);
	XtAddCallback(cancel, XtNcallback, cancel_popup, 
		(XtPointer) &FileColorData);

	/* put the dialog box in it for the filename */

	FileDialog = XtVaCreateManagedWidget("file", dialogWidgetClass,
		outerform, XtNfromVert, apply, NULL);
}

PRIVATE void init_bonds(Widget parent)
{
	Widget	outerform, form;
	Widget	apply, cancel, last;
	char	formname[80], fglabel[80], sizelabel[80];
	char	fgname[80], sizename[80], vname[80], iname[80];
	int	i, j;
	size_t	size;

	/* create a popup shell for bond colors */

	BondPopup = XtCreatePopupShell("bondcolors",
			transientShellWidgetClass,
			TopLevel, NULL, ZERO);

	/* put a form widget in it for layout */

	outerform = XtCreateManagedWidget("form", formWidgetClass, BondPopup,
		NULL, ZERO);
	
	/* add the buttons */

	apply = XtCreateManagedWidget("apply", commandWidgetClass, 
		outerform, NULL, ZERO);
	XtAddCallback(apply, XtNcallback, bond_apply, 
		(XtPointer) NULL);

	cancel = XtVaCreateManagedWidget("cancel", commandWidgetClass,
		outerform, XtNfromHoriz, apply, NULL);
	XtAddCallback(cancel, XtNcallback, cancel_popup, 
		(XtPointer) &BondData);

	/* ************************************************ */
	/* allocate all the widgets                         */
	/* ************************************************ */

	size = Common.nbondcolors * sizeof(Widget);

	BFGDialog 	= (Widget *) XtMalloc(size);
	ThicknessDialog	= (Widget *) XtMalloc(size);
	BVisible	= (Widget *) XtMalloc(size);
	BInvisible	= (Widget *) XtMalloc(size);

	/* ************************************************ */
	/* build the rest of the panel                      */
	/* ************************************************ */

	last = apply;

	for(i = 0; i < Common.nbondcolors; i++) {
		j = i+1;

		sprintf(formname, 	"bondform%i",		j);
		sprintf(fgname,		"color%i",		j);
		sprintf(fglabel,	"Bond Type %i Color",	j);
		sprintf(sizename,	"thickness%i",		j);
		sprintf(sizelabel,	"Bond Type %i Thickness",j);
		sprintf(vname,		"visible%i",		j);
		sprintf(iname,		"invisible%i",		j);

		form = XtVaCreateManagedWidget(formname, 
			formWidgetClass, outerform, XtNfromVert, 
			last, NULL);

		BFGDialog[i] = XtVaCreateManagedWidget(fgname, 
			dialogWidgetClass, form, XtNlabel, fglabel, 
			XtNvalue, "white", NULL);

		ThicknessDialog[i] = XtVaCreateManagedWidget(sizename,
			dialogWidgetClass, form, XtNlabel, sizelabel,
			XtNvalue, "1", XtNfromHoriz, BFGDialog[i], NULL);

		BVisible[i] = XtVaCreateManagedWidget(vname, 
			toggleWidgetClass, form, XtNlabel, "Visible", 
			XtNstate, True, XtNfromVert, BFGDialog[i], NULL);

		BInvisible[i] = XtVaCreateManagedWidget(iname,	
			toggleWidgetClass, form, XtNlabel, "Invisible",
			XtNstate, False, XtNradioGroup, BVisible[i],
			XtNfromHoriz, BVisible[i], XtNfromVert, BFGDialog[i],
			NULL);

		XtAddCallback(BVisible[i], XtNcallback, bvisible_toggle, 
			(XtPointer) i);

		XtAddCallback(BInvisible[i], XtNcallback, bvisible_toggle,
			(XtPointer) i);

		Common.bonds_visible[i] = True;
		last = form;
	}
}

/* **************************************************************** */

PRIVATE void color_apply(Widget w, XtPointer client_data, XtPointer call_data)
{
	static String	*fg, *size;
	static Pixel	*fgpixel;
	static int	*diam;
	static Bool	first = True;

	Boolean		ok;
	int		i, n;


	if (first) {
		n = Common.natomcolors;

		fg 	= (String *)	XtMalloc( n * sizeof(String) );
		size	= (String *)	XtMalloc( n * sizeof(String) );
		fgpixel	= (Pixel *)	XtMalloc( n * sizeof(Pixel) );
		diam	= (int *)	XtMalloc( n * sizeof(int) );

		first = False;
	}
	
	ok = True;

	for(i = 0; i < Common.natomcolors; i++) {
		fg[i] = XawDialogGetValueString(FGDialog[i]);
		fgpixel[i] = ConvertColorToPixel(FGDialog[i], fg[i]);

		if (fgpixel[i] == (Pixel) -1) {
			ok = False;
			XtVaSetValues(FGDialog[i], XtNvalue,
				"-Invalid Color-", NULL);
		}

		size[i] = XawDialogGetValueString(SizeDialog[i]);
		diam[i] = atoi(size[i]);

		if (diam[i] < 1 || diam[i] > MAXDIAM) {
			ok = False;
			XtVaSetValues(SizeDialog[i], XtNvalue,
				"-Invalid Size-", NULL);
		}
	}

	if (!ok){
		XBell(XtDisplay(BGDialog), 0);
		return;
	}

	SetAtomColors(fgpixel);

	for(i = 0; i < Common.natomcolors; i++)
		Common.diameter[i] = diam[i];
	
	SceneUpdate();
}

/* **************************************************************** */

PRIVATE void bond_apply(Widget w, XtPointer client_data, XtPointer call_data)
{
	static String	*fg, *size;
	static Pixel	*fgpixel;
	static Dimension *thick;
	static Bool	first = True;

	Boolean		ok;
	int		i, n, ithick;

	if (first) {
		n = Common.nbondcolors;

		fg 	= (String *)	XtMalloc( n * sizeof(String) );
		size	= (String *)	XtMalloc( n * sizeof(String) );
		fgpixel	= (Pixel *)	XtMalloc( n * sizeof(Pixel) );
		thick	= (Dimension *)	XtMalloc( n * sizeof(Dimension) );

		first = False;
	}
	
	ok = True;

	for(i = 0; i < Common.nbondcolors; i++) {
		fg[i] = XawDialogGetValueString(BFGDialog[i]);
		fgpixel[i] = ConvertColorToPixel(BFGDialog[i], fg[i]);

		if (fgpixel[i] == (Pixel) -1) {
			ok = False;
			XtVaSetValues(BFGDialog[i], XtNvalue,
				"-Invalid Color-", NULL);
		}

		size[i] = XawDialogGetValueString(ThicknessDialog[i]);
		ithick = atoi(size[i]);

		if (ithick < 1 || ithick > MAXTHICK) {
			ok = False;
			XtVaSetValues(ThicknessDialog[i], XtNvalue,
				"-Invalid Size-", NULL);
		}
		thick[i] = ithick;		
	}

	if (!ok){
		XBell(XtDisplay(BGDialog), 0);
		return;
	}

	SetBondColors(fgpixel, thick);

	SceneUpdate();
}

/* **************************************************************** */

PRIVATE void bg_apply(Widget w, XtPointer client_data, XtPointer call_data)
{
	String		bg;
	Pixel		bgpixel;

	bg = XawDialogGetValueString(BGDialog);
	bgpixel = ConvertColorToPixel(BGDialog, bg);

	if (bgpixel == (Pixel) -1) {
		XtVaSetValues(BGDialog, XtNvalue, "-Invalid Color-", NULL);
		XBell(XtDisplay(BGDialog), 0);
		return;
	}

	SetBGColor(bgpixel);

	SceneUpdate();
}

PRIVATE Pixel ConvertColorToPixel(Widget w, String color)
{
	XrmValue	from, to;
	Pixel		pixel;

	from.size = strlen(color) + 1;
	from.addr = color;
	/* SJP - 10/02, new Linux couldn't fine caddr_t */
	/*	to.addr = (caddr_t) &pixel; */
	to.addr = (char *) &pixel;
	to.size = sizeof(pixel);

	if (XtConvertAndStore(w, XtRString, (XrmValuePtr) &from, 
			XtRPixel,  (XrmValuePtr) &to))
		return(pixel);
	else
		return( (Pixel) -1);
}


PRIVATE void file_apply(Widget w, XtPointer client_data, XtPointer call_data)
{
	String		bg;
	Pixel		bgpixel;
	static String	*fg, *size;
	static Pixel	*fgpixel;
	static int	*diam;
	static Bool	first = True;

	int		i, n;

	if (first) {
		n = Common.natomcolors;

		fg 	= (String *)	XtMalloc( n * sizeof(String) );
		size	= (String *)	XtMalloc( n * sizeof(String) );
		fgpixel	= (Pixel *)	XtMalloc( n * sizeof(Pixel) );
		diam	= (int *)	XtMalloc( n * sizeof(int) );

		first = False;
	}
	
	for(i = 0; i < Common.natomcolors; i++) {
	  fg[i] = "red";
	  fgpixel[i] = ConvertColorToPixel(FGDialog[i], fg[i]);
	  
	  if (fgpixel[i] == (Pixel) -1) return;

	  XtVaSetValues(FGDialog[i], XtNvalue, fg[i], NULL);

	  size[i] = "20";
	  diam[i] = atoi(size[i]);

	  if (diam[i] < 1 || diam[i] > MAXDIAM)  return;
	  XtVaSetValues(SizeDialog[i], XtNvalue, size[i], NULL);
	}

	SetAtomColors(fgpixel);

	for(i = 0; i < Common.natomcolors; i++)
		Common.diameter[i] = diam[i];

	/*
	Common.atoms_visible[0] = 0;
	Common.atoms_visible[1] = 0;
	Common.atoms_visible[2] = 0;
	Common.atoms_visible[3] = 0;
	*/

	bg = "green";
	bgpixel = ConvertColorToPixel(BGDialog, bg);

	if (bgpixel == (Pixel) -1) {
		XtVaSetValues(BGDialog, XtNvalue, "-Invalid Color-", NULL);
		XBell(XtDisplay(BGDialog), 0);
		return;
	}

	SetBGColor(bgpixel);

	SceneUpdate();
}

PRIVATE void visible_toggle(Widget w, XtPointer client_data, 
	XtPointer call_data)
{
	int		i;
	XtPointer	data;
	char		c;

	i = (int) client_data;

	XtVaGetValues(w, XtNradioData, &data, NULL);

	c = *((char *) data);

	Common.atoms_visible[i] = c == 'v';
}

PRIVATE void bvisible_toggle(Widget w, XtPointer client_data, 
	XtPointer call_data)
{
	int		i;
	XtPointer	data;
	char		c;

	i = (int) client_data;

	XtVaGetValues(w, XtNradioData, &data, NULL);

	c = *((char *) data);

	Common.bonds_visible[i] = c == 'v';
}
