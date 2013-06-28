/* **************************************************************** */
/* Resource database values */
/* Can be overridden from .Xdefaults */

String FallbackResources[] = {
	"*.foreground:				black",
	"*.background:				white",

	"*.horizDistance:			4",
	"*.vertDistance:			4",

	"*Dialog*Translations:			#override \\n\
		<Key>Return:			NoOp()\\n\
		<Key>Up:			NoOp()\\n\
		<Key>Down:			NoOp()\\n\
		Ctrl<Key>C:			NoOp()\\n\
		Ctrl<Key>I:			NoOp()\\n\
		Ctrl<Key>J:			NoOp()\\n\
		Ctrl<Key>M:			NoOp()\\n\
		Ctrl<Key>N:			NoOp()\\n\
		Ctrl<Key>O:			NoOp()\\n\
		Ctrl<Key>P:			NoOp()\\n\
		Ctrl<Key>Q:			NoOp()\\n\
		Ctrl<Key>R:			NoOp()\\n\
		Ctrl<Key>S:			NoOp()\\n\
		Ctrl<Key>V:			NoOp()\\n\
		Ctrl<Key>X:			NoOp()\\n\
		Ctrl<Key>Z:			NoOp()\\n\
		Meta<Key>V:			NoOp()",

	"*.Toggle.Translations:	#replace \\n\
		<EnterWindow>:			highlight(Always) \\n\
		<LeaveWindow>:			unhighlight() \\n\
		<Btn1Down>,<Btn1Up>:		set()notify() \\n\
		<Btn2Down>,<Btn2Up>:		set()notify() \\n\
		<Btn3Down>,<Btn3Up>:		set()notify()",

	"*.quit.label:				Quit",
	"*.start.label:				Start",
	"*.stop.label:				Stop",
	"*.restart.label:			Restart",
	"*.step.label:				+Step",
	"*.back.label:				-Step",
	"*.save.label:				Save",
	"*.color.label:				Color",

	"*.start.fromHoriz:			quit",
	"*.stop.fromHoriz:			start",
	"*.restart.fromHoriz:			stop",
	"*.step.fromHoriz:			restart",
	"*.back.fromHoriz:			step",
	"*.save.fromHoriz:			back",
	"*.color.fromHoriz:			save",

	"*.quit.popup.dialog.label:		Do you really want to Quit?",
	"*.quit.popup.dialog.ok.label:		Ok",
	"*.quit.popup.dialog.cancel.label:	Cancel",
	"*.quit.popup.overrideRedirect:		True",

	"*.save.menu.atoms.label:		Animaton Off",
	"*.save.menu.bonds.label:		Animation On",
	"*.save.menu.background.label:		This Snapshot",
	      
	"*.color.menu.atoms.label:		Atoms",
	"*.color.menu.bonds.label:		Bonds",
	"*.color.menu.background.label:		Background",
	"*.color.menu.filecolor.label:		from File",

	"*.atomcolors.title:			xmovie atom colors",

	"*.atomcolors.*.apply.label:		Apply",
	"*.atomcolors.*.cancel.label:		Dismiss",

	"*.atomcolors.*.Form.Dialog.borderWidth:	2",
	"*.atomcolors.*.Dialog.*.resizable:	True",
	"*.atomcolors.*.Dialog.Text.width:	200",
	"*.atomcolors.*.Form.Form.Dialog.borderWidth:	0",
	"*.atomcolors.*.Form.Form.borderWidth:	2",

	"*.bondcolors.title:			xmovie bond colors",

	"*.bondcolors.*.apply.label:		Apply",
	"*.bondcolors.*.cancel.label:		Dismiss",

	"*.bondcolors.*.Form.Dialog.borderWidth:	2",
	"*.bondcolors.*.Dialog.*.resizable:	True",
	"*.bondcolors.*.Dialog.Text.width:	200",
	"*.bondcolors.*.Form.Form.Dialog.borderWidth:	0",
	"*.bondcolors.*.Form.Form.borderWidth:	2",

	"*.backcolors.title:			xmovie background colors",

	"*.backcolors.*.apply.label:		Apply",
	"*.backcolors.*.cancel.label:		Dismiss",

	"*.backcolors.*.bg.label:		Background Color",
	"*.backcolors.*.bg.value:		black",
	"*.backcolors.*.Dialog.Text.width:	200",
	"*.backcolors.*.Dialog.*.resizable:	True",

	"*.filecolors.title:			xmovie colors from file",

	"*.filecolors.*.apply.label:		Apply",
	"*.filecolors.*.cancel.label:		Dismiss",

	"*.filecolors.*.file.label:		File of color settings",
	"*.filecolors.*.file.value:		xmovie.colors",
	"*.filecolors.*.Dialog.Text.width:	200",
	"*.filecolors.*.Dialog.*.resizable:	True",

	"*.sceneshell.title:			xmovie scene",
	"*.scene.foreground:			white",
	"*.scene.background:			black",
	"*.scene.width:				400",
	"*.scene.height:			400",
	"*.scene.Translations:			#override \\n\
		<Expose>:			ExposeScene()",

	"*.Scrollbar.orientation:		horizontal",
	"*.Scrollbar.height:			25",
	"*.Scrollbar.width:			200",
	"*.Scrollbar.shown:			0.06",

	"*.Label.resize:			False",

	"*.time.label:				Time: ?",
	"*.time.width:				406",
	"*.time.fromVert:			quit",

	"*.speedbar.fromVert:			time",
	"*.speedbar.topOfThumb:			0.25",

	"*.speed.fromHoriz:			speedbar",
	"*.speed.fromVert:			time",
	"*.speed.label:				Speed: ?",
	"*.speed.width:				200",
	"*.speed.height:			25",
	"*.speed.resize:			False",

	"*.positionbar.fromVert:		speedbar",
	"*.positionbar.topOfThumb:		0.5",

	"*.position.fromHoriz:			positionbar",
	"*.position.fromVert:			speed",
	"*.position.label:			Position: ?",
	"*.position.width:			200",
	"*.position.height:			25",
	"*.position.resize:			False",

	"*.thicknessbar.fromVert:		positionbar",
	"*.thicknessbar.topOfThumb:		0.1",

	"*.thickness.fromHoriz:			thicknessbar",
	"*.thickness.fromVert:			position",
	"*.thickness.label:			Speed: ?",
	"*.thickness.width:			200",
	"*.thickness.height:			25",
	"*.thickness.resize:			False",

	"*.Scrollbar.Translations:	#override \\n\
	<Btn1Down>:	StartScroll(Continuous)MoveThumb()NotifyThumb() \\n\
	<Btn3Down>:	StartScroll(Continuous)MoveThumb()NotifyThumb() \\n\
	<Btn1Motion>:	MoveThumb()NotifyThumb() \\n\
	<Btn3Motion>:	MoveThumb()NotifyThumb()",

	"*.axis.label:				Viewing Axis",
	"*.axis.fromVert:			thicknessbar",
	"*.axis.width:				200",
	"*.axis.height:				25",

	"*.x.fromVert:				axis",
	"*.x.label:				X",
	"*.x.width:				26",
	"*.x.height:				25",
	"*.x.horizDistance:			61",

	"*.y.fromVert:				axis",
	"*.y.fromHoriz:				x",
	"*.y.label:				Y",
	"*.y.width:				26",
	"*.y.height:				25",
	"*.y.radioGroup:			x",
		
	"*.z.fromVert:				axis",
	"*.z.fromHoriz:				y",
	"*.z.label:				Z",
	"*.z.width:				26",
	"*.z.height:				25",
	"*.z.radioGroup:			x",
	"*.z.state:				True",
	
	"*.viewdir.label:			Viewing Direction",
	"*.viewdir.fromVert:			x",
	"*.viewdir.width:			200",
	"*.viewdir.height:			25",

	"*.vplus.fromVert:			viewdir",
	"*.vplus.horizDistance:			76",
	"*.vplus.label:				+",
	"*.vplus.height:			25",
	"*.vplus.width:				26",

	"*.vminus.fromVert:			viewdir",
	"*.vminus.fromHoriz:			vplus",
	"*.vminus.label:			-",
	"*.vminus.width:			26",
	"*.vminus.height:			25",
	"*.vminus.radioGroup:			vplus",
	"*.vminus.state:			True",
	
	"*.motiondir.label:			Movie Direction",
	"*.motiondir.fromVert:			vplus",
	"*.motiondir.width:			200",
	"*.motiondir.height:			25",

	"*.mplus.fromVert:			motiondir",
	"*.mplus.horizDistance:			12",
	"*.mplus.label:				Forward",
	"*.mplus.width:				90",
	"*.mplus.height:			25",
	"*.mplus.state:				True",

	"*.mminus.fromVert:			motiondir",
	"*.mminus.fromHoriz:			mplus",
	"*.mminus.label:			Backward",
	"*.mminus.width:			90",
	"*.mminus.height:			25",
	"*.mminus.radioGroup:			mplus",
	
	"*.axes.fromVert:			thickness",
	"*.axes.fromHoriz:			axis",
	"*.axes.width:				200",
	"*.axes.height:				177",
	"*.axes.resize:				False",
	"*.axes.Translations:		#override \\n\
		<Expose>:			ExposeAxes()",

	"*.read.fromVert:			mplus",
	"*.read.height:				100",
	"*.read.width:				406",
	"*.read.resize:				False",

	NULL,
	};
