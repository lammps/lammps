#include <stdio.h>
#include <stdlib.h>

#include <X11/StringDefs.h>
#include <X11/Intrinsic.h>

#include "xmovie.h"

void Version(void)
{
	fprintf(stderr,"%s %s, compiled %s %s\n"
			"Copyright 1992 Michael Uttormark,  All rights reserved.\n",
			Progname, VERSION, __DATE__, __TIME__);
	exit(EXIT_SUCCESS);
}
