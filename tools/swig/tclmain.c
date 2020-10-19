#include <tcl.h>
int AppInit(Tcl_Interp *interp) {
     if (Tcl_Init(interp) == TCL_ERROR) return TCL_ERROR;
     Tcl_SetVar(interp,"tcl_rcFileName","~/.wishrc",TCL_GLOBAL_ONLY);
     return TCL_OK;
}

int main(int argc, char *argv[]) {
    Tcl_Main(argc, argv, AppInit);
    return 0;
}
