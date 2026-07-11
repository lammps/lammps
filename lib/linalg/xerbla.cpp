
#include "lmp_f2c.h"

#undef abs
#include <cstdio>
#include <cstdlib>
#include <cstring>

extern "C" {

static constexpr int BUFSZ = 1024;

integer xerbla_(const char *srname, integer *info)
{
    char buf[BUFSZ];
    buf[0] = '\0';

    strcat(buf, " ** On entry to ");
    for (int i = 0; i < BUFSZ - 16; ++i) {
        if ((srname[i] == '\0') || (srname[i] == ' ')) {
            buf[i + 16] = '\0';
            break;
        }
        buf[i + 16] = srname[i];
    }
    int len = strlen(buf);
    snprintf(buf + len, BUFSZ - len, " parameter number %d had an illegal value\n", *info);
    exit(1);
    return 0;
}
}
