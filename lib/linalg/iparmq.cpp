#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
integer iparmq_(integer *ispec, char *name__, char *opts, integer *n, integer *ilo, integer *ihi,
                integer *lwork, ftnlen name_len, ftnlen opts_len)
{
    integer ret_val, i__1, i__2;
    real r__1;
    double log(doublereal);
    integer i_lmp_nint(real *);
    int s_lmp_copy(char *, char *, ftnlen, ftnlen);
    integer s_lmp_cmp(char *, char *, ftnlen, ftnlen);
    integer i__, ic, nh, ns, iz;
    char subnam[6];
    if (*ispec == 15 || *ispec == 13 || *ispec == 16) {
        nh = *ihi - *ilo + 1;
        ns = 2;
        if (nh >= 30) {
            ns = 4;
        }
        if (nh >= 60) {
            ns = 10;
        }
        if (nh >= 150) {
            r__1 = log((real)nh) / log((float)2.);
            i__1 = 10, i__2 = nh / i_lmp_nint(&r__1);
            ns = max(i__1, i__2);
        }
        if (nh >= 590) {
            ns = 64;
        }
        if (nh >= 3000) {
            ns = 128;
        }
        if (nh >= 6000) {
            ns = 256;
        }
        i__1 = 2, i__2 = ns - ns % 2;
        ns = max(i__1, i__2);
    }
    if (*ispec == 12) {
        ret_val = 75;
    } else if (*ispec == 14) {
        ret_val = 14;
    } else if (*ispec == 15) {
        ret_val = ns;
    } else if (*ispec == 13) {
        if (nh <= 500) {
            ret_val = ns;
        } else {
            ret_val = ns * 3 / 2;
        }
    } else if (*ispec == 16) {
        ret_val = 0;
        s_lmp_copy(subnam, name__, (ftnlen)6, name_len);
        ic = *(unsigned char *)subnam;
        iz = 'Z';
        if (iz == 90 || iz == 122) {
            if (ic >= 97 && ic <= 122) {
                *(unsigned char *)subnam = (char)(ic - 32);
                for (i__ = 2; i__ <= 6; ++i__) {
                    ic = *(unsigned char *)&subnam[i__ - 1];
                    if (ic >= 97 && ic <= 122) {
                        *(unsigned char *)&subnam[i__ - 1] = (char)(ic - 32);
                    }
                }
            }
        } else if (iz == 233 || iz == 169) {
            if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 162 && ic <= 169) {
                *(unsigned char *)subnam = (char)(ic + 64);
                for (i__ = 2; i__ <= 6; ++i__) {
                    ic = *(unsigned char *)&subnam[i__ - 1];
                    if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 ||
                        ic >= 162 && ic <= 169) {
                        *(unsigned char *)&subnam[i__ - 1] = (char)(ic + 64);
                    }
                }
            }
        } else if (iz == 218 || iz == 250) {
            if (ic >= 225 && ic <= 250) {
                *(unsigned char *)subnam = (char)(ic - 32);
                for (i__ = 2; i__ <= 6; ++i__) {
                    ic = *(unsigned char *)&subnam[i__ - 1];
                    if (ic >= 225 && ic <= 250) {
                        *(unsigned char *)&subnam[i__ - 1] = (char)(ic - 32);
                    }
                }
            }
        }
        if (s_lmp_cmp(subnam + 1, (char *)"GGHRD", (ftnlen)5, (ftnlen)5) == 0 ||
            s_lmp_cmp(subnam + 1, (char *)"GGHD3", (ftnlen)5, (ftnlen)5) == 0) {
            ret_val = 1;
            if (nh >= 14) {
                ret_val = 2;
            }
        } else if (s_lmp_cmp(subnam + 3, (char *)"EXC", (ftnlen)3, (ftnlen)3) == 0) {
            if (nh >= 14) {
                ret_val = 1;
            }
            if (nh >= 14) {
                ret_val = 2;
            }
        } else if (s_lmp_cmp(subnam + 1, (char *)"HSEQR", (ftnlen)5, (ftnlen)5) == 0 ||
                   s_lmp_cmp(subnam + 1, (char *)"LAQR", (ftnlen)4, (ftnlen)4) == 0) {
            if (ns >= 14) {
                ret_val = 1;
            }
            if (ns >= 14) {
                ret_val = 2;
            }
        }
    } else if (*ispec == 17) {
        ret_val = 10;
    } else {
        ret_val = -1;
    }
    return ret_val;
}
#ifdef __cplusplus
}
#endif
