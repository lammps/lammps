#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static real c_b176 = (float)0.;
static real c_b177 = (float)1.;
static integer c__0 = 0;
integer ilaenv_(integer *ispec, char *name__, char *opts, integer *n1, integer *n2, integer *n3,
                integer *n4, ftnlen name_len, ftnlen opts_len)
{
    integer ret_val, i__1, i__2, i__3;
    int s_lmp_copy(char *, char *, ftnlen, ftnlen);
    integer i_lmp_len(char *, ftnlen), s_lmp_cmp(char *, char *, ftnlen, ftnlen);
    logical twostage;
    integer i__;
    char c1[1], c2[2], c3[3], c4[2];
    integer ic, nb, iz, nx;
    logical cname;
    integer nbmin;
    logical sname;
    extern integer ieeeck_(integer *, real *, real *);
    char subnam[16];
    extern integer iparmq_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    switch (*ispec) {
        case 1:
            goto L10;
        case 2:
            goto L10;
        case 3:
            goto L10;
        case 4:
            goto L80;
        case 5:
            goto L90;
        case 6:
            goto L100;
        case 7:
            goto L110;
        case 8:
            goto L120;
        case 9:
            goto L130;
        case 10:
            goto L140;
        case 11:
            goto L150;
        case 12:
            goto L160;
        case 13:
            goto L160;
        case 14:
            goto L160;
        case 15:
            goto L160;
        case 16:
            goto L160;
        case 17:
            goto L160;
    }
    ret_val = -1;
    return ret_val;
L10:
    ret_val = 1;
    s_lmp_copy(subnam, name__, (ftnlen)16, name_len);
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
                if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 162 && ic <= 169) {
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
    *(unsigned char *)c1 = *(unsigned char *)subnam;
    sname = *(unsigned char *)c1 == 'S' || *(unsigned char *)c1 == 'D';
    cname = *(unsigned char *)c1 == 'C' || *(unsigned char *)c1 == 'Z';
    if (!(cname || sname)) {
        return ret_val;
    }
    s_lmp_copy(c2, subnam + 1, (ftnlen)2, (ftnlen)2);
    s_lmp_copy(c3, subnam + 3, (ftnlen)3, (ftnlen)3);
    s_lmp_copy(c4, c3 + 1, (ftnlen)2, (ftnlen)2);
    twostage = i_lmp_len(subnam, (ftnlen)16) >= 11 && *(unsigned char *)&subnam[10] == '2';
    switch (*ispec) {
        case 1:
            goto L50;
        case 2:
            goto L60;
        case 3:
            goto L70;
    }
L50:
    nb = 1;
    if (s_lmp_cmp(subnam + 1, (char *)"LAORH", (ftnlen)5, (ftnlen)5) == 0) {
        if (sname) {
            nb = 32;
        } else {
            nb = 32;
        }
    } else if (s_lmp_cmp(c2, (char *)"GE", (ftnlen)2, (ftnlen)2) == 0) {
        if (s_lmp_cmp(c3, (char *)"TRF", (ftnlen)3, (ftnlen)3) == 0) {
            if (sname) {
                nb = 64;
            } else {
                nb = 64;
            }
        } else if (s_lmp_cmp(c3, (char *)"QRF", (ftnlen)3, (ftnlen)3) == 0 ||
                   s_lmp_cmp(c3, (char *)"RQF", (ftnlen)3, (ftnlen)3) == 0 ||
                   s_lmp_cmp(c3, (char *)"LQF", (ftnlen)3, (ftnlen)3) == 0 ||
                   s_lmp_cmp(c3, (char *)"QLF", (ftnlen)3, (ftnlen)3) == 0) {
            if (sname) {
                nb = 32;
            } else {
                nb = 32;
            }
        } else if (s_lmp_cmp(c3, (char *)"QR ", (ftnlen)3, (ftnlen)3) == 0) {
            if (*n3 == 1) {
                if (sname) {
                    if (*n1 * *n2 <= 131072 || *n1 <= 8192) {
                        nb = *n1;
                    } else {
                        nb = 32768 / *n2;
                    }
                } else {
                    if (*n1 * *n2 <= 131072 || *n1 <= 8192) {
                        nb = *n1;
                    } else {
                        nb = 32768 / *n2;
                    }
                }
            } else {
                if (sname) {
                    nb = 1;
                } else {
                    nb = 1;
                }
            }
        } else if (s_lmp_cmp(c3, (char *)"LQ ", (ftnlen)3, (ftnlen)3) == 0) {
            if (*n3 == 2) {
                if (sname) {
                    if (*n1 * *n2 <= 131072 || *n1 <= 8192) {
                        nb = *n1;
                    } else {
                        nb = 32768 / *n2;
                    }
                } else {
                    if (*n1 * *n2 <= 131072 || *n1 <= 8192) {
                        nb = *n1;
                    } else {
                        nb = 32768 / *n2;
                    }
                }
            } else {
                if (sname) {
                    nb = 1;
                } else {
                    nb = 1;
                }
            }
        } else if (s_lmp_cmp(c3, (char *)"HRD", (ftnlen)3, (ftnlen)3) == 0) {
            if (sname) {
                nb = 32;
            } else {
                nb = 32;
            }
        } else if (s_lmp_cmp(c3, (char *)"BRD", (ftnlen)3, (ftnlen)3) == 0) {
            if (sname) {
                nb = 32;
            } else {
                nb = 32;
            }
        } else if (s_lmp_cmp(c3, (char *)"TRI", (ftnlen)3, (ftnlen)3) == 0) {
            if (sname) {
                nb = 64;
            } else {
                nb = 64;
            }
        }
    } else if (s_lmp_cmp(c2, (char *)"PO", (ftnlen)2, (ftnlen)2) == 0) {
        if (s_lmp_cmp(c3, (char *)"TRF", (ftnlen)3, (ftnlen)3) == 0) {
            if (sname) {
                nb = 64;
            } else {
                nb = 64;
            }
        }
    } else if (s_lmp_cmp(c2, (char *)"SY", (ftnlen)2, (ftnlen)2) == 0) {
        if (s_lmp_cmp(c3, (char *)"TRF", (ftnlen)3, (ftnlen)3) == 0) {
            if (sname) {
                if (twostage) {
                    nb = 192;
                } else {
                    nb = 64;
                }
            } else {
                if (twostage) {
                    nb = 192;
                } else {
                    nb = 64;
                }
            }
        } else if (sname && s_lmp_cmp(c3, (char *)"TRD", (ftnlen)3, (ftnlen)3) == 0) {
            nb = 32;
        } else if (sname && s_lmp_cmp(c3, (char *)"GST", (ftnlen)3, (ftnlen)3) == 0) {
            nb = 64;
        }
    } else if (cname && s_lmp_cmp(c2, (char *)"HE", (ftnlen)2, (ftnlen)2) == 0) {
        if (s_lmp_cmp(c3, (char *)"TRF", (ftnlen)3, (ftnlen)3) == 0) {
            if (twostage) {
                nb = 192;
            } else {
                nb = 64;
            }
        } else if (s_lmp_cmp(c3, (char *)"TRD", (ftnlen)3, (ftnlen)3) == 0) {
            nb = 32;
        } else if (s_lmp_cmp(c3, (char *)"GST", (ftnlen)3, (ftnlen)3) == 0) {
            nb = 64;
        }
    } else if (sname && s_lmp_cmp(c2, (char *)"OR", (ftnlen)2, (ftnlen)2) == 0) {
        if (*(unsigned char *)c3 == 'G') {
            if (s_lmp_cmp(c4, (char *)"QR", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"RQ", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"LQ", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"QL", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"HR", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"TR", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"BR", (ftnlen)2, (ftnlen)2) == 0) {
                nb = 32;
            }
        } else if (*(unsigned char *)c3 == 'M') {
            if (s_lmp_cmp(c4, (char *)"QR", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"RQ", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"LQ", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"QL", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"HR", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"TR", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"BR", (ftnlen)2, (ftnlen)2) == 0) {
                nb = 32;
            }
        }
    } else if (cname && s_lmp_cmp(c2, (char *)"UN", (ftnlen)2, (ftnlen)2) == 0) {
        if (*(unsigned char *)c3 == 'G') {
            if (s_lmp_cmp(c4, (char *)"QR", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"RQ", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"LQ", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"QL", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"HR", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"TR", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"BR", (ftnlen)2, (ftnlen)2) == 0) {
                nb = 32;
            }
        } else if (*(unsigned char *)c3 == 'M') {
            if (s_lmp_cmp(c4, (char *)"QR", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"RQ", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"LQ", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"QL", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"HR", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"TR", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"BR", (ftnlen)2, (ftnlen)2) == 0) {
                nb = 32;
            }
        }
    } else if (s_lmp_cmp(c2, (char *)"GB", (ftnlen)2, (ftnlen)2) == 0) {
        if (s_lmp_cmp(c3, (char *)"TRF", (ftnlen)3, (ftnlen)3) == 0) {
            if (sname) {
                if (*n4 <= 64) {
                    nb = 1;
                } else {
                    nb = 32;
                }
            } else {
                if (*n4 <= 64) {
                    nb = 1;
                } else {
                    nb = 32;
                }
            }
        }
    } else if (s_lmp_cmp(c2, (char *)"PB", (ftnlen)2, (ftnlen)2) == 0) {
        if (s_lmp_cmp(c3, (char *)"TRF", (ftnlen)3, (ftnlen)3) == 0) {
            if (sname) {
                if (*n2 <= 64) {
                    nb = 1;
                } else {
                    nb = 32;
                }
            } else {
                if (*n2 <= 64) {
                    nb = 1;
                } else {
                    nb = 32;
                }
            }
        }
    } else if (s_lmp_cmp(c2, (char *)"TR", (ftnlen)2, (ftnlen)2) == 0) {
        if (s_lmp_cmp(c3, (char *)"TRI", (ftnlen)3, (ftnlen)3) == 0) {
            if (sname) {
                nb = 64;
            } else {
                nb = 64;
            }
        } else if (s_lmp_cmp(c3, (char *)"EVC", (ftnlen)3, (ftnlen)3) == 0) {
            if (sname) {
                nb = 64;
            } else {
                nb = 64;
            }
        } else if (s_lmp_cmp(c3, (char *)"SYL", (ftnlen)3, (ftnlen)3) == 0) {
            if (sname) {
                i__2 = 48, i__3 = (min(*n1, *n2) << 4) / 100;
                i__1 = max(i__2, i__3);
                nb = min(i__1, 240);
            } else {
                i__2 = 24, i__3 = (min(*n1, *n2) << 3) / 100;
                i__1 = max(i__2, i__3);
                nb = min(i__1, 80);
            }
        }
    } else if (s_lmp_cmp(c2, (char *)"LA", (ftnlen)2, (ftnlen)2) == 0) {
        if (s_lmp_cmp(c3, (char *)"UUM", (ftnlen)3, (ftnlen)3) == 0) {
            if (sname) {
                nb = 64;
            } else {
                nb = 64;
            }
        } else if (s_lmp_cmp(c3, (char *)"TRS", (ftnlen)3, (ftnlen)3) == 0) {
            if (sname) {
                nb = 32;
            } else {
                nb = 32;
            }
        }
    } else if (sname && s_lmp_cmp(c2, (char *)"ST", (ftnlen)2, (ftnlen)2) == 0) {
        if (s_lmp_cmp(c3, (char *)"EBZ", (ftnlen)3, (ftnlen)3) == 0) {
            nb = 1;
        }
    } else if (s_lmp_cmp(c2, (char *)"GG", (ftnlen)2, (ftnlen)2) == 0) {
        nb = 32;
        if (s_lmp_cmp(c3, (char *)"HD3", (ftnlen)3, (ftnlen)3) == 0) {
            if (sname) {
                nb = 32;
            } else {
                nb = 32;
            }
        }
    }
    ret_val = nb;
    return ret_val;
L60:
    nbmin = 2;
    if (s_lmp_cmp(c2, (char *)"GE", (ftnlen)2, (ftnlen)2) == 0) {
        if (s_lmp_cmp(c3, (char *)"QRF", (ftnlen)3, (ftnlen)3) == 0 ||
            s_lmp_cmp(c3, (char *)"RQF", (ftnlen)3, (ftnlen)3) == 0 ||
            s_lmp_cmp(c3, (char *)"LQF", (ftnlen)3, (ftnlen)3) == 0 ||
            s_lmp_cmp(c3, (char *)"QLF", (ftnlen)3, (ftnlen)3) == 0) {
            if (sname) {
                nbmin = 2;
            } else {
                nbmin = 2;
            }
        } else if (s_lmp_cmp(c3, (char *)"HRD", (ftnlen)3, (ftnlen)3) == 0) {
            if (sname) {
                nbmin = 2;
            } else {
                nbmin = 2;
            }
        } else if (s_lmp_cmp(c3, (char *)"BRD", (ftnlen)3, (ftnlen)3) == 0) {
            if (sname) {
                nbmin = 2;
            } else {
                nbmin = 2;
            }
        } else if (s_lmp_cmp(c3, (char *)"TRI", (ftnlen)3, (ftnlen)3) == 0) {
            if (sname) {
                nbmin = 2;
            } else {
                nbmin = 2;
            }
        }
    } else if (s_lmp_cmp(c2, (char *)"SY", (ftnlen)2, (ftnlen)2) == 0) {
        if (s_lmp_cmp(c3, (char *)"TRF", (ftnlen)3, (ftnlen)3) == 0) {
            if (sname) {
                nbmin = 8;
            } else {
                nbmin = 8;
            }
        } else if (sname && s_lmp_cmp(c3, (char *)"TRD", (ftnlen)3, (ftnlen)3) == 0) {
            nbmin = 2;
        }
    } else if (cname && s_lmp_cmp(c2, (char *)"HE", (ftnlen)2, (ftnlen)2) == 0) {
        if (s_lmp_cmp(c3, (char *)"TRD", (ftnlen)3, (ftnlen)3) == 0) {
            nbmin = 2;
        }
    } else if (sname && s_lmp_cmp(c2, (char *)"OR", (ftnlen)2, (ftnlen)2) == 0) {
        if (*(unsigned char *)c3 == 'G') {
            if (s_lmp_cmp(c4, (char *)"QR", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"RQ", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"LQ", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"QL", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"HR", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"TR", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"BR", (ftnlen)2, (ftnlen)2) == 0) {
                nbmin = 2;
            }
        } else if (*(unsigned char *)c3 == 'M') {
            if (s_lmp_cmp(c4, (char *)"QR", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"RQ", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"LQ", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"QL", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"HR", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"TR", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"BR", (ftnlen)2, (ftnlen)2) == 0) {
                nbmin = 2;
            }
        }
    } else if (cname && s_lmp_cmp(c2, (char *)"UN", (ftnlen)2, (ftnlen)2) == 0) {
        if (*(unsigned char *)c3 == 'G') {
            if (s_lmp_cmp(c4, (char *)"QR", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"RQ", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"LQ", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"QL", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"HR", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"TR", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"BR", (ftnlen)2, (ftnlen)2) == 0) {
                nbmin = 2;
            }
        } else if (*(unsigned char *)c3 == 'M') {
            if (s_lmp_cmp(c4, (char *)"QR", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"RQ", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"LQ", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"QL", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"HR", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"TR", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"BR", (ftnlen)2, (ftnlen)2) == 0) {
                nbmin = 2;
            }
        }
    } else if (s_lmp_cmp(c2, (char *)"GG", (ftnlen)2, (ftnlen)2) == 0) {
        nbmin = 2;
        if (s_lmp_cmp(c3, (char *)"HD3", (ftnlen)3, (ftnlen)3) == 0) {
            nbmin = 2;
        }
    }
    ret_val = nbmin;
    return ret_val;
L70:
    nx = 0;
    if (s_lmp_cmp(c2, (char *)"GE", (ftnlen)2, (ftnlen)2) == 0) {
        if (s_lmp_cmp(c3, (char *)"QRF", (ftnlen)3, (ftnlen)3) == 0 ||
            s_lmp_cmp(c3, (char *)"RQF", (ftnlen)3, (ftnlen)3) == 0 ||
            s_lmp_cmp(c3, (char *)"LQF", (ftnlen)3, (ftnlen)3) == 0 ||
            s_lmp_cmp(c3, (char *)"QLF", (ftnlen)3, (ftnlen)3) == 0) {
            if (sname) {
                nx = 128;
            } else {
                nx = 128;
            }
        } else if (s_lmp_cmp(c3, (char *)"HRD", (ftnlen)3, (ftnlen)3) == 0) {
            if (sname) {
                nx = 128;
            } else {
                nx = 128;
            }
        } else if (s_lmp_cmp(c3, (char *)"BRD", (ftnlen)3, (ftnlen)3) == 0) {
            if (sname) {
                nx = 128;
            } else {
                nx = 128;
            }
        }
    } else if (s_lmp_cmp(c2, (char *)"SY", (ftnlen)2, (ftnlen)2) == 0) {
        if (sname && s_lmp_cmp(c3, (char *)"TRD", (ftnlen)3, (ftnlen)3) == 0) {
            nx = 32;
        }
    } else if (cname && s_lmp_cmp(c2, (char *)"HE", (ftnlen)2, (ftnlen)2) == 0) {
        if (s_lmp_cmp(c3, (char *)"TRD", (ftnlen)3, (ftnlen)3) == 0) {
            nx = 32;
        }
    } else if (sname && s_lmp_cmp(c2, (char *)"OR", (ftnlen)2, (ftnlen)2) == 0) {
        if (*(unsigned char *)c3 == 'G') {
            if (s_lmp_cmp(c4, (char *)"QR", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"RQ", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"LQ", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"QL", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"HR", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"TR", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"BR", (ftnlen)2, (ftnlen)2) == 0) {
                nx = 128;
            }
        }
    } else if (cname && s_lmp_cmp(c2, (char *)"UN", (ftnlen)2, (ftnlen)2) == 0) {
        if (*(unsigned char *)c3 == 'G') {
            if (s_lmp_cmp(c4, (char *)"QR", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"RQ", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"LQ", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"QL", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"HR", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"TR", (ftnlen)2, (ftnlen)2) == 0 ||
                s_lmp_cmp(c4, (char *)"BR", (ftnlen)2, (ftnlen)2) == 0) {
                nx = 128;
            }
        }
    } else if (s_lmp_cmp(c2, (char *)"GG", (ftnlen)2, (ftnlen)2) == 0) {
        nx = 128;
        if (s_lmp_cmp(c3, (char *)"HD3", (ftnlen)3, (ftnlen)3) == 0) {
            nx = 128;
        }
    }
    ret_val = nx;
    return ret_val;
L80:
    ret_val = 6;
    return ret_val;
L90:
    ret_val = 2;
    return ret_val;
L100:
    ret_val = (integer)((real)min(*n1, *n2) * (float)1.6);
    return ret_val;
L110:
    ret_val = 1;
    return ret_val;
L120:
    ret_val = 50;
    return ret_val;
L130:
    ret_val = 25;
    return ret_val;
L140:
    ret_val = 1;
    if (ret_val == 1) {
        ret_val = ieeeck_(&c__1, &c_b176, &c_b177);
    }
    return ret_val;
L150:
    ret_val = 1;
    if (ret_val == 1) {
        ret_val = ieeeck_(&c__0, &c_b176, &c_b177);
    }
    return ret_val;
L160:
    ret_val = iparmq_(ispec, name__, opts, n1, n2, n3, n4, name_len, opts_len);
    return ret_val;
}
#ifdef __cplusplus
}
#endif
