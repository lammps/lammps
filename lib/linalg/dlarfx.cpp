#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
int dlarfx_(char *side, integer *m, integer *n, doublereal *v, doublereal *tau, doublereal *c__,
            integer *ldc, doublereal *work, ftnlen side_len)
{
    integer c_dim1, c_offset, i__1;
    integer j;
    doublereal t1, t2, t3, t4, t5, t6, t7, t8, t9, v1, v2, v3, v4, v5, v6, v7, v8, v9, t10, v10,
        sum;
    extern int dlarf_(char *, integer *, integer *, doublereal *, integer *, doublereal *,
                      doublereal *, integer *, doublereal *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    --v;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    if (*tau == 0.) {
        return 0;
    }
    if (lsame_(side, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        switch (*m) {
            case 1:
                goto L10;
            case 2:
                goto L30;
            case 3:
                goto L50;
            case 4:
                goto L70;
            case 5:
                goto L90;
            case 6:
                goto L110;
            case 7:
                goto L130;
            case 8:
                goto L150;
            case 9:
                goto L170;
            case 10:
                goto L190;
        }
        dlarf_(side, m, n, &v[1], &c__1, tau, &c__[c_offset], ldc, &work[1], (ftnlen)1);
        goto L410;
    L10:
        t1 = 1. - *tau * v[1] * v[1];
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            c__[j * c_dim1 + 1] = t1 * c__[j * c_dim1 + 1];
        }
        goto L410;
    L30:
        v1 = v[1];
        t1 = *tau * v1;
        v2 = v[2];
        t2 = *tau * v2;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2];
            c__[j * c_dim1 + 1] -= sum * t1;
            c__[j * c_dim1 + 2] -= sum * t2;
        }
        goto L410;
    L50:
        v1 = v[1];
        t1 = *tau * v1;
        v2 = v[2];
        t2 = *tau * v2;
        v3 = v[3];
        t3 = *tau * v3;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * c__[j * c_dim1 + 3];
            c__[j * c_dim1 + 1] -= sum * t1;
            c__[j * c_dim1 + 2] -= sum * t2;
            c__[j * c_dim1 + 3] -= sum * t3;
        }
        goto L410;
    L70:
        v1 = v[1];
        t1 = *tau * v1;
        v2 = v[2];
        t2 = *tau * v2;
        v3 = v[3];
        t3 = *tau * v3;
        v4 = v[4];
        t4 = *tau * v4;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * c__[j * c_dim1 + 3] +
                  v4 * c__[j * c_dim1 + 4];
            c__[j * c_dim1 + 1] -= sum * t1;
            c__[j * c_dim1 + 2] -= sum * t2;
            c__[j * c_dim1 + 3] -= sum * t3;
            c__[j * c_dim1 + 4] -= sum * t4;
        }
        goto L410;
    L90:
        v1 = v[1];
        t1 = *tau * v1;
        v2 = v[2];
        t2 = *tau * v2;
        v3 = v[3];
        t3 = *tau * v3;
        v4 = v[4];
        t4 = *tau * v4;
        v5 = v[5];
        t5 = *tau * v5;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * c__[j * c_dim1 + 3] +
                  v4 * c__[j * c_dim1 + 4] + v5 * c__[j * c_dim1 + 5];
            c__[j * c_dim1 + 1] -= sum * t1;
            c__[j * c_dim1 + 2] -= sum * t2;
            c__[j * c_dim1 + 3] -= sum * t3;
            c__[j * c_dim1 + 4] -= sum * t4;
            c__[j * c_dim1 + 5] -= sum * t5;
        }
        goto L410;
    L110:
        v1 = v[1];
        t1 = *tau * v1;
        v2 = v[2];
        t2 = *tau * v2;
        v3 = v[3];
        t3 = *tau * v3;
        v4 = v[4];
        t4 = *tau * v4;
        v5 = v[5];
        t5 = *tau * v5;
        v6 = v[6];
        t6 = *tau * v6;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * c__[j * c_dim1 + 3] +
                  v4 * c__[j * c_dim1 + 4] + v5 * c__[j * c_dim1 + 5] + v6 * c__[j * c_dim1 + 6];
            c__[j * c_dim1 + 1] -= sum * t1;
            c__[j * c_dim1 + 2] -= sum * t2;
            c__[j * c_dim1 + 3] -= sum * t3;
            c__[j * c_dim1 + 4] -= sum * t4;
            c__[j * c_dim1 + 5] -= sum * t5;
            c__[j * c_dim1 + 6] -= sum * t6;
        }
        goto L410;
    L130:
        v1 = v[1];
        t1 = *tau * v1;
        v2 = v[2];
        t2 = *tau * v2;
        v3 = v[3];
        t3 = *tau * v3;
        v4 = v[4];
        t4 = *tau * v4;
        v5 = v[5];
        t5 = *tau * v5;
        v6 = v[6];
        t6 = *tau * v6;
        v7 = v[7];
        t7 = *tau * v7;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * c__[j * c_dim1 + 3] +
                  v4 * c__[j * c_dim1 + 4] + v5 * c__[j * c_dim1 + 5] + v6 * c__[j * c_dim1 + 6] +
                  v7 * c__[j * c_dim1 + 7];
            c__[j * c_dim1 + 1] -= sum * t1;
            c__[j * c_dim1 + 2] -= sum * t2;
            c__[j * c_dim1 + 3] -= sum * t3;
            c__[j * c_dim1 + 4] -= sum * t4;
            c__[j * c_dim1 + 5] -= sum * t5;
            c__[j * c_dim1 + 6] -= sum * t6;
            c__[j * c_dim1 + 7] -= sum * t7;
        }
        goto L410;
    L150:
        v1 = v[1];
        t1 = *tau * v1;
        v2 = v[2];
        t2 = *tau * v2;
        v3 = v[3];
        t3 = *tau * v3;
        v4 = v[4];
        t4 = *tau * v4;
        v5 = v[5];
        t5 = *tau * v5;
        v6 = v[6];
        t6 = *tau * v6;
        v7 = v[7];
        t7 = *tau * v7;
        v8 = v[8];
        t8 = *tau * v8;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * c__[j * c_dim1 + 3] +
                  v4 * c__[j * c_dim1 + 4] + v5 * c__[j * c_dim1 + 5] + v6 * c__[j * c_dim1 + 6] +
                  v7 * c__[j * c_dim1 + 7] + v8 * c__[j * c_dim1 + 8];
            c__[j * c_dim1 + 1] -= sum * t1;
            c__[j * c_dim1 + 2] -= sum * t2;
            c__[j * c_dim1 + 3] -= sum * t3;
            c__[j * c_dim1 + 4] -= sum * t4;
            c__[j * c_dim1 + 5] -= sum * t5;
            c__[j * c_dim1 + 6] -= sum * t6;
            c__[j * c_dim1 + 7] -= sum * t7;
            c__[j * c_dim1 + 8] -= sum * t8;
        }
        goto L410;
    L170:
        v1 = v[1];
        t1 = *tau * v1;
        v2 = v[2];
        t2 = *tau * v2;
        v3 = v[3];
        t3 = *tau * v3;
        v4 = v[4];
        t4 = *tau * v4;
        v5 = v[5];
        t5 = *tau * v5;
        v6 = v[6];
        t6 = *tau * v6;
        v7 = v[7];
        t7 = *tau * v7;
        v8 = v[8];
        t8 = *tau * v8;
        v9 = v[9];
        t9 = *tau * v9;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * c__[j * c_dim1 + 3] +
                  v4 * c__[j * c_dim1 + 4] + v5 * c__[j * c_dim1 + 5] + v6 * c__[j * c_dim1 + 6] +
                  v7 * c__[j * c_dim1 + 7] + v8 * c__[j * c_dim1 + 8] + v9 * c__[j * c_dim1 + 9];
            c__[j * c_dim1 + 1] -= sum * t1;
            c__[j * c_dim1 + 2] -= sum * t2;
            c__[j * c_dim1 + 3] -= sum * t3;
            c__[j * c_dim1 + 4] -= sum * t4;
            c__[j * c_dim1 + 5] -= sum * t5;
            c__[j * c_dim1 + 6] -= sum * t6;
            c__[j * c_dim1 + 7] -= sum * t7;
            c__[j * c_dim1 + 8] -= sum * t8;
            c__[j * c_dim1 + 9] -= sum * t9;
        }
        goto L410;
    L190:
        v1 = v[1];
        t1 = *tau * v1;
        v2 = v[2];
        t2 = *tau * v2;
        v3 = v[3];
        t3 = *tau * v3;
        v4 = v[4];
        t4 = *tau * v4;
        v5 = v[5];
        t5 = *tau * v5;
        v6 = v[6];
        t6 = *tau * v6;
        v7 = v[7];
        t7 = *tau * v7;
        v8 = v[8];
        t8 = *tau * v8;
        v9 = v[9];
        t9 = *tau * v9;
        v10 = v[10];
        t10 = *tau * v10;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * c__[j * c_dim1 + 3] +
                  v4 * c__[j * c_dim1 + 4] + v5 * c__[j * c_dim1 + 5] + v6 * c__[j * c_dim1 + 6] +
                  v7 * c__[j * c_dim1 + 7] + v8 * c__[j * c_dim1 + 8] + v9 * c__[j * c_dim1 + 9] +
                  v10 * c__[j * c_dim1 + 10];
            c__[j * c_dim1 + 1] -= sum * t1;
            c__[j * c_dim1 + 2] -= sum * t2;
            c__[j * c_dim1 + 3] -= sum * t3;
            c__[j * c_dim1 + 4] -= sum * t4;
            c__[j * c_dim1 + 5] -= sum * t5;
            c__[j * c_dim1 + 6] -= sum * t6;
            c__[j * c_dim1 + 7] -= sum * t7;
            c__[j * c_dim1 + 8] -= sum * t8;
            c__[j * c_dim1 + 9] -= sum * t9;
            c__[j * c_dim1 + 10] -= sum * t10;
        }
        goto L410;
    } else {
        switch (*n) {
            case 1:
                goto L210;
            case 2:
                goto L230;
            case 3:
                goto L250;
            case 4:
                goto L270;
            case 5:
                goto L290;
            case 6:
                goto L310;
            case 7:
                goto L330;
            case 8:
                goto L350;
            case 9:
                goto L370;
            case 10:
                goto L390;
        }
        dlarf_(side, m, n, &v[1], &c__1, tau, &c__[c_offset], ldc, &work[1], (ftnlen)1);
        goto L410;
    L210:
        t1 = 1. - *tau * v[1] * v[1];
        i__1 = *m;
        for (j = 1; j <= i__1; ++j) {
            c__[j + c_dim1] = t1 * c__[j + c_dim1];
        }
        goto L410;
    L230:
        v1 = v[1];
        t1 = *tau * v1;
        v2 = v[2];
        t2 = *tau * v2;
        i__1 = *m;
        for (j = 1; j <= i__1; ++j) {
            sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)];
            c__[j + c_dim1] -= sum * t1;
            c__[j + (c_dim1 << 1)] -= sum * t2;
        }
        goto L410;
    L250:
        v1 = v[1];
        t1 = *tau * v1;
        v2 = v[2];
        t2 = *tau * v2;
        v3 = v[3];
        t3 = *tau * v3;
        i__1 = *m;
        for (j = 1; j <= i__1; ++j) {
            sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * c__[j + c_dim1 * 3];
            c__[j + c_dim1] -= sum * t1;
            c__[j + (c_dim1 << 1)] -= sum * t2;
            c__[j + c_dim1 * 3] -= sum * t3;
        }
        goto L410;
    L270:
        v1 = v[1];
        t1 = *tau * v1;
        v2 = v[2];
        t2 = *tau * v2;
        v3 = v[3];
        t3 = *tau * v3;
        v4 = v[4];
        t4 = *tau * v4;
        i__1 = *m;
        for (j = 1; j <= i__1; ++j) {
            sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * c__[j + c_dim1 * 3] +
                  v4 * c__[j + (c_dim1 << 2)];
            c__[j + c_dim1] -= sum * t1;
            c__[j + (c_dim1 << 1)] -= sum * t2;
            c__[j + c_dim1 * 3] -= sum * t3;
            c__[j + (c_dim1 << 2)] -= sum * t4;
        }
        goto L410;
    L290:
        v1 = v[1];
        t1 = *tau * v1;
        v2 = v[2];
        t2 = *tau * v2;
        v3 = v[3];
        t3 = *tau * v3;
        v4 = v[4];
        t4 = *tau * v4;
        v5 = v[5];
        t5 = *tau * v5;
        i__1 = *m;
        for (j = 1; j <= i__1; ++j) {
            sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * c__[j + c_dim1 * 3] +
                  v4 * c__[j + (c_dim1 << 2)] + v5 * c__[j + c_dim1 * 5];
            c__[j + c_dim1] -= sum * t1;
            c__[j + (c_dim1 << 1)] -= sum * t2;
            c__[j + c_dim1 * 3] -= sum * t3;
            c__[j + (c_dim1 << 2)] -= sum * t4;
            c__[j + c_dim1 * 5] -= sum * t5;
        }
        goto L410;
    L310:
        v1 = v[1];
        t1 = *tau * v1;
        v2 = v[2];
        t2 = *tau * v2;
        v3 = v[3];
        t3 = *tau * v3;
        v4 = v[4];
        t4 = *tau * v4;
        v5 = v[5];
        t5 = *tau * v5;
        v6 = v[6];
        t6 = *tau * v6;
        i__1 = *m;
        for (j = 1; j <= i__1; ++j) {
            sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * c__[j + c_dim1 * 3] +
                  v4 * c__[j + (c_dim1 << 2)] + v5 * c__[j + c_dim1 * 5] + v6 * c__[j + c_dim1 * 6];
            c__[j + c_dim1] -= sum * t1;
            c__[j + (c_dim1 << 1)] -= sum * t2;
            c__[j + c_dim1 * 3] -= sum * t3;
            c__[j + (c_dim1 << 2)] -= sum * t4;
            c__[j + c_dim1 * 5] -= sum * t5;
            c__[j + c_dim1 * 6] -= sum * t6;
        }
        goto L410;
    L330:
        v1 = v[1];
        t1 = *tau * v1;
        v2 = v[2];
        t2 = *tau * v2;
        v3 = v[3];
        t3 = *tau * v3;
        v4 = v[4];
        t4 = *tau * v4;
        v5 = v[5];
        t5 = *tau * v5;
        v6 = v[6];
        t6 = *tau * v6;
        v7 = v[7];
        t7 = *tau * v7;
        i__1 = *m;
        for (j = 1; j <= i__1; ++j) {
            sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * c__[j + c_dim1 * 3] +
                  v4 * c__[j + (c_dim1 << 2)] + v5 * c__[j + c_dim1 * 5] +
                  v6 * c__[j + c_dim1 * 6] + v7 * c__[j + c_dim1 * 7];
            c__[j + c_dim1] -= sum * t1;
            c__[j + (c_dim1 << 1)] -= sum * t2;
            c__[j + c_dim1 * 3] -= sum * t3;
            c__[j + (c_dim1 << 2)] -= sum * t4;
            c__[j + c_dim1 * 5] -= sum * t5;
            c__[j + c_dim1 * 6] -= sum * t6;
            c__[j + c_dim1 * 7] -= sum * t7;
        }
        goto L410;
    L350:
        v1 = v[1];
        t1 = *tau * v1;
        v2 = v[2];
        t2 = *tau * v2;
        v3 = v[3];
        t3 = *tau * v3;
        v4 = v[4];
        t4 = *tau * v4;
        v5 = v[5];
        t5 = *tau * v5;
        v6 = v[6];
        t6 = *tau * v6;
        v7 = v[7];
        t7 = *tau * v7;
        v8 = v[8];
        t8 = *tau * v8;
        i__1 = *m;
        for (j = 1; j <= i__1; ++j) {
            sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * c__[j + c_dim1 * 3] +
                  v4 * c__[j + (c_dim1 << 2)] + v5 * c__[j + c_dim1 * 5] +
                  v6 * c__[j + c_dim1 * 6] + v7 * c__[j + c_dim1 * 7] + v8 * c__[j + (c_dim1 << 3)];
            c__[j + c_dim1] -= sum * t1;
            c__[j + (c_dim1 << 1)] -= sum * t2;
            c__[j + c_dim1 * 3] -= sum * t3;
            c__[j + (c_dim1 << 2)] -= sum * t4;
            c__[j + c_dim1 * 5] -= sum * t5;
            c__[j + c_dim1 * 6] -= sum * t6;
            c__[j + c_dim1 * 7] -= sum * t7;
            c__[j + (c_dim1 << 3)] -= sum * t8;
        }
        goto L410;
    L370:
        v1 = v[1];
        t1 = *tau * v1;
        v2 = v[2];
        t2 = *tau * v2;
        v3 = v[3];
        t3 = *tau * v3;
        v4 = v[4];
        t4 = *tau * v4;
        v5 = v[5];
        t5 = *tau * v5;
        v6 = v[6];
        t6 = *tau * v6;
        v7 = v[7];
        t7 = *tau * v7;
        v8 = v[8];
        t8 = *tau * v8;
        v9 = v[9];
        t9 = *tau * v9;
        i__1 = *m;
        for (j = 1; j <= i__1; ++j) {
            sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * c__[j + c_dim1 * 3] +
                  v4 * c__[j + (c_dim1 << 2)] + v5 * c__[j + c_dim1 * 5] +
                  v6 * c__[j + c_dim1 * 6] + v7 * c__[j + c_dim1 * 7] +
                  v8 * c__[j + (c_dim1 << 3)] + v9 * c__[j + c_dim1 * 9];
            c__[j + c_dim1] -= sum * t1;
            c__[j + (c_dim1 << 1)] -= sum * t2;
            c__[j + c_dim1 * 3] -= sum * t3;
            c__[j + (c_dim1 << 2)] -= sum * t4;
            c__[j + c_dim1 * 5] -= sum * t5;
            c__[j + c_dim1 * 6] -= sum * t6;
            c__[j + c_dim1 * 7] -= sum * t7;
            c__[j + (c_dim1 << 3)] -= sum * t8;
            c__[j + c_dim1 * 9] -= sum * t9;
        }
        goto L410;
    L390:
        v1 = v[1];
        t1 = *tau * v1;
        v2 = v[2];
        t2 = *tau * v2;
        v3 = v[3];
        t3 = *tau * v3;
        v4 = v[4];
        t4 = *tau * v4;
        v5 = v[5];
        t5 = *tau * v5;
        v6 = v[6];
        t6 = *tau * v6;
        v7 = v[7];
        t7 = *tau * v7;
        v8 = v[8];
        t8 = *tau * v8;
        v9 = v[9];
        t9 = *tau * v9;
        v10 = v[10];
        t10 = *tau * v10;
        i__1 = *m;
        for (j = 1; j <= i__1; ++j) {
            sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * c__[j + c_dim1 * 3] +
                  v4 * c__[j + (c_dim1 << 2)] + v5 * c__[j + c_dim1 * 5] +
                  v6 * c__[j + c_dim1 * 6] + v7 * c__[j + c_dim1 * 7] +
                  v8 * c__[j + (c_dim1 << 3)] + v9 * c__[j + c_dim1 * 9] +
                  v10 * c__[j + c_dim1 * 10];
            c__[j + c_dim1] -= sum * t1;
            c__[j + (c_dim1 << 1)] -= sum * t2;
            c__[j + c_dim1 * 3] -= sum * t3;
            c__[j + (c_dim1 << 2)] -= sum * t4;
            c__[j + c_dim1 * 5] -= sum * t5;
            c__[j + c_dim1 * 6] -= sum * t6;
            c__[j + c_dim1 * 7] -= sum * t7;
            c__[j + (c_dim1 << 3)] -= sum * t8;
            c__[j + c_dim1 * 9] -= sum * t9;
            c__[j + c_dim1 * 10] -= sum * t10;
        }
        goto L410;
    }
L410:
    return 0;
}
#ifdef __cplusplus
}
#endif
