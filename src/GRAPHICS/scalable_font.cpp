/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// adapted from ssfn.h
//
// Copyright (C) 2019 bzt (bztsrc@gitlab)
// https://gitlab.com/bztsrc/scalable-font

// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use, copy,
// modify, merge, publish, distribute, sublicense, and/or sell copies
// of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

// Scalable Screen Font renderer in a single ANSI C/C++ file

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>

#include "fmt/format.h"
#include "image.h"
#include "scalable_font.h"

static constexpr int SSFN_DATA_MAX = 65536;

namespace SSFN {

using ssfn_font_t = struct _ssfn_font_t {
  uint8_t magic[4];  /* SSFN magic bytes */
  uint32_t size;     /* total size in bytes */
  uint8_t family;    /* font family group */
  uint8_t style;     /* font style, zero or OR'd SSFN_STYLE_BOLD and SSFN_STYLE_ITALIC */
  uint8_t quality;   /* quality, defines grid size, 0 - 8 */
  uint8_t features;  /* feature flags, OR'd SSFN_FEAT_* */
  uint8_t revision;  /* format revision, must be zero */
  uint8_t reserved0; /* must be zero */
  uint16_t reserved1;
  uint16_t baseline;  /* horizontal baseline in grid pixels */
  uint16_t underline; /* position of under line in grid pixels */
  uint16_t bbox_left; /* overall bounding box for all glyphs in grid pixels */
  uint16_t bbox_top;
  uint16_t bbox_right;
  uint16_t bbox_bottom;
  uint32_t fragments_offs;  /* offset of fragments table relative to magic */
  uint32_t characters_offs; /* offset of characters tables relative to magic */
  uint32_t kerning_offs;    /* kerning table offset relative to magic */
};

/* returned bitmap struct */
using ssfn_glyph_t = struct _ssfn_glyph_t {
  uint8_t mode;                /* returned glyph's data format */
  uint8_t baseline;            /* baseline of glyph, scaled to size */
  uint8_t w;                   /* width */
  uint8_t h;                   /* height */
  uint8_t adv_x;               /* advance x */
  uint8_t adv_y;               /* advance y */
  uint16_t pitch;              /* data buffer bytes per line */
  uint32_t *cmap;              /* pointer to color map */
  uint8_t data[SSFN_DATA_MAX]; /* data buffer */
};

/* renderer context */

using ssfn_t = struct _ssfn_t {
  const ssfn_font_t **fnt[5];             /* font registry */
  const ssfn_font_t *s;                   /* explicitly selected font */
  const ssfn_font_t *f;                   /* font selected by best match */
  ssfn_glyph_t *ret;                      /* glyph to return */
  uint16_t *p;                            /* outline points */
  uint16_t *r[256];                       /* raster for scanlines */
  uint16_t *h;                            /* auto hinting grid */
  int len[5];                             /* number of fonts in registry */
  int mp;                                 /* memory allocated for points */
  int np;                                 /* how many points actually are there */
  int nr[256];                            /* number of coordinates in each raster line */
  int family;                             /* required family */
  int style;                              /* required style */
  int size;                               /* required size */
  int mode;                               /* required mode */
  int g;                                  /* shift value for grid size */
  int m, ix, u, uix, uax, lx, ly, mx, my; /* helper variables */
};

/***** file format *****/

/* magic bytes */
#define SSFN_MAGIC "SSFN"
#define SSFN_COLLECTION "SFNC"
#define SSFN_ENDMAGIC "NFSS"

/* file format features */
#define SSFN_FEAT_HASBMAP 1 /* there's at least one bitmap fragment */
#define SSFN_FEAT_HASCMAP \
  2 /* there's at least one pixmap fragment or one color command, so a color map too */
#define SSFN_FEAT_HASHINT 4  /* there's at least one hinting fragment */
#define SSFN_FEAT_KBIGLKP 8  /* big offsets in kerning look up table */
#define SSFN_FEAT_KBIGCHR 16 /* big characters in kerning look up tables */
#define SSFN_FEAT_KBIGCRD 32 /* big coordinates in kerning groups */
#define SSFN_FEAT_HBIGCRD 64 /* bit coordinates in autohinting fragments */

/* contour commands for vector fragments */
#define SSFN_CONTOUR_MOVE 0
#define SSFN_CONTOUR_LINE 1
#define SSFN_CONTOUR_QUAD 2
#define SSFN_CONTOUR_CUBIC 3
#define SSFN_CONTOUR_COLOR 4

/* bitmap and pixmap fragments and hinting grid info */
#define SSFN_FRAG_BITMAP 0
#define SSFN_FRAG_LBITMAP 1
#define SSFN_FRAG_PIXMAP 2
#define SSFN_FRAG_HINTING 3

/* main SSFN header */

/***** renderer API *****/
#define SSFN_FAMILY_ANY 0xff    /* select the first loaded font */
#define SSFN_FAMILY_BYNAME 0xfe /* select font by its unique name */

#define SSFN_STYLE_UNDERLINE 4    /* under line glyph */
#define SSFN_STYLE_STHROUGH 8     /* strike through glyph */
#define SSFN_STYLE_NOHINTING 0x40 /* no auto hinting grid */
#define SSFN_STYLE_ABS_SIZE 0x80  /* use absolute size value */

#define SSFN_FRAG_CONTOUR 255

/* grid fitting */
#define SSFN_HINTING_THRESHOLD 16 /* don't change unless you really know what you're doing */

/*** normal renderer (ca. 22k, fully featured with error checking) ***/

namespace {
  /* error codes */
  enum {
    SSFN_OK = 0,       /* success */
    SSFN_ERR_ALLOC,    /* allocation error */
    SSFN_ERR_NOFACE,   /* no font face selected */
    SSFN_ERR_INVINP,   /* invalid input */
    SSFN_ERR_BADFILE,  /* bad SSFN file format */
    SSFN_ERR_BADSTYLE, /* bad style */
    SSFN_ERR_BADSIZE,  /* bad size */
    SSFN_ERR_BADMODE,  /* bad mode */
    SSFN_ERR_NOGLYPH,  /* glyph (or kerning info) not found */
    SSFN_ERR_NOMAP,    /* could not access color map information */
    SSFN_ERR_LAST
  };

  /**
   * Error code strings
   */
  const char *ssfn_errstr[] = {"",
                               "Memory allocation error",
                               "No font face found",
                               "Invalid input value",
                               "Bad file format",
                               "Invalid style",
                               "Invalid size",
                               "Invalid mode",
                               "Glyph not found",
                               "Color map not available"};

// include font data as constant in memory byte sequence
#include "scalable_sans_font.h"
  const ssfn_font_t *const ssfn_sans_font = (ssfn_font_t *) VeraR_sfn;

  /*** Private functions ***/

/* f = file scale, g = grid 4095.15, o = screen point 255.255, i = screen pixel 255, c = ceil */
#define _ssfn_i2g(x) ((x) ? (((x) << 16) - (1 << 15)) / ctx->m : 0)
#define _ssfn_g2o(x) (((x) * ctx->m + (1 << 7)) >> 8)
#define _ssfn_g2i(x) (((x) * ctx->m + (1 << 15)) >> 16)
#define _ssfn_g2ic(x) (((x) * ctx->m + (1 << 16) - 1) >> 16)
#define _ssfn_f2i(x) ((((x) << s) * ctx->m + (1 << 15)) >> 16)
#define _ssfn_o2i(x) (((x) + (1 << 7)) >> 8)
#define _ssfn_o2ic(x) ((x + (1 << 8) - 1) >> 8)
#define _ssfn_g2ox(x) ((x) >= (4095 << 4) ? _ssfn_g2o(x) : ctx->h[((x) >> 4)])
#define _ssfn_g2ix(x) ((x) >= (4095 << 4) ? _ssfn_g2i(x) : _ssfn_o2i(ctx->h[((x) >> 4)]))
#define _ssfn_g2ixc(x) ((x) >= (4095 << 4) ? _ssfn_g2ic(x) : _ssfn_o2ic(ctx->h[((x) >> 4)]))
#define _ssfn_g2oy(y) (_ssfn_g2o(y))
#define _ssfn_g2iy(y) (_ssfn_g2i(y))
#define _ssfn_g2iyc(y) (_ssfn_g2ic(y))
#define _ssfn_igg(y) (((4096 << 4) - (y)) >> (2))
#define _ssfn_igi(y) ((((4096 << 4) - (y)) * ctx->m + (1 << (15 + 3))) >> (16 + 3))

  /* parse character table */
  uint8_t *_ssfn_c(const ssfn_font_t *font, uint32_t unicode)
  {
    uint32_t i, l;
    uint8_t *ptr;

    if (!font->characters_offs) return nullptr;

    ptr = (uint8_t *) font + font->characters_offs;
    l = (font->quality < 5 && font->characters_offs < 65536)
        ? 4
        : (font->characters_offs < 1048576 ? 5 : 6);

    for (i = 0; i < 0x110000; i++) {
      if (ptr[0] & 0x80) {
        if (ptr[0] & 0x40) {
          i += ptr[1] | ((ptr[0] & 0x3f) << 8);
          ptr += 2;
        } else {
          i += ptr[0] & 0x3f;
          ptr++;
        }
      } else {
        if (i == unicode) return ptr;
        ptr += ptr[0] * l + 10;
      }
    }
    return nullptr;
  }

  /* add a line to contour */
  void _ssfn_l(ssfn_t *ctx, int x, int y, int l)
  {
    if (x > (4096 << 4) - 16) x = (4096 << 4) - 16;
    if (y > (4096 << 4) - 16) y = (4096 << 4) - 16;
    if (x < -1 || y < -1 || (x == ctx->lx && y == ctx->ly)) return;

    if (ctx->np + 2 >= ctx->mp) {
      ctx->mp += 512;
      ctx->p = (uint16_t *) realloc(ctx->p, ctx->mp * sizeof(uint16_t));
      if (!ctx->p) throw SSFNException(__FILE__, __LINE__, SSFN_ERR_ALLOC);
    }
    if (!ctx->np || !l || _ssfn_g2i(ctx->p[ctx->np - 2]) != _ssfn_g2i(x) ||
        _ssfn_g2i(ctx->p[ctx->np - 1]) != _ssfn_g2i(y)) {
      ctx->p[ctx->np++] = x;
      ctx->p[ctx->np++] = y;
      ctx->lx = x;
      ctx->ly = y;
    }
    if ((ctx->style & 0x200) && x >= 0 && ctx->ix > x) ctx->ix = x;
  }

  /* add a Bezier curve to contour */
  void _ssfn_b(ssfn_t *ctx, int x0, int y0, int x1, int y1, int x2, int y2, int x3, int y3, int l)
  {
    int m0x, m0y, m1x, m1y, m2x, m2y, m3x, m3y, m4x, m4y, m5x, m5y;
    if (l < 8 && (x0 != x3 || y0 != y3)) {
      m0x = ((x1 - x0) / 2) + x0;
      m0y = ((y1 - y0) / 2) + y0;
      m1x = ((x2 - x1) / 2) + x1;
      m1y = ((y2 - y1) / 2) + y1;
      m2x = ((x3 - x2) / 2) + x2;
      m2y = ((y3 - y2) / 2) + y2;
      m3x = ((m1x - m0x) / 2) + m0x;
      m3y = ((m1y - m0y) / 2) + m0y;
      m4x = ((m2x - m1x) / 2) + m1x;
      m4y = ((m2y - m1y) / 2) + m1y;
      m5x = ((m4x - m3x) / 2) + m3x;
      m5y = ((m4y - m3y) / 2) + m3y;
      _ssfn_b(ctx, x0, y0, m0x, m0y, m3x, m3y, m5x, m5y, l + 1);
      _ssfn_b(ctx, m5x, m5y, m4x, m4y, m2x, m2y, x3, y3, l + 1);
    }
    _ssfn_l(ctx, x3, y3, l);
  }

  /* rasterize contour */
  void _ssfn_r(ssfn_t *ctx)
  {
    int i, k, l, m, n = 0, x, y, Y, M = 0;
    uint16_t *r;
    uint8_t *pix = ctx->ret->data;

    for (y = 0; y < ctx->ret->h; y++) {
      Y = _ssfn_i2g(y);
      r = ctx->r[y];
      for (n = 0, i = 0; i < ctx->np - 3; i += 2) {
        if ((ctx->p[i] == 0xffff && ctx->p[i + 1] == 0xffff) ||
            (ctx->p[i + 2] == 0xffff && ctx->p[i + 3] == 0xffff))
          continue;
        if ((ctx->p[i + 1] < Y && ctx->p[i + 3] >= Y) ||
            (ctx->p[i + 3] < Y && ctx->p[i + 1] >= Y)) {
          if (_ssfn_g2iy(ctx->p[i + 1]) == _ssfn_g2iy(ctx->p[i + 3]))
            x = (((int) ctx->p[i] + (int) ctx->p[i + 2]) >> 1);
          else
            x = ((int) ctx->p[i]) +
                ((Y - (int) ctx->p[i + 1]) * ((int) ctx->p[i + 2] - (int) ctx->p[i]) /
                 ((int) ctx->p[i + 3] - (int) ctx->p[i + 1]));
          if (y == ctx->u) {
            if (x < ctx->uix) { ctx->uix = x; }
            if (x > ctx->uax) { ctx->uax = x; }
          }
          x = _ssfn_g2ox(x - ctx->ix);
          for (k = 0; k < n && x > r[k]; k++);
          if (n >= ctx->nr[y]) {
            ctx->nr[y] = (n < ctx->np) ? ctx->np : (n + 1) << 1;
            ctx->r[y] = (uint16_t *) realloc(ctx->r[y], (ctx->nr[y] << 1));
            if (!ctx->r[y]) throw SSFNException(__FILE__, __LINE__, SSFN_ERR_ALLOC);
            r = ctx->r[y];
          }
          for (l = n; l > k; l--) r[l] = r[l - 1];
          r[k] = x;
          n++;
        }
      }
      if (n > 1 && n & 1) {
        r[n - 2] = r[n - 1];
        n--;
      }
      ctx->nr[y] = n;
      if (n) {
        if (y > M) M = y;
        k = y * ctx->ret->pitch;
        for (i = 0; i < n - 1; i += 2) {
          if (ctx->style & 0x100) {
            l = (r[i] + r[i + 1]) >> 9;
            if (pix[k + (l >> 3)] & (1 << (l & 7)) ? true : false) {
              if (r[i] + 256 > r[i + 1])
                r[i] = r[i + 1] - 128;
              else
                r[i] += 256;
              if (r[i + 1] - 256 > r[i]) r[i + 1] -= 256;
            } else {
              if (i >= n - 2 || r[i + 1] + 256 < r[i + 2]) r[i + 1] += 256;
            }
          }
          l = ((r[i] + 128) >> 8);
          m = ((r[i + 1] + 128) >> 8);
          for (; l < m; l++) pix[k + (l >> 3)] ^= 1 << (l & 7);
          if (l + 1 > ctx->ret->w) ctx->ret->w = l + 1;
        }
      }
    }
    /* fix rounding errors */
    if (M + 1 == ctx->ret->baseline) {
      ctx->ret->baseline--;
      ctx->u--;
    }
  }

  /* parse a glyph */
  void _ssfn_g(ssfn_t *ctx, uint8_t *rg, int render)
  {
    int i, j, nf, m, n, o, ox, oy, ol, t, x, y, a, b, c, d, w, h, s;
    uint8_t *raw, *ra, *re, *pix = ctx->ret->data;

    ctx->lx = ctx->ly = ctx->mx = ctx->my = -1;
    ol = (ctx->f->quality < 5 && ctx->f->characters_offs < 65536)
        ? 4
        : ((ctx->f->characters_offs < 1048576) ? 5 : 6);
    s = 16 - ctx->g;
    nf = rg[0];
    rg += 10;

    for (h = 0; nf--; rg += ol) {
      switch (ol) {
        case 4:
          o = ((rg[1] << 8) | rg[0]);
          ox = rg[2];
          oy = rg[3];
          break;
        case 5:
          o = (((rg[2] & 0xF) << 16) | (rg[1] << 8) | rg[0]);
          ox = (((rg[2] >> 4) & 3) << 8) | rg[3];
          oy = (((rg[2] >> 6) & 3) << 8) | rg[4];
          break;
        default:
          o = (((rg[3] & 0xF) << 24) | (rg[2] << 16) | (rg[1] << 8) | rg[0]);
          ox = ((rg[3] & 0xF) << 8) | rg[4];
          oy = (((rg[3] >> 4) & 0xF) << 8) | rg[5];
          break;
      }
      ox <<= s;
      oy <<= s;
      raw = (uint8_t *) ctx->f + o;
      if (raw[0] & 0x80) {
        t = (raw[0] & 0x60) >> 5;
        if (!render && t != SSFN_FRAG_HINTING) break;
        switch (t) {
          case SSFN_FRAG_LBITMAP:
            x = ((((raw[0] >> 2) & 3) << 8) + raw[1]) + 1;
            y = (((raw[0] & 3) << 8) | raw[2]) + 1;
            raw += 3;
            goto bitmap;

          case SSFN_FRAG_BITMAP:
            x = (raw[0] & 0x1F) + 1;
            y = raw[1] + 1;
            raw += 2;
          bitmap:
            if (ctx->mode == MODE_OUTLINE) {
              x <<= 3;
            outline:
              ctx->lx = ctx->ly = -1;
              x <<= s;
              y <<= s;
              if (ctx->style & 0x200) {
                a = (((4096 << 4) - (oy)) >> (3));
                b = (((4096 << 4) - (oy + y)) >> (3));
                if (ctx->ix > ox + a) ctx->ix = ox + a;
              } else
                a = b = 0;
              if (ctx->np) _ssfn_l(ctx, -1, -1, 0);
              _ssfn_l(ctx, ox + a, oy, 0);
              _ssfn_l(ctx, ox + x + a, oy, 0);
              _ssfn_l(ctx, ox + x + b, oy + y, 0);
              _ssfn_l(ctx, ox + b, oy + y, 0);
              _ssfn_l(ctx, ox + a, oy, 0);
            } else {
              a = x << 3;
              b = y << s;
              c = _ssfn_g2i(oy);
              n = _ssfn_g2i(ox);
              w = _ssfn_g2i(x << (3 + s));
              h = _ssfn_g2i(b);
              if (c + h >= ctx->ret->h) c = ctx->ret->h - h; /* due to rounding */
              c = t = c * ctx->ret->pitch;
              for (j = 0; j < h; j++) {
                o = (((j << 8) * (y << 8) / (h << 8)) >> 8) * x;
                for (i = 0; i < w; i++) {
                  m = ((i << 8) * (a << 8) / (w << 8)) >> 8;
                  if (raw[o + (m >> 3)] & (1 << (m & 7))) {
                    d = n + ((ctx->style & 0x200) ? _ssfn_igi(oy + (j << s) + 127) : 0) + i;
                    pix[c + (d >> 3)] |= 1 << (d & 7);
                    d++;
                    if ((ctx->style & 0x100) && (d >> 3) < ctx->ret->pitch) {
                      pix[c + (d >> 3)] |= 1 << (d & 7);
                      d++;
                      if (ctx->size > 127 && (d >> 3) < ctx->ret->pitch) {
                        pix[c + (d >> 3)] |= 1 << (d & 7);
                      }
                    }
                    if (d > ctx->ret->w) ctx->ret->w = d;
                    d = ox + _ssfn_i2g(d);
                    if ((ctx->style & 0x200) && ctx->ix > d) ctx->ix = d;
                    if (_ssfn_g2i(oy) + j == ctx->u) {
                      if (d < ctx->uix) { ctx->uix = d; }
                      if (d > ctx->uax) { ctx->uax = d; }
                    }
                  }
                }
                c += ctx->ret->pitch;
              }
            }
            h = 0;
            break;

          case SSFN_FRAG_PIXMAP:
            x = (((raw[0] & 12) << 6) | raw[1]) + 1;
            y = (((raw[0] & 3) << 8) | raw[2]) + 1;
            n = ((raw[4] << 8) | raw[3]) + 1;
            raw += 5;
            if (ctx->mode == MODE_OUTLINE) goto outline;
            if (raw[-5] & 0x10) { /* todo: direct ARGB values in pixmap fragment */
            }
            a = x * y;
            if (a >= (ctx->nr[0] << 1)) {
              ctx->nr[0] = (a + 1) >> 1;
              ctx->r[0] = (uint16_t *) realloc(ctx->r[0], (ctx->nr[0] << 1));
              if (!ctx->r[0]) throw SSFNException(__FILE__, __LINE__, SSFN_ERR_ALLOC);
            }
            ctx->ret->cmap = (uint32_t *) ((uint8_t *) ctx->f + ctx->f->size - 964);
            for (re = raw + n, ra = (uint8_t *) ctx->r[0], i = 0; i < a && raw < re;) {
              c = (raw[0] & 0x7F) + 1;
              if (raw[0] & 0x80) {
                for (j = 0; j < c; j++) ra[i++] = raw[1];
                raw += 2;
              } else {
                raw++;
                for (j = 0; j < c; j++) ra[i++] = *raw++;
              }
            }
            b = y << s;
            c = _ssfn_g2i(oy);
            n = _ssfn_g2i(ox);
            w = _ssfn_g2i(x << s);
            h = _ssfn_g2i(b);
            if (c + h >= ctx->ret->h) c = ctx->ret->h - h; /* due to rounding */
            c = t = c * ctx->ret->pitch;
            for (j = 0; j < h; j++) {
              o = (((j << 8) * (y << 8) / (h << 8)) >> 8) * x;
              for (i = 0; i < w; i++) {
                m = ((i << 8) * (x << 8) / (w << 8)) >> 8;
                if (ra[o + m] < 0xF0) {
                  d = n + ((ctx->style & 0x200) ? _ssfn_igi(oy + (j << s) + 127) : 0) + i;
                  re = (uint8_t *) &ctx->ret->cmap[ra[o + m]];
                  a = (re[0] + re[1] + re[2] + 255) >> 2;
                  if (a > 127) pix[c + (d >> 3)] |= 1 << (d & 7);
                  if (d > ctx->ret->w) ctx->ret->w = d;
                }
              }
              c += ctx->ret->pitch;
            }
            h = 0;
            break;

          case SSFN_FRAG_HINTING:
            if (raw[0] & 0x10) {
              n = ((raw[0] & 0xF) << 8) | raw[1];
              raw += 2;
            } else {
              n = raw[0] & 0xF;
              raw++;
            }
            if (render || !ox) {
              raw += n << (ctx->f->features & SSFN_FEAT_HBIGCRD ? 1 : 0);
              continue;
            }
            y = 4096;
            x = ((ox >> s) - 1) << (s - 4);
            ctx->h[y++] = x;
            for (n++; n-- && x < 4096;) {
              x = raw[0];
              raw++;
              if (ctx->f->features & SSFN_FEAT_HBIGCRD) {
                x |= (raw[0] << 8);
                raw++;
              }
              x <<= (s - 4);
              ctx->h[y++] = x;
            }
            if (y < 4096) ctx->h[y++] = 65535;
            h = 1;
            break;
        }
      } else {
        if (!render && h) break;
        if (raw[0] & 0x40) {
          n = ((raw[0] & 0x3F) << 8) | raw[1];
          raw += 2;
        } else {
          n = raw[0] & 0x3F;
          raw++;
        }
        if (ctx->f->quality < 5) {
          x = raw[0];
          y = raw[1];
          raw += 2;
        } else {
          x = ((raw[0] & 3) << 8) | raw[1];
          y = ((raw[0] & 0x30) << 4) | raw[2];
          raw += 3;
        }
        x <<= s;
        y <<= s;
        y += oy;
        x += ox + (ctx->style & 0x200 ? _ssfn_igg(y) : 0);
        if (render) {
          if (ctx->np) {
            _ssfn_l(ctx, ctx->mx, ctx->my, 0);
            _ssfn_l(ctx, -1, -1, 0);
          }
          _ssfn_l(ctx, x, y, 0);
        }
        ctx->lx = ctx->mx = x;
        ctx->ly = ctx->my = y;
        for (n++; n--;) {
          t = ctx->g < 8 ? (raw[0] >> 7) | ((raw[1] >> 6) & 2) : raw[0] & 3;
          x = y = a = b = c = d = j = 0;
          switch (ctx->g) {
            case 4:
            case 5:
            case 6:
            case 7:
              x = raw[0] & 0x7F;
              y = raw[1] & 0x7F;
              switch (t) {
                case 0:
                  raw += raw[0] & 4 ? 5 : 2;
                  break;
                case 1:
                  raw += 2;
                  break;
                case 2:
                  a = raw[2] & 0x7F;
                  b = raw[3] & 0x7F;
                  raw += 4;
                  break;
                case 3:
                  a = raw[2] & 0x7F;
                  b = raw[3] & 0x7F;
                  c = raw[4] & 0x7F;
                  d = raw[5] & 0x7F;
                  raw += 6;
                  break;
              }
              break;

            case 8:
              x = raw[1];
              y = raw[2];
              switch (t) {
                case 0:
                  raw += raw[0] & 4 ? 5 : 2;
                  break;
                case 1:
                  raw += 3;
                  break;
                case 2:
                  a = raw[3];
                  b = raw[4];
                  raw += 5;
                  break;
                case 3:
                  a = raw[3];
                  b = raw[4];
                  c = raw[5];
                  d = raw[6];
                  raw += 7;
                  break;
              }
              break;

            case 9:
              x = ((raw[0] & 4) << 6) | raw[1];
              y = ((raw[0] & 8) << 5) | raw[2];
              switch (t) {
                case 0:
                  raw += raw[0] & 4 ? 5 : 2;
                  break;
                case 1:
                  raw += 3;
                  break;
                case 2:
                  a = ((raw[0] & 16) << 4) | raw[3];
                  b = ((raw[0] & 32) << 3) | raw[4];
                  raw += 5;
                  break;
                case 3:
                  a = ((raw[0] & 16) << 4) | raw[3];
                  b = ((raw[0] & 32) << 3) | raw[4];
                  c = ((raw[0] & 64) << 2) | raw[5];
                  d = ((raw[0] & 128) << 1) | raw[6];
                  raw += 7;
                  break;
              }
              break;

            default:
              x = ((raw[0] & 12) << 6) | raw[1];
              y = ((raw[0] & 48) << 4) | raw[2];
              switch (t) {
                case 0:
                  raw += raw[0] & 4 ? 5 : 2;
                  break;
                case 1:
                  raw += 3;
                  break;
                case 2:
                  a = ((raw[3] & 3) << 8) | raw[4];
                  b = ((raw[3] & 12) << 6) | raw[5];
                  raw += 6;
                  break;
                case 3:
                  a = ((raw[3] & 3) << 8) | raw[4];
                  b = ((raw[3] & 12) << 6) | raw[5];
                  c = ((raw[3] & 48) << 4) | raw[6];
                  d = ((raw[3] & 192) << 2) | raw[7];
                  raw += 8;
                  break;
              }
              break;
          }
          x <<= s;
          y <<= s;
          a <<= s;
          b <<= s;
          c <<= s;
          d <<= s;
          x += ox;
          y += oy;
          a += ox;
          b += oy;
          c += ox;
          d += oy;
          if (ctx->style & 0x200) {
            x += _ssfn_igg(y);
            a += _ssfn_igg(b);
            c += _ssfn_igg(d);
          }
          if (render) {
            switch (t) {
              case 0: /* this v1.0 renderer does not support colored contours */
                break;
              case 1:
                _ssfn_l(ctx, x, y, 0);
                break;
              case 2:
                _ssfn_b(ctx, ctx->lx, ctx->ly, ((a - ctx->lx) >> 1) + ctx->lx,
                        ((b - ctx->ly) >> 1) + ctx->ly, ((x - a) >> 1) + a, ((y - b) >> 1) + b, x,
                        y, 0);
                break;
              case 3:
                _ssfn_b(ctx, ctx->lx, ctx->ly, a, b, c, d, x, y, 0);
                break;
            }
          } else if (t == 1 && x >= 0 && y >= 0) {
            a = ((ctx->lx < x) ? x - ctx->lx : ctx->lx - x) >> 4;
            b = ((ctx->ly < y) ? y - ctx->ly : ctx->ly - y) >> 4;
            c = (ctx->lx + x) >> 5;
            if (a < 2)
              ctx->h[4096 + (!ctx->h[4096 + c] && c && ctx->h[4096 + c - 1] ? c - 1 : c)] += b;
          }
          ctx->lx = x;
          ctx->ly = y;
        }
      }
    }

    if (!render && !h) {
      for (j = m = x = y = 0; j < 4096; j++) {
        if (ctx->h[4096 + j] >= 4096 / SSFN_HINTING_THRESHOLD) {
          if (!j)
            m++;
          else {
            ctx->h[4096 + m++] = j - x;
            x = j;
          }
        }
      }
      if (m < 4096) ctx->h[4096 + m] = 65535;
    }
  }

  /**
 * Load a font or font collection into renderer context
 *
 * @param ctx rendering context
 * @param font SSFN font or font collection in memory
 */
  void _ssfn_load(ssfn_t *ctx, const ssfn_font_t *font)
  {
    ssfn_font_t *ptr, *end;

    if (!ctx || !font) throw SSFNException(__FILE__, __LINE__, SSFN_ERR_INVINP);

    if (!memcmp(font->magic, SSFN_COLLECTION, 4)) {
      end = (ssfn_font_t *) ((uint8_t *) font + font->size);
      for (ptr = (ssfn_font_t *) ((uint8_t *) font + 8); ptr < end;
           ptr = (ssfn_font_t *) ((uint8_t *) ptr + ptr->size))
        _ssfn_load(ctx, ptr);
    } else {
      if (memcmp(font->magic, SSFN_MAGIC, 4) ||
          memcmp((uint8_t *) font + font->size - 4, SSFN_ENDMAGIC, 4) ||
          font->family > FAMILY_HAND || font->fragments_offs > font->size ||
          font->characters_offs > font->size || font->kerning_offs > font->size ||
          font->fragments_offs >= font->characters_offs || font->quality > 8) {
        throw SSFNException(__FILE__, __LINE__, SSFN_ERR_BADFILE);
      } else {
        ctx->len[font->family]++;
        ctx->fnt[font->family] = (const ssfn_font_t **) realloc(
            ctx->fnt[font->family], ctx->len[font->family] * sizeof(void *));
        if (!ctx->fnt[font->family])
          throw SSFNException(__FILE__, __LINE__, SSFN_ERR_ALLOC);
        else
          ctx->fnt[font->family][ctx->len[font->family] - 1] = font;
      }
    }
  }

  /**
 * Set up rendering parameters
 *
 * @param ctx rendering context
 * @param family one of SSFN_FAMILY_*
 * @param name NULL or UTF-8 string if family is SSFN_FAMILY_BYNAME
 * @param style OR'd values of SSFN_STYLE_*
 * @param size how big glyph it should render, 8 - 255
 */
  void _ssfn_select(ssfn_t *ctx, int family, int style, int size)
  {
    if (!ctx) throw SSFNException(__FILE__, __LINE__, SSFN_ERR_INVINP);
    if ((style & ~0xCF)) throw SSFNException(__FILE__, __LINE__, SSFN_ERR_BADSTYLE);
    if (size < 8 || size > 255) throw SSFNException(__FILE__, __LINE__, SSFN_ERR_BADSIZE);

    ctx->np = ctx->mp = 0;
    if (ctx->p) {
      free(ctx->p);
      ctx->p = nullptr;
    }
    ctx->f = nullptr;
    ctx->family = family;
    ctx->style = style;
    ctx->size = size;
    ctx->mode = MODE_BITMAP;
  }

  /**
 * Glyph renderer
 *
 * @param ctx rendering context
 * @param unicode character to render
 * @return newly allocated rasterized glyph
 */
  ssfn_glyph_t *_ssfn_render(ssfn_t *ctx, uint32_t unicode)
  {
    ssfn_font_t **fl;
    int i, j, s, h, p, m, n, bt, bl;
    int l, x, y;
    uint8_t *rg = nullptr, c, d;

    if (!ctx) return nullptr;
    if (ctx->size < 8) throw SSFNException(__FILE__, __LINE__, SSFN_ERR_NOFACE);

    if (ctx->s) {
      ctx->f = (ssfn_font_t *) ctx->s;
      rg = _ssfn_c(ctx->f, unicode);
      if (!rg) rg = _ssfn_c(ctx->f, unicode);
    } else {
      p = ctx->family;
    again:
      if (p == SSFN_FAMILY_ANY) {
        n = 0;
        m = 4;
      } else
        n = m = p;
      for (; n <= m; n++) {
        fl = (ssfn_font_t **) ctx->fnt[n];
        if (ctx->style & 3) {
          /* check if we have a specific ctx->f for the requested style */
          for (i = 0; i < ctx->len[n]; i++)
            if ((fl[i]->style & 3) == (ctx->style & 3) && (rg = _ssfn_c(fl[i], unicode))) {
              ctx->f = fl[i];
              break;
            }
          /* if bold italic was requested, check if we have at least bold or italic */
          if (!rg && (ctx->style & 3) == 3)
            for (i = 0; i < ctx->len[n]; i++)
              if ((fl[i]->style & 3) && (rg = _ssfn_c(fl[i], unicode))) {
                ctx->f = fl[i];
                break;
              }
        }
        /* last resort, get the first ctx->f which has a glyph for this unicode, no matter style */
        if (!rg) {
          for (i = 0; i < ctx->len[n]; i++)
            if ((rg = _ssfn_c(fl[i], unicode))) {
              ctx->f = fl[i];
              break;
            }
        }
      }
      /* if glyph still not found, try any family group */
      if (!rg) {
        if (p != SSFN_FAMILY_ANY) {
          p = SSFN_FAMILY_ANY;
          goto again;
        }
      }
    }
    if (!rg) throw SSFNException(__FILE__, __LINE__, SSFN_ERR_NOGLYPH);

    ctx->style &= 0xFF;
    if ((ctx->style & 1) && !(ctx->f->style & 1)) ctx->style |= 0x100;
    if ((ctx->style & 2) && !(ctx->f->style & 2)) ctx->style |= 0x200;
    if (ctx->f->family == FAMILY_MONOSPACE) ctx->style |= SSFN_STYLE_ABS_SIZE;

    ctx->g = 4 + ctx->f->quality;
    ctx->np = 0;

    s = 16 - ctx->g;
    if (ctx->mode == MODE_OUTLINE) {
      h = ctx->size;
      p = 0;
    } else {
      if (!(((rg[2] & 0x0F) << 8) | rg[6]) || ctx->style & SSFN_STYLE_ABS_SIZE)
        h = ctx->size;
      else
        h = (4096 << 4) * ctx->size / ((ctx->f->baseline - ctx->f->bbox_top) << s);
      p = (h + (ctx->style & 0x100 ? 2 : 0) + (ctx->style & 0x200 ? h >> 2 : 0));
      if (p > 255) {
        p = h = 255;
        if (ctx->style & 0x100) h -= 2;
        if (ctx->style & 0x200) h = h * 4 / 5;
      }
      if (ctx->mode == MODE_BITMAP) p = (p + 7) >> 3;
    }
    ctx->m = h;

    if (!ctx->h) ctx->h = (uint16_t *) malloc(4096 * 2 * sizeof(uint16_t));
    if (!ctx->h) throw SSFNException(__FILE__, __LINE__, SSFN_ERR_ALLOC);

    if (!(ctx->style & SSFN_STYLE_NOHINTING)) {
      memset(&ctx->h[4096], 0, 4096 * sizeof(uint16_t));
      _ssfn_g(ctx, rg, 0);
    } else
      ctx->h[4096] = 65535;

    ctx->h[0] = 0;
    for (i = j = x = y = m = 0; i < 4096 && j < 4095; i++) {
      y = ctx->h[4096 + i] == 65535 ? 4095 - j : ctx->h[4096 + i];
      j += y;
      if (j == x) {
        ctx->h[j] = (((j << 4) * h + (1 << 15)) >> 16) << 8;
        if (!y) j++;
      } else {
        y = _ssfn_g2o(y << 4);
        m += y;
        if ((i & 1) || y < 256) { continue; }
        m &= ~0xFF;
        n = ctx->h[x];
        for (l = 0; l + x <= j; l++) ctx->h[x + l] = n + ((m - n) * l / (j - x));
      }
      x = j;
    }
    ctx->uix = _ssfn_g2ixc((((rg[2] & 0x0F) << 8) | rg[6]) << s);
    ctx->uax = _ssfn_g2iyc((((rg[2] & 0xF0) << 4) | rg[7]) << s);
    if (ctx->mode == MODE_NONE) return nullptr;

    bl = ((((rg[3] & 0x0F) << 8) | rg[8])) << s;
    bt = ((((rg[3] & 0xF0) << 4) | rg[9])) << s;
    i = p * h;
    ctx->ret = (ssfn_glyph_t *) malloc(i + 8 + sizeof(uint8_t *));
    if (!ctx->ret) return nullptr;

    memset(&ctx->ret->data, 0, i);
    ctx->ret->cmap = nullptr;
    ctx->ret->mode = ctx->mode;
    ctx->ret->pitch = p;
    ctx->ret->w = 0;
    ctx->ret->h = h;
    ctx->ret->baseline = (((ctx->f->baseline << s) - bt) * h + (1 << 16) - 1) >> 16;
    ctx->u = ctx->ret->baseline + ((((ctx->f->underline - ctx->f->baseline) << s) * h) >> 16);

    ctx->ret->adv_x = ctx->uix;
    ctx->ret->adv_y = ctx->uax;

    ctx->ix = ctx->uix = 4096 << 4;
    ctx->uax = 0;
    _ssfn_g(ctx, rg, 1);

    if (!(ctx->style & 0x200) || ctx->ix == 4096 << 4) ctx->ix = 0;
    if (ctx->mode == MODE_OUTLINE) {
      if (ctx->np > p * h) {
        ctx->ret = (ssfn_glyph_t *) realloc(ctx->ret, ctx->np + 8 + sizeof(uint8_t *));
        if (!ctx->ret) throw SSFNException(__FILE__, __LINE__, SSFN_ERR_ALLOC);
      }
      for (s = i = 0; i < ctx->np; i += 2) {
        c = ctx->p[i + 0] == 0xffff ? 0xff : _ssfn_g2ix(ctx->p[i + 0] + bl);
        d = ctx->p[i + 1] == 0xffff ? 0xff : _ssfn_g2iy(ctx->p[i + 1] + bt);
        if (s < 2 || ctx->ret->data[s - 2] != c || ctx->ret->data[s - 1] != d) {
          ctx->ret->data[s++] = c;
          ctx->ret->data[s++] = d;
        }
      }
      ctx->ret->pitch = s;
      ctx->ret = (ssfn_glyph_t *) realloc(ctx->ret, s + 8 + sizeof(uint8_t *));
      if (!ctx->ret) throw SSFNException(__FILE__, __LINE__, SSFN_ERR_ALLOC);
    } else {
      _ssfn_r(ctx);
      if (ctx->style & SSFN_STYLE_STHROUGH) {
        if (ctx->ret->w < ctx->ret->adv_x) ctx->ret->w = ctx->ret->adv_x;
        memset(&ctx->ret->data[(ctx->ret->baseline - (ctx->size >> 2)) * p], 0xFF,
               (ctx->size / 64 + 2) * p);
      }
      if (ctx->style & SSFN_STYLE_UNDERLINE) {
        if (ctx->ret->w < ctx->ret->adv_x) ctx->ret->w = ctx->ret->adv_x;
        if (ctx->uax > ctx->ix) ctx->uax = _ssfn_g2i(ctx->uax - ctx->ix);
        if (ctx->uix != 4096 << 4)
          ctx->uix = _ssfn_g2i(ctx->uix - ctx->ix);
        else
          ctx->uix = ctx->ret->w + 3;
        m = ctx->u * p;
        n = ctx->size > 127 ? 2 : 1;
        while (n--) {
          if (ctx->uix > 3) {
            j = ctx->uix - 3;
            if (ctx->mode == MODE_BITMAP)
              for (i = 0; i < j; i++) ctx->ret->data[m + (i >> 3)] |= 1 << (i & 7);
            else
              memset(&ctx->ret->data[m], 0xFF, j);
          }
          if (ctx->uax) {
            j = ctx->uax + 2;
            if (ctx->mode == MODE_BITMAP)
              for (i = j; i < ctx->ret->w; i++) ctx->ret->data[m + (i >> 3)] |= 1 << (i & 7);
            else
              memset(&ctx->ret->data[m + j], 0xFF, p - j);
          }
          m += p;
        }
      }
      if (ctx->ret->adv_y) ctx->ret->baseline = ctx->ret->w >> 1;
    }
    return ctx->ret;
  }

  /**
 * Free renderer context
 *
 * @param ctx rendering context
 */
  void _ssfn_free(ssfn_t *ctx)
  {
    int i;

    if (!ctx) return;

    for (i = 0; i < 5; i++)
      if (ctx->fnt[i]) free(ctx->fnt[i]);
    if (ctx->p) free(ctx->p);
    for (i = 0; i < 256; i++)
      if (ctx->r[i]) free(ctx->r[i]);
    if (ctx->h) free(ctx->h);
    memset(ctx, 0, sizeof(ssfn_t));
  }
  // helper function to truncate a string to a segment starting with "src/";

  std::string truncpath(const std::string &path)
  {
    std::size_t found = path.find("src/");
    if (found != std::string::npos)
      return path.substr(found);
    else
      return path;
  }
}    // namespace

/* public API implementation */

// constructor: create and initialize context and load font
ScalableFont::ScalableFont()
{
  auto *ctxptr = new ssfn_t;
  memset(ctxptr, 0, sizeof(ssfn_t));
  _ssfn_load(ctxptr, ssfn_sans_font);
  ctx = (void *) ctxptr;
}

// clean up
ScalableFont::~ScalableFont()
{
  auto *ctxptr = (ssfn_t *) ctx;
  _ssfn_free(ctxptr);
  delete ctxptr;
}

// set font settings
void ScalableFont::select_font(int family, int style, int size)
{
  _ssfn_select((ssfn_t *) ctx, family, style, size);
}

unsigned char *ScalableFont::create_label(const std::string &text, int &width, int &height,
                                          const unsigned char *font, const unsigned char *frame,
                                          const unsigned char *back, bool horizontal)
{
  auto *ctxptr = (ssfn_t *) ctx;
  ssfn_glyph_t *g;

  // get a font size specific spacing for a border
  g = _ssfn_render(ctxptr, ' ');
  int xspace = g->adv_x;
  free(g);

  // dry run to determine size of pixmap
  width = 0;
  int miny = 1073741824;
  int maxy = 0;
  for (auto c : text + "gll") {    // append these characters for consistent spacing
    if (c == '_') c = ' ';         // ugly hack to work around font issue

    // render character and apply its width
    g = _ssfn_render(ctxptr, c);
    width += g->adv_x;

    // loop over bitmap to find minimum and maximum y position
    for (int y = 0; y < g->h; ++y) {
      const int ypos = g->h - 1 - y + g->baseline;
      for (int x = 0, i = 0, m = 1; x < g->w; ++x, m <<= 1) {
        if (m > 0x80) {
          m = 1;
          ++i;
        }
        if (g->data[y * g->pitch + i] & m) {
          miny = std::min(miny, ypos);
          maxy = std::max(maxy, ypos);
        }
      }
    }
    free(g);
  }

  int xhalf = xspace / 2;
  height = maxy - miny + 1 + 3 * xspace;

  // allocate and fill pixmap with background and frame color

  auto *pixmap = new unsigned char[width * height * 3];
  for (int y = 0; y < height; ++y) {
    int yoffs = 3 * y * width;
    for (int x = 0; x < width; ++x) {
      if ((y < xhalf) || (y >= height - xhalf) || (x < xhalf) || (x >= width - xhalf)) {
        pixmap[yoffs + 3 * x] = frame[0];
        pixmap[yoffs + 3 * x + 1] = frame[1];
        pixmap[yoffs + 3 * x + 2] = frame[2];
      } else {
        pixmap[yoffs + 3 * x] = back[0];
        pixmap[yoffs + 3 * x + 1] = back[1];
        pixmap[yoffs + 3 * x + 2] = back[2];
      }
    }
  }

  // now render each character again and change the pixels in the pixmap accordingly
  int penx = 2 * xspace;
  for (auto c : text) {
    if (c == '_') c = ' ';    // ugly hack to work around font issue

    g = _ssfn_render(ctxptr, c);
    for (int y = 0; y < g->h; ++y) {
      const int yoffs = (g->h - 1 - y + g->baseline - miny + xspace + xhalf / 2) * width * 3;
      for (int x = 0, i = 0, m = 1; x < g->w; ++x, m <<= 1) {
        if (m > 0x80) {
          m = 1;
          ++i;
        }
        const int xoffs = (penx + x) * 3;
        if (g->data[y * g->pitch + i] & m) {
          pixmap[yoffs + xoffs] = font[0];
          pixmap[yoffs + xoffs + 1] = font[1];
          pixmap[yoffs + xoffs + 2] = font[2];
        }
      }
    }
    penx += g->adv_x;
    free(g);
  }

  // for vertical text copy the rotated pixmap accordingly and swap width and height
  if (!horizontal) {
    auto *flipped = new unsigned char[width * height * 3];
    for (int y = 0; y < height; ++y) {
      for (int x = 0; x < width; ++x) {
        for (int i = 0; i < 3; ++i) {
          flipped[x * 3 * height + 3 * (height - 1 - y) + i] = pixmap[y * 3 * width + 3 * x + i];
        }
      }
    }
    delete[] pixmap;
    pixmap = flipped;
    std::swap(width, height);
  }
  return pixmap;
}

unsigned char *ScalableFont::create_colorscale(const std::string &text, int &width, int &height,
                                               const unsigned char *font,
                                               const unsigned char *frame,
                                               const unsigned char *back, bool horizontal,
                                               int minwidth, LAMMPS_NS::Image *image, int mapidx,
                                               int tics)
{
  auto *ctxptr = (ssfn_t *) ctx;
  ssfn_glyph_t *g;

  // get a font size specific spacing for a border
  g = _ssfn_render(ctxptr, ' ');
  int xspace = g->adv_x;
  free(g);
  // get a font size specific spacing for a angular bracket
  g = _ssfn_render(ctxptr, '<');
  int xfill = g->adv_x;
  free(g);

  double lo, hi;
  image->map_info(mapidx, lo, hi);
  auto newtext = fmt::format("{:3.3}  {}  {:3.3}", lo, text, hi);

  // dry run to determine size of pixmap
  width = 0;
  int miny = 1073741824;
  int maxy = 0;
  for (auto c : newtext + "gll") {    // append these characters for consistent spacing
    if (c == '_') c = ' ';            // ugly hack to work around font issue

    // render character and apply its width
    g = _ssfn_render(ctxptr, c);
    width += g->adv_x;

    // loop over bitmap to find minimum and maximum y position
    for (int y = 0; y < g->h; ++y) {
      const int ypos = g->h - 1 - y + g->baseline;
      for (int x = 0, i = 0, m = 1; x < g->w; ++x, m <<= 1) {
        if (m > 0x80) {
          m = 1;
          ++i;
        }
        if (g->data[y * g->pitch + i] & m) {
          miny = std::min(miny, ypos);
          maxy = std::max(maxy, ypos);
        }
      }
    }
    free(g);
  }

  int xhalf = xspace / 2;
  height = maxy - miny + 1 + 6 * xspace;

  if (minwidth > width) {
    int wextra = (minwidth - width) / xfill / 4;
    newtext =
        fmt::format("{:3.3} {:<<{}} {} {:>>{}} {:3.3}", lo, '<', wextra, text, '>', wextra, hi);
    width = minwidth;
  }

  // re-do dry run to determine size of the pixmap with the padded range string
  width = 0;
  miny = 1073741824;
  maxy = 0;
  for (auto c : newtext + "gll") {    // append these characters for consistent spacing
    if (c == '_') c = ' ';            // ugly hack to work around font issue

    // render character and apply its width
    g = _ssfn_render(ctxptr, c);
    width += g->adv_x;

    // loop over bitmap to find minimum and maximum y position
    for (int y = 0; y < g->h; ++y) {
      const int ypos = g->h - 1 - y + g->baseline;
      for (int x = 0, i = 0, m = 1; x < g->w; ++x, m <<= 1) {
        if (m > 0x80) {
          m = 1;
          ++i;
        }
        if (g->data[y * g->pitch + i] & m) {
          miny = std::min(miny, ypos);
          maxy = std::max(maxy, ypos);
        }
      }
    }
    free(g);
  }

  // allocate and fill pixmap with background and frame color

  auto *pixmap = new unsigned char[width * height * 3];
  for (int y = 0; y < height; ++y) {
    int yoffs = 3 * y * width;
    for (int x = 0; x < width; ++x) {
      if ((y < xhalf) || (y >= height - xhalf) || (x < xhalf) || (x >= width - xhalf)) {
        pixmap[yoffs + 3 * x] = frame[0];
        pixmap[yoffs + 3 * x + 1] = frame[1];
        pixmap[yoffs + 3 * x + 2] = frame[2];
      } else {
        pixmap[yoffs + 3 * x] = back[0];
        pixmap[yoffs + 3 * x + 1] = back[1];
        pixmap[yoffs + 3 * x + 2] = back[2];
      }
    }
  }

  // draw colormap
  double delta = (hi - lo) / (width - 2 * xspace);
  int ticinc = tics ? (width - 2 * xspace) / (tics + 1) : 1 << 31;
  int ticmin = ticinc + xspace + 1;
  int ticmax = ticmin + xhalf;
  for (int x = xspace; x < width - xspace; ++x) {
    for (int y = xhalf + xhalf / 2; y < height - 4 * xspace - xhalf; ++y) {
      int offs = 3 * y * width + 3 * x;
      double val = lo + delta * static_cast<double>(x - xspace);
      auto *color = image->map_value2color(mapidx, val);
      if (color) {
        if (tics > 0) {
          if ((((x >= ticmin) && (x <= ticmax)) || (x <= xspace + xhalf) ||
               (x >= width - xspace - xhalf)) ||
              (y <= xspace + xhalf / 2) || (y >= height - 5 * xspace)) {

            pixmap[offs] = font[0];
            pixmap[offs + 1] = font[1];
            pixmap[offs + 2] = font[2];
          } else {
            pixmap[offs] = static_cast<unsigned char>(color[0] * 255.0);
            pixmap[offs + 1] = static_cast<unsigned char>(color[1] * 255.0);
            pixmap[offs + 2] = static_cast<unsigned char>(color[2] * 255.0);
          }
          if (x == ticmax) {
            ticmin += ticinc;
            ticmax += ticinc;
          }
        } else {
          if ((x <= xspace + xhalf) || (x >= width - xspace - xhalf) || (y <= xspace + xhalf / 2) ||
              (y >= height - 5 * xspace)) {
            pixmap[offs] = font[0];
            pixmap[offs + 1] = font[1];
            pixmap[offs + 2] = font[2];
          } else {
            pixmap[offs] = static_cast<unsigned char>(color[0] * 255.0);
            pixmap[offs + 1] = static_cast<unsigned char>(color[1] * 255.0);
            pixmap[offs + 2] = static_cast<unsigned char>(color[2] * 255.0);
          }
        }
      }
    }
  }

  // now render each character again and change the pixels in the pixmap accordingly
  int penx = 2 * xspace;
  int peny = 3 * xspace + xhalf;
  for (auto c : newtext) {
    if (c == '_') c = ' ';    // ugly hack to work around font issue

    g = _ssfn_render(ctxptr, c);
    for (int y = 0; y < g->h; ++y) {
      const int yoffs = (g->h - 1 - y + peny + g->baseline - miny + xspace + xhalf / 2) * width * 3;
      for (int x = 0, i = 0, m = 1; x < g->w; ++x, m <<= 1) {
        if (m > 0x80) {
          m = 1;
          ++i;
        }
        const int xoffs = (penx + x) * 3;
        if (g->data[y * g->pitch + i] & m) {
          pixmap[yoffs + xoffs] = font[0];
          pixmap[yoffs + xoffs + 1] = font[1];
          pixmap[yoffs + xoffs + 2] = font[2];
        }
      }
    }
    penx += g->adv_x;
    free(g);
  }

  // for vertical text copy the rotated pixmap accordingly and swap width and height
  if (!horizontal) {
    auto *flipped = new unsigned char[width * height * 3];
    for (int y = 0; y < height; ++y) {
      for (int x = 0; x < width; ++x) {
        for (int i = 0; i < 3; ++i) {
          flipped[x * 3 * height + 3 * (height - 1 - y) + i] = pixmap[y * 3 * width + 3 * x + i];
        }
      }
    }
    delete[] pixmap;
    pixmap = flipped;
    std::swap(width, height);
  }
  return pixmap;
}

SSFNException::SSFNException(const std::string &file, int line, int flag)
{
  message = fmt::format("In file {}:{} ", truncpath(file), line);
  if ((flag < SSFN_OK) || (flag >= SSFN_ERR_LAST)) {
    message.append("Unknown Error");
  } else {
    message.append(ssfn_errstr[flag]);
  }
}
}    // namespace SSFN
