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

#include "fix_graphics_labels.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "graphics.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "respa.h"
#include "text_file_reader.h"
#include "tokenizer.h"
#include "update.h"
#include "variable.h"

#include <algorithm>
#include <array>
#include <cstring>

#ifdef LAMMPS_JPEG
#include <jpeglib.h>
#endif

#ifdef LAMMPS_PNG
#include <csetjmp>
#include <png.h>
#include <zlib.h>
#endif

#include "scalable_font.h"

using namespace LAMMPS_NS;
using namespace FixConst;

namespace {

void get_color(const std::string &color, unsigned char *rgb)
{
  if (color == "white") {
    rgb[0] = 255;
    rgb[1] = 255;
    rgb[2] = 255;
  } else if (color == "black") {
    rgb[0] = 0;
    rgb[1] = 0;
    rgb[2] = 0;
  } else if (color == "silver") {
    rgb[0] = 192;
    rgb[1] = 192;
    rgb[2] = 192;
  } else if (color == "darkgray") {
    rgb[0] = 64;
    rgb[1] = 64;
    rgb[2] = 64;
  } else {
    auto val = ValueTokenizer(color, "/");
    auto tmp = val.next_int();
    if ((tmp < 0) || (tmp > 255)) throw TokenizerException("Invalid RGB value", color);
    rgb[0] = (unsigned char) tmp;
    tmp = val.next_int();
    if ((tmp < 0) || (tmp > 255)) throw TokenizerException("Invalid RGB value", color);
    rgb[1] = (unsigned char) tmp;
    tmp = val.next_int();
    if ((tmp < 0) || (tmp > 255)) throw TokenizerException("Invalid RGB value", color);
    rgb[2] = (unsigned char) tmp;
    if (val.has_next()) throw TokenizerException("Extra token", val.next_string());
  }
}

// read image into buffer that is locally allocated with new
// return null pointer if incompatible format or not supported

unsigned char *read_image(FILE *fp, int &width, int &height, const std::string &filename,
                          std::string &info)
{
  if (!fp) return nullptr;
  unsigned char *pixmap = nullptr;

  if (utils::strmatch(filename, "\\.jpg$") || utils::strmatch(filename, "\\.JPG$") ||
      utils::strmatch(filename, "\\.jpeg$") || utils::strmatch(filename, "\\.JPEG$")) {

#if defined(LAMMPS_JPEG)
    struct jpeg_decompress_struct cinfo;
    struct jpeg_error_mgr jerr;

    // initialize for reading from input stream
    jpeg_create_decompress(&cinfo);
    cinfo.err = jpeg_std_error(&jerr);
    cinfo.out_color_space = JCS_RGB;
    jpeg_stdio_src(&cinfo, fp);

    int rv = jpeg_read_header(&cinfo, TRUE);
    if (rv != JPEG_HEADER_OK) return nullptr;

    jpeg_start_decompress(&cinfo);

    // we currently only can handle 8-bit 3-component images
    if ((cinfo.data_precision != 8) || (cinfo.output_components != 3)) return nullptr;

    // read file line-by-line and convert to RGB buffer
    width = cinfo.output_width;
    height = cinfo.output_height;
    pixmap = new unsigned char[3 * width * height];
    auto **scanline = new unsigned char *[height];
    for (int i = 0; i < height; ++i) {
      scanline[0] = &pixmap[(height - 1 - i) * 3 * width];
      jpeg_read_scanlines(&cinfo, scanline, 1);
    }
    delete[] scanline;
    jpeg_destroy_decompress(&cinfo);

    info = fmt::format("{}x{} JPEG file, 8-bit RGB", width, height);
    return pixmap;
#else
    info = "JPEG image format not supported in this LAMMPS binary";
    return nullptr;
#endif

  } else if (utils::strmatch(filename, "\\.png$") || utils::strmatch(filename, "\\.PNG$")) {

#if defined(LAMMPS_PNG)
    png_structp png_ptr = nullptr;
    png_infop info_ptr = nullptr;
    unsigned char sig[8];

    // read and check PNG file signature
    if ((fread(sig, sizeof(unsigned char), 8, fp) != 8) || !png_check_sig(sig, 8)) return nullptr;

    // set up reading from file
    png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    if (!png_ptr) return nullptr;

    info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) {
      png_destroy_read_struct(&png_ptr, nullptr, nullptr);
      return nullptr;
    }

    // set up error handling
    if (setjmp(png_jmpbuf(png_ptr))) {    // NOLINT
      png_destroy_read_struct(&png_ptr, &info_ptr, nullptr);
      delete[] pixmap;
      return nullptr;
    }

    png_uint_32 pngwidth, pngheight;
    int bit_depth, color_type, interlace_type;

    png_init_io(png_ptr, fp);
    png_set_sig_bytes(png_ptr, 8);    /* we already read the 8 signature bytes */
    png_read_info(png_ptr, info_ptr); /* read all PNG info up to image data */
    png_get_IHDR(png_ptr, info_ptr, &pngwidth, &pngheight, &bit_depth, &color_type, &interlace_type,
                 nullptr, nullptr);
    width = pngwidth;
    height = pngheight;

    if ((bit_depth != 8) || (color_type != PNG_COLOR_TYPE_RGB))
      info = fmt::format("{}x{} PNG file, Converted to 8-bit RGB", width, height);
    else
      info = fmt::format("{}x{} PNG file, 8-bit RGB", width, height);

    // convert data to compatible RGB data while reading
    if (color_type == PNG_COLOR_TYPE_PALETTE) png_set_expand(png_ptr);
    if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8) png_set_expand(png_ptr);
    if (png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS)) png_set_expand(png_ptr);
    if (bit_depth == 16) png_set_strip_16(png_ptr);
    if ((color_type == PNG_COLOR_TYPE_GRAY) || (color_type == PNG_COLOR_TYPE_GRAY_ALPHA))
      png_set_gray_to_rgb(png_ptr);
    png_set_strip_alpha(png_ptr);
    png_set_interlace_handling(png_ptr);

    png_read_update_info(png_ptr, info_ptr);

    int rowbytes = png_get_rowbytes(png_ptr, info_ptr);
    int channels = (int) png_get_channels(png_ptr, info_ptr);
    // sanity check
    if ((channels != 3) || (rowbytes != 3 * width)) return nullptr;

    pixmap = new unsigned char[height * width * 3];
    auto *row_pointers = new png_bytep[height];
    for (int i = 0; i < height; ++i) row_pointers[i] = pixmap + (height - 1 - i) * rowbytes;

    png_read_image(png_ptr, row_pointers);
    png_read_end(png_ptr, nullptr);

    // cleanup
    delete[] row_pointers;
    png_destroy_read_struct(&png_ptr, &info_ptr, nullptr);

    return pixmap;
#else
    info = "PNG image format not supported in this LAMMPS binary";
    return nullptr;
#endif

  } else {

    // read file in NetPBM binary or ASCII format
    char buffer[128];
    char *ptr = fgets(buffer, 128, fp);
    if (!ptr || (strlen(buffer) < 3) || (buffer[0] != 'P')) return nullptr;

    // detect binary versus ASCII variant
    bool binary = true;
    if (buffer[1] == '3')
      binary = false;
    else if (buffer[1] != '6')
      return nullptr;

    // skip over optional comments
    do {
      (void) fgets(buffer, 128, fp);
    } while (buffer[0] == '#');

    int rv = sscanf(buffer, "%d%d", &width, &height);
    // don't read invalid or oversize images
    if ((rv != 2) || (width < 1) || (height < 1) || ((width * height) > (1 << 30))) return nullptr;

    int tmp = 0;
    ptr = fgets(buffer, 128, fp);
    rv = sscanf(buffer, "%d", &tmp);
    if ((rv != 1) || (tmp != 255)) return nullptr;

    info = fmt::format("{}x{} PPM {} file, 8-bit RGB", width, height, binary ? "binary" : "text");
    pixmap = new unsigned char[3 * width * height];
    if (binary) {
      // read raw data directly into buffer in the expected order of lines
      // this is the inverse of what Image::write_PPM() does
      for (int y = height - 1; y >= 0; --y) {
        rv = fread(&pixmap[y * width * 3], 3, width, fp);
        if (rv != width) {
          delete[] pixmap;
          return nullptr;
        }
      }
    } else {
      // read file line-by-line and store three RGB values at a time
      auto reader = TextFileReader(fp, "NetPBM ASCII pixmap");
      try {
        int y = height - 1;
        int x = 0;
        auto *line = reader.next_line();
        while (line) {
          auto values = ValueTokenizer(line);

          while (values.has_next()) {
            pixmap[y * width * 3 + 3 * x] = values.next_int();
            pixmap[y * width * 3 + 3 * x + 1] = values.next_int();
            pixmap[y * width * 3 + 3 * x + 2] = values.next_int();
            ++x;
            // next line of pixels
            if (x >= width) {
              --y;
              x = 0;
            }
          }
          line = reader.next_line();
        }
      } catch (std::exception &e) {
        delete[] pixmap;
        return nullptr;
      }
    }

    return pixmap;
  }
}
}    // namespace

/* ---------------------------------------------------------------------- */

#define PARSE_VARIABLE(value, name, index)      \
  if (strstr(arg[index], "v_") == arg[index]) { \
    varflag = 1;                                \
    delete[] name;                              \
    name = utils::strdup(arg[index] + 2);       \
  } else                                        \
    value = utils::numeric(FLERR, arg[index], false, lmp)

FixGraphicsLabels::FixGraphicsLabels(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), imgobjs(nullptr), imgparms(nullptr)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, "fix graphics/labels", error);

  // parse mandatory arg

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  if (nevery <= 0) error->all(FLERR, 3, "Illegal fix graphics/labels nevery value");
  global_freq = nevery;
  dynamic_group_allow = 1;

  // set defaults
  numobjs = 0;
  varflag = 0;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "image") == 0) {
      if (iarg + 5 > narg) utils::missing_cmd_args(FLERR, "fix graphics/labels image", error);

      // clang-format off
      PixmapInfo pix{"", 0.0, {0.0, 0.0, 0.0}, 0, 0, nullptr, {-1.0, -1.0, -1.0}, 1.0,
        -1, -1, -1, -1, nullptr, nullptr, nullptr, nullptr};
      // clang-format on

      // read and store image file with pixmap only on MPI rank 0.
      // must always open in binary mode to avoid data corruption on Windows
      if (comm->me == 0) {
        pix.filename = arg[iarg + 1];
        FILE *fp = fopen(pix.filename.c_str(), "rb");
        if (!fp)
          error->one(FLERR, iarg + 1, "Cannot open fix graphics/labels image file {}: {}",
                     pix.filename, utils::getsyserror());
        pix.timestamp = platform::file_write_time(pix.filename);
        std::string info;
        pix.pixmap = read_image(fp, pix.width, pix.height, pix.filename, info);
        fclose(fp);
        if (!pix.pixmap)
          error->one(FLERR, iarg + 1, "Reading fix graphics/labels image file {} failed: {}",
                     pix.filename, info);

        utils::logmesg(lmp, "Read image from {} file: {} format\n", pix.filename, info);
      }
      PARSE_VARIABLE(pix.pos[0], pix.xstr, iarg + 2);
      PARSE_VARIABLE(pix.pos[1], pix.ystr, iarg + 3);
      PARSE_VARIABLE(pix.pos[2], pix.zstr, iarg + 4);
      iarg += 5;

      // check remaining arguments for optional image arguments
      while (iarg < narg) {
        // next argument is next keyword; exit loop
        if ((strcmp(arg[iarg], "image") == 0) || (strcmp(arg[iarg], "text") == 0)) break;

        if (strcmp(arg[iarg], "scale") == 0) {
          if (iarg + 2 > narg)
            utils::missing_cmd_args(FLERR, "fix graphics/labels image scale", error);
          PARSE_VARIABLE(pix.scale, pix.sstr, iarg + 1);
          iarg += 2;
        } else if (strcmp(arg[iarg], "transcolor") == 0) {
          if (iarg + 2 > narg)
            utils::missing_cmd_args(FLERR, "fix graphics/labels image transcolor", error);
          if (strcmp(arg[iarg + 1], "auto") == 0) {
            if (pix.pixmap) {
              pix.transcolor[0] = pix.pixmap[0];
              pix.transcolor[1] = pix.pixmap[1];
              pix.transcolor[2] = pix.pixmap[2];
            }
          } else if (strcmp(arg[iarg + 1], "none") == 0) {
            pix.transcolor[0] = -255.0;
            pix.transcolor[1] = -255.0;
            pix.transcolor[2] = -255.0;
          } else {
            try {
              unsigned char rgb[3];
              get_color(arg[iarg + 1], rgb);
              pix.transcolor[0] = rgb[0];
              pix.transcolor[1] = rgb[1];
              pix.transcolor[2] = rgb[2];
            } catch (TokenizerException &e) {
              error->all(FLERR, iarg + 1, "Error parsing RGB color value {}: {}", arg[iarg + 1],
                         e.what());
            }
          }
          iarg += 2;
        } else {
          error->all(FLERR, iarg, "Unknown fix graphics/labels image keyword: {}", arg[iarg]);
        }
      }
      pixmaps.emplace_back(pix);

    } else if (strcmp(arg[iarg], "text") == 0) {
      if (iarg + 5 > narg) utils::missing_cmd_args(FLERR, "fix graphics/labels text", error);

      // clang-format off
      TextInfo txt{"", {0.0, 0.0, 0.0}, 0, 0, nullptr, {255, 255, 255}, {192, 192, 192},
                   {192, 192, 192}, {192, 192, 192}, false, 48.0, 0.5,
                   -1, -1, -1, -1, nullptr, nullptr, nullptr, nullptr};
      // clang-format on

      txt.text = arg[iarg + 1];
      if (txt.text.find('$') != std::string::npos) varflag = 1;

      PARSE_VARIABLE(txt.pos[0], txt.xstr, iarg + 2);
      PARSE_VARIABLE(txt.pos[1], txt.ystr, iarg + 3);
      PARSE_VARIABLE(txt.pos[2], txt.zstr, iarg + 4);
      iarg += 5;

      // check remaining arguments for optional image arguments
      while (iarg < narg) {
        // next argument is next keyword; exit loop
        if ((strcmp(arg[iarg], "image") == 0) || (strcmp(arg[iarg], "text") == 0)) break;

        if (strcmp(arg[iarg], "size") == 0) {
          if (iarg + 2 > narg)
            utils::missing_cmd_args(FLERR, "fix graphics/labels text size", error);
          PARSE_VARIABLE(txt.size, txt.sstr, iarg + 1);
          // for sizes 4 to 64, text is rendered at 2x2 size and scaled down for anti-aliasing.
          // for larger sizes, the image is rendered at max supported size and scaled as needed.
          txt.size *= 2.0;
          if ((txt.size < 8.0) || (txt.size > 1024.0))
            error->all(FLERR, iarg + 1, "Invalid fix graphics/labels text size value: {}",
                       txt.size * 0.5);
          if (txt.size > 128.0) {
            txt.scale = txt.size / 256.0;
            txt.size = 128.0;
          } else {
            txt.scale = 0.5;
          }
          iarg += 2;
        } else if (strcmp(arg[iarg], "fontcolor") == 0) {
          if (iarg + 2 > narg)
            utils::missing_cmd_args(FLERR, "fix graphics/labels text fontcolor", error);
          try {
            get_color(arg[iarg + 1], txt.fontcolor);
          } catch (TokenizerException &e) {
            error->all(FLERR, iarg + 1, "Error parsing RGB font color value {}: {}", arg[iarg + 1],
                       e.what());
          }
          iarg += 2;
        } else if (strcmp(arg[iarg], "backcolor") == 0) {
          if (iarg + 2 > narg)
            utils::missing_cmd_args(FLERR, "fix graphics/labels text backcolor", error);
          try {
            get_color(arg[iarg + 1], txt.backcolor);
          } catch (TokenizerException &e) {
            error->all(FLERR, iarg + 1, "Error parsing RGB font color value {}: {}", arg[iarg + 1],
                       e.what());
          }
          iarg += 2;
        } else if (strcmp(arg[iarg], "framecolor") == 0) {
          if (iarg + 2 > narg)
            utils::missing_cmd_args(FLERR, "fix graphics/labels text framecolor", error);
          try {
            get_color(arg[iarg + 1], txt.framecolor);
          } catch (TokenizerException &e) {
            error->all(FLERR, iarg + 1, "Error parsing RGB font color value {}: {}", arg[iarg + 1],
                       e.what());
          }
          iarg += 2;
        } else if (strcmp(arg[iarg], "transcolor") == 0) {
          if (iarg + 2 > narg)
            utils::missing_cmd_args(FLERR, "fix graphics/labels text transcolor", error);
          if (strcmp(arg[iarg + 1], "none") == 0) {
            txt.notrans = true;
          } else {
            try {
              get_color(arg[iarg + 1], txt.transcolor);
            } catch (TokenizerException &e) {
              error->all(FLERR, iarg + 1, "Error parsing RGB font color value {}: {}",
                         arg[iarg + 1], e.what());
            }
          }
          iarg += 2;
        } else {
          error->all(FLERR, iarg, "Unknown fix graphics/labels text keyword: {}", arg[iarg]);
        }
      }
      texts.emplace_back(txt);
    } else {
      error->all(FLERR, iarg, "Unknown fix graphics/labels keyword: {}", arg[iarg]);
    }
  }
}

/* ---------------------------------------------------------------------- */

FixGraphicsLabels::~FixGraphicsLabels()
{
  for (auto &pix : pixmaps) {
    delete[] pix.pixmap;
    delete[] pix.xstr;
    delete[] pix.ystr;
    delete[] pix.zstr;
    delete[] pix.sstr;
  }

  for (auto &txt : texts) {
    delete[] txt.pixmap;
    delete[] txt.xstr;
    delete[] txt.ystr;
    delete[] txt.zstr;
    delete[] txt.sstr;
  }

  memory->destroy(imgobjs);
  memory->destroy(imgparms);
}

/* ---------------------------------------------------------------------- */

int FixGraphicsLabels::setmask()
{
  return END_OF_STEP;
}

/* ---------------------------------------------------------------------- */

void FixGraphicsLabels::init()
{
  for (auto &pix : pixmaps) {
    if (pix.xstr) {
      int ivar = input->variable->find(pix.xstr);
      if (ivar < 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "Variable name {} for fix graphics/labels does not exist", pix.xstr);
      if (input->variable->equalstyle(ivar) == 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "fix graphics/labels variable {} is not equal-style variable", pix.xstr);
      pix.xvar = ivar;
    }
    if (pix.ystr) {
      int ivar = input->variable->find(pix.ystr);
      if (ivar < 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "Variable name {} for fix graphics/labels does not exist", pix.ystr);
      if (input->variable->equalstyle(ivar) == 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "fix graphics/labels variable {} is not equal-style variable", pix.ystr);
      pix.yvar = ivar;
    }
    if (pix.zstr) {
      int ivar = input->variable->find(pix.zstr);
      if (ivar < 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "Variable name {} for fix graphics/labels does not exist", pix.zstr);
      if (input->variable->equalstyle(ivar) == 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "fix graphics/labels variable {} is not equal-style variable", pix.zstr);
      pix.zvar = ivar;
    }
    if (pix.sstr) {
      int ivar = input->variable->find(pix.sstr);
      if (ivar < 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "Variable name {} for fix graphics/labels does not exist", pix.sstr);
      if (input->variable->equalstyle(ivar) == 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "fix graphics/labels variable {} is not equal-style variable", pix.sstr);
      pix.svar = ivar;
    }
  }

  for (auto &txt : texts) {
    if (txt.xstr) {
      int ivar = input->variable->find(txt.xstr);
      if (ivar < 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "Variable name {} for fix graphics/labels does not exist", txt.xstr);
      if (input->variable->equalstyle(ivar) == 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "fix graphics/labels variable {} is not equal-style variable", txt.xstr);
      txt.xvar = ivar;
    }
    if (txt.ystr) {
      int ivar = input->variable->find(txt.ystr);
      if (ivar < 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "Variable name {} for fix graphics/labels does not exist", txt.ystr);
      if (input->variable->equalstyle(ivar) == 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "fix graphics/labels variable {} is not equal-style variable", txt.ystr);
      txt.yvar = ivar;
    }
    if (txt.zstr) {
      int ivar = input->variable->find(txt.zstr);
      if (ivar < 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "Variable name {} for fix graphics/labels does not exist", txt.zstr);
      if (input->variable->equalstyle(ivar) == 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "fix graphics/labels variable {} is not equal-style variable", txt.zstr);
      txt.zvar = ivar;
    }
    if (txt.sstr) {
      int ivar = input->variable->find(txt.sstr);
      if (ivar < 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "Variable name {} for fix graphics/labels does not exist", txt.sstr);
      if (input->variable->equalstyle(ivar) == 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "fix graphics/labels variable {} is not equal-style variable", txt.sstr);
      txt.svar = ivar;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixGraphicsLabels::setup(int)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixGraphicsLabels::end_of_step()
{
  numobjs = pixmaps.size() + texts.size();
  if (numobjs == 0) return;

  if (varflag) modify->clearstep_compute();

  memory->destroy(imgobjs);
  memory->destroy(imgparms);
  memory->create(imgobjs, numobjs, "fix_graphics_labels:imgobjs");
  memory->create(imgparms, numobjs, 11, "fix_graphics_labels:imgparms");

  int n = 0;
  for (auto &pix : pixmaps) {
    if (pix.xstr) pix.pos[0] = input->variable->compute_equal(pix.xvar);
    if (pix.ystr) pix.pos[1] = input->variable->compute_equal(pix.yvar);
    if (pix.zstr) pix.pos[2] = input->variable->compute_equal(pix.zvar);
    if (pix.sstr) pix.scale = input->variable->compute_equal(pix.svar);

    // if image file has been changed since, free old data and re-read file

    if (comm->me == 0) {
      auto timestamp = platform::file_write_time(pix.filename);
      if (pix.timestamp != timestamp) {
        pix.timestamp = timestamp;

        FILE *fp = fopen(pix.filename.c_str(), "rb");
        if (!fp)
          error->one(FLERR, Error::NOLASTLINE, "Cannot open fix graphics/labels image file {}: {}",
                     pix.filename, utils::getsyserror());
        std::string info;
        pix.pixmap = read_image(fp, pix.width, pix.height, pix.filename, info);
        fclose(fp);
        if (!pix.pixmap)
          error->one(FLERR, Error::NOLASTLINE,
                     "Reading fix graphics/labels image file {} failed: {}", pix.filename, info);

        utils::logmesg(lmp, "Re-read image from {} file: {} format\n", pix.filename, info);
      }
    }

    imgobjs[n] = Graphics::PIXMAP;
    imgparms[n][0] = 1;
    imgparms[n][1] = pix.pos[0];
    imgparms[n][2] = pix.pos[1];
    imgparms[n][3] = pix.pos[2];
    imgparms[n][4] = pix.width;
    imgparms[n][5] = pix.height;
    imgparms[n][6] = ubuf((int64_t) pix.pixmap).d;
    imgparms[n][7] = pix.transcolor[0] / 255.0;
    imgparms[n][8] = pix.transcolor[1] / 255.0;
    imgparms[n][9] = pix.transcolor[2] / 255.0;
    imgparms[n][10] = pix.scale;
    ++n;
  }

  // initialize font renderer and load in-memory font

  try {
    SSFN::ScalableFont renderfont;

    for (auto &txt : texts) {
      if (txt.xstr) txt.pos[0] = input->variable->compute_equal(txt.xvar);
      if (txt.ystr) txt.pos[1] = input->variable->compute_equal(txt.yvar);
      if (txt.zstr) txt.pos[2] = input->variable->compute_equal(txt.zvar);

      // text is rasterized at twice the size for some anti-aliasing. clamp to avoid crashes.
      if (txt.sstr) {
        txt.size = 2.0 * input->variable->compute_equal(txt.svar);
        if (txt.size > 128.0) {
          txt.scale = txt.size / 256.0;
          txt.size = 128.0;
        } else {
          txt.size = MAX(txt.size, 8.0);
          txt.scale = 0.5;
        }
      }

      renderfont.select_font(SSFN::FAMILY_SANS, SSFN::STYLE_REGULAR, (int) txt.size);

      // need to render the pixmap if NULL, the size is a variable, or we need to substitute the text
      if (txt.sstr || !txt.pixmap || (txt.text.find('$') != std::string::npos)) {
        auto expanded = txt.text;

        // substitute variables in text
        if (expanded.find('$') != std::string::npos) {
          int ncopy = expanded.length() + 1;
          int nwork = ncopy;
          char *copy = (char *) memory->smalloc(ncopy * sizeof(char), "fix/graphics/labels:copy");
          char *work = (char *) memory->smalloc(nwork * sizeof(char), "fix/graphics/labels:work");
          strncpy(copy, expanded.c_str(), ncopy);
          input->substitute(copy, work, ncopy, nwork, 0);
          expanded = copy;
          memory->sfree(copy);
          memory->sfree(work);
        }

        delete[] txt.pixmap;
        txt.pixmap = renderfont.create_pixmap(expanded, txt.width, txt.height, txt.fontcolor,
                                              txt.framecolor, txt.backcolor);
      }
      imgobjs[n] = Graphics::PIXMAP;
      imgparms[n][0] = 1;
      imgparms[n][1] = txt.pos[0];
      imgparms[n][2] = txt.pos[1];
      imgparms[n][3] = txt.pos[2];
      imgparms[n][4] = txt.width;
      imgparms[n][5] = txt.height;
      imgparms[n][6] = ubuf((int64_t) txt.pixmap).d;
      if (txt.notrans) {
        imgparms[n][7] = imgparms[n][8] = imgparms[n][9] = -1.0;
      } else {
        imgparms[n][7] = (double) txt.transcolor[0] / 255.0;
        imgparms[n][8] = (double) txt.transcolor[1] / 255.0;
        imgparms[n][9] = (double) txt.transcolor[2] / 255.0;
      }

      imgparms[n][10] = txt.scale;
      ++n;
    }
  } catch (const SSFN::SSFNException &e) {
    error->all(FLERR, Error::NOLASTLINE, "Error during font rendering: {}", e.what());
  }

  if (varflag) modify->addstep_compute((update->ntimestep / nevery) * nevery + nevery);
}

/* ----------------------------------------------------------------------
   provide graphics information to dump image
------------------------------------------------------------------------- */

int FixGraphicsLabels::image(int *&objs, double **&parms)
{
  if (comm->me == 0) {
    objs = imgobjs;
    parms = imgparms;
    return numobjs;
  }
  return 0;
}
