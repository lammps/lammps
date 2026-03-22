/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_SCALABLE_FONT_H
#define LMP_SCALABLE_FONT_H

#include <exception>
#include <string>

namespace SSFN {

/* rendering modes */
enum { MODE_NONE, MODE_OUTLINE, MODE_BITMAP };

/* font family group */
enum { FAMILY_SERIF, FAMILY_SANS, FAMILY_DECOR, FAMILY_MONOSPACE, FAMILY_HAND };

/* font style flags */
enum { STYLE_REGULAR, STYLE_BOLD, STYLE_ITALIC };

/** text to pixmap renderer class */
class ScalableFont {
 public:
  ScalableFont();
  ~ScalableFont();

  void select_font(int family, int style, int size);
  unsigned char *create_label(const std::string &text, int &width, int &height,
                              const unsigned char *font, const unsigned char *frame,
                              const unsigned char *back, bool horizontal = true);
  unsigned char *create_colorscale(const std::string &text, int &width, int &height,
                                   const unsigned char *font, const unsigned char *frame,
                                   const unsigned char *back, bool horizontal, int minwidth,
                                   class LAMMPS_NS::Image *image, int mapidx, int tics);

 private:
  void *ctx;
};

/** Font renderer exception class */
class SSFNException : public std::exception {
 public:
  SSFNException() = delete;
  /** Thrown during font processing
   *
   * \param   file    source file where exception was thrown
   * \param   line    line in source file where exception was thrown
   * \param   flag    select error message */
  explicit SSFNException(const std::string &file, int line, int flag);

  /** Retrieve message describing the thrown exception
   *
   * This function provides the message that can be retrieved when the corresponding
   * exception is caught.
   *
   * \return  String with error message */
  [[nodiscard]] const char *what() const noexcept override { return message.c_str(); }

 private:
  std::string message;
};
}    // namespace SSFN

#endif
