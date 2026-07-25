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

/* ----------------------------------------------------------------------
   Writer for the subset of the VTK file formats used by LAMMPS.  The
   layouts implemented here were determined from the output of VTK 9.2 so
   that files written without the VTK library are read back the same way.
   The relevant details, none of which are obvious from the format
   documentation:

   legacy format
   * the version line says 5.1.  that revision replaced the old inline cell
     lists with separate OFFSETS and CONNECTIVITY arrays.
   * data arrays that are not the active scalars end up in a single
     "FIELD FieldData <n>" section, one header line plus values per array.
   * binary payloads follow the newline of their keyword line as raw
     BIG endian values and are followed by another newline.
   * strings carry a length prefix.  lengths below 64 use the single byte
     0xc0|len; longer strings use undocumented multi byte forms, so we
     refuse to write those rather than guess.

   XML format
   * the file version is 0.1 and the byte count header type is UInt32.
   * inline binary data is base64 encoded.  without compression the header
     and the data form a single encoded stream, but with compression the
     header is encoded separately and the two encodings are concatenated.
   * the compressed header is [nblocks][block size][size of last partial
     block][compressed size per block].  we always use a single block.
   * cell offsets are END offsets and have one entry per cell, unlike the
     legacy format where they start at zero and have one extra entry.
   * string arrays use <Array>, not <DataArray>, and their ASCII form is a
     list of decimal character codes with a zero terminating each string.
------------------------------------------------------------------------- */

#include "vtk_writer.h"

#include "safe_pointers.h"
#include "utils.h"

#include <cstdint>
#include <cstring>

#ifdef LAMMPS_ZLIB
#include <zlib.h>
#endif

using namespace LAMMPS_NS;

static constexpr int PER_LINE = 9;
static constexpr int VTK_VERTEX = 1;
static constexpr int VTK_HEXAHEDRON = 12;

// longest string we can write into a binary legacy file, see file header

static constexpr std::size_t MAX_BINARY_STRING = 63;

namespace {

const char b64chars[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

bool little_endian()
{
  const int one = 1;
  return *((const char *) &one) == 1;
}

std::string base64(const void *data, std::size_t len)
{
  const auto *d = static_cast<const unsigned char *>(data);
  std::string out;
  out.reserve(((len + 2) / 3) * 4);

  std::size_t i = 0;
  for (; i + 3 <= len; i += 3) {
    const unsigned v = (d[i] << 16) | (d[i + 1] << 8) | d[i + 2];
    out += b64chars[(v >> 18) & 63];
    out += b64chars[(v >> 12) & 63];
    out += b64chars[(v >> 6) & 63];
    out += b64chars[v & 63];
  }
  if (i + 1 == len) {
    const unsigned v = d[i] << 16;
    out += b64chars[(v >> 18) & 63];
    out += b64chars[(v >> 12) & 63];
    out += "==";
  } else if (i + 2 == len) {
    const unsigned v = (d[i] << 16) | (d[i + 1] << 8);
    out += b64chars[(v >> 18) & 63];
    out += b64chars[(v >> 12) & 63];
    out += b64chars[(v >> 6) & 63];
    out += '=';
  }
  return out;
}

// encode a block of raw bytes the way the VTK library does for inline binary
// XML data.  see the notes at the top of this file for the layout.

std::string encode_bytes(const std::string &raw)
{
#ifdef LAMMPS_ZLIB
  uLongf bound = compressBound(static_cast<uLong>(raw.size()));
  std::vector<unsigned char> buf(bound ? bound : 1);
  if (compress(buf.data(), &bound, reinterpret_cast<const Bytef *>(raw.data()),
               static_cast<uLong>(raw.size())) == Z_OK) {
    const std::uint32_t header[4] = {1u, static_cast<std::uint32_t>(raw.size()), 0u,
                                     static_cast<std::uint32_t>(bound)};
    return base64(header, sizeof(header)) + base64(buf.data(), bound);
  }
  // fall through to writing the data uncompressed if zlib fails
#endif
  std::string block;
  block.reserve(raw.size() + sizeof(std::uint32_t));
  const auto nbytes = static_cast<std::uint32_t>(raw.size());
  block.append(reinterpret_cast<const char *>(&nbytes), sizeof(nbytes));
  block.append(raw);
  return base64(block.data(), block.size());
}

// append the raw bytes of a value in the host byte order

template <typename T> void append_raw(std::string &out, T value)
{
  out.append(reinterpret_cast<const char *>(&value), sizeof(T));
}

// write raw values in big endian byte order, as legacy binary files require

template <typename T> void fwrite_be(FILE *fp, const std::vector<T> &values)
{
  if (values.empty()) return;
  std::vector<char> buf(values.size() * sizeof(T));
  if (little_endian()) {
    for (std::size_t i = 0; i < values.size(); ++i) {
      const auto *src = reinterpret_cast<const char *>(&values[i]);
      for (std::size_t b = 0; b < sizeof(T); ++b) buf[i * sizeof(T) + b] = src[sizeof(T) - 1 - b];
    }
  } else {
    std::memcpy(buf.data(), values.data(), buf.size());
  }
  fwrite(buf.data(), 1, buf.size(), fp);
}

// write values of a legacy data section, either as text or as binary

template <typename T, typename F>
void write_legacy_values(FILE *fp, bool binary, const std::vector<T> &values, F &&tostr)
{
  if (binary) {
    fwrite_be(fp, values);
    fputc('\n', fp);
    return;
  }
  for (std::size_t i = 0; i < values.size(); ++i) {
    fputs(tostr(values[i]).c_str(), fp);
    fputc(((i + 1) % PER_LINE == 0) ? '\n' : ' ', fp);
  }
  if (values.size() % PER_LINE != 0) fputc('\n', fp);
}

// build the content of an XML data array, either as text or base64 encoded

template <typename T, typename F>
std::string xml_values(bool binary, const std::vector<T> &values, F &&tostr, int indent)
{
  if (binary) {
    std::string raw;
    raw.reserve(values.size() * sizeof(T));
    for (const auto &value : values) append_raw(raw, value);
    return std::string(indent, ' ') + encode_bytes(raw) + "\n";
  }

  const std::string pad(indent, ' ');
  std::string out;
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i % PER_LINE == 0) out += pad;
    out += tostr(values[i]);
    out += ((i + 1) % PER_LINE == 0) ? '\n' : ' ';
  }
  if (values.size() % PER_LINE != 0) out += '\n';
  return out;
}

std::string fmt_int(int value)
{
  return fmt::format("{}", value);
}

std::string fmt_int64(std::int64_t value)
{
  return fmt::format("{}", value);
}

std::string fmt_uint8(std::uint8_t value)
{
  return fmt::format("{}", static_cast<int>(value));
}

// 11 significant digits is what the VTK library writes for double values

std::string fmt_double(double value)
{
  return fmt::format("{:.11g}", value);
}

}    // namespace

/* ---------------------------------------------------------------------- */

VTKWriter::VTKWriter(Flavor _flavor, bool _binary) :
    flavor(_flavor), binary(_binary), dataset(NONE), title("Generated by LAMMPS"), npoints(0),
    ncells(0), celltype(VTK_VERTEX)
{
  dims[0] = dims[1] = dims[2] = 0;
  origin[0] = origin[1] = origin[2] = 0.0;
  spacing[0] = spacing[1] = spacing[2] = 1.0;
}

/* ---------------------------------------------------------------------- */

void VTKWriter::set_title(const std::string &_title)
{
  // the legacy format reserves a single line of at most 256 characters

  title = _title.substr(0, 256);
  for (auto &c : title)
    if (c == '\n' || c == '\r') c = ' ';
}

/* ---------------------------------------------------------------------- */

void VTKWriter::set_vertex_cells(const std::vector<double> &xyz, Dataset type)
{
  if (dataset != NONE) throw VTKWriterException("VTK writer already has a dataset");
  if (xyz.size() % 3) throw VTKWriterException("VTK point list is not a multiple of 3");

  points = xyz;
  npoints = static_cast<int>(xyz.size() / 3);
  ncells = npoints;
  celltype = VTK_VERTEX;
  dataset = type;
}

void VTKWriter::set_polydata(const std::vector<double> &xyz)
{
  set_vertex_cells(xyz, POLYDATA);
}

void VTKWriter::set_unstructured_grid(const std::vector<double> &xyz)
{
  set_vertex_cells(xyz, UNSTRUCTURED);
}

/* ---------------------------------------------------------------------- */

void VTKWriter::set_hexahedron(const double corners[8][3])
{
  if (dataset != NONE) throw VTKWriterException("VTK writer already has a dataset");

  points.resize(24);
  for (int i = 0; i < 8; ++i)
    for (int j = 0; j < 3; ++j) points[3 * i + j] = corners[i][j];

  npoints = 8;
  ncells = 1;
  celltype = VTK_HEXAHEDRON;
  dataset = UNSTRUCTURED;
}

/* ---------------------------------------------------------------------- */

void VTKWriter::set_rectilinear_grid(const std::vector<double> &xc, const std::vector<double> &yc,
                                     const std::vector<double> &zc)
{
  if (dataset != NONE) throw VTKWriterException("VTK writer already has a dataset");
  if (xc.empty() || yc.empty() || zc.empty())
    throw VTKWriterException("VTK rectilinear grid needs coordinates in all 3 dimensions");

  xcoord = xc;
  ycoord = yc;
  zcoord = zc;
  dims[0] = static_cast<int>(xc.size());
  dims[1] = static_cast<int>(yc.size());
  dims[2] = static_cast<int>(zc.size());
  npoints = dims[0] * dims[1] * dims[2];
  ncells = (dims[0] - 1) * (dims[1] - 1) * (dims[2] - 1);
  dataset = RECTILINEAR;
}

/* ---------------------------------------------------------------------- */

void VTKWriter::set_image_data(const int _dim[3], const double _origin[3], const double _spacing[3])
{
  if (dataset != NONE) throw VTKWriterException("VTK writer already has a dataset");
  if (_dim[0] < 1 || _dim[1] < 1 || _dim[2] < 1)
    throw VTKWriterException("VTK image data needs at least one point per dimension");

  for (int i = 0; i < 3; ++i) {
    dims[i] = _dim[i];
    origin[i] = _origin[i];
    spacing[i] = _spacing[i];
  }
  npoints = dims[0] * dims[1] * dims[2];
  ncells = (dims[0] - 1) * (dims[1] - 1) * (dims[2] - 1);
  if (ncells < 0) ncells = 0;
  dataset = IMAGE;
}

/* ---------------------------------------------------------------------- */

void VTKWriter::add_array(std::vector<DataArray> &arrays, int nitems, const char *kind,
                          DataArray &&array)
{
  if (dataset == NONE) throw VTKWriterException("VTK writer has no dataset selected");
  if (array.ncomp < 1) throw VTKWriterException("VTK data array needs at least one component");

  std::size_t nvalues = 0;
  switch (array.type) {
    case INT:
      nvalues = array.ivalues.size();
      break;
    case DOUBLE:
      nvalues = array.dvalues.size();
      break;
    case STRING:
      nvalues = array.svalues.size();
      break;
  }

  if (nvalues != static_cast<std::size_t>(nitems) * array.ncomp)
    throw VTKWriterException(fmt::format("VTK data array {} has {} values but {} {} times {} "
                                         "components are required",
                                         array.name, nvalues, nitems, kind, array.ncomp));

  arrays.push_back(std::move(array));
}

/* ---------------------------------------------------------------------- */

void VTKWriter::add_point_array(const std::string &name, int ncomp, const std::vector<int> &data)
{
  add_array(point_arrays, npoints, "points", DataArray{name, INT, ncomp, data, {}, {}});
}

void VTKWriter::add_point_array(const std::string &name, int ncomp, const std::vector<double> &data)
{
  add_array(point_arrays, npoints, "points", DataArray{name, DOUBLE, ncomp, {}, data, {}});
}

void VTKWriter::add_point_array(const std::string &name, const std::vector<std::string> &data)
{
  add_array(point_arrays, npoints, "points", DataArray{name, STRING, 1, {}, {}, data});
}

void VTKWriter::add_cell_array(const std::string &name, int ncomp, const std::vector<int> &data)
{
  add_array(cell_arrays, ncells, "cells", DataArray{name, INT, ncomp, data, {}, {}});
}

void VTKWriter::add_cell_array(const std::string &name, int ncomp, const std::vector<double> &data)
{
  add_array(cell_arrays, ncells, "cells", DataArray{name, DOUBLE, ncomp, {}, data, {}});
}

void VTKWriter::add_cell_array(const std::string &name, const std::vector<std::string> &data)
{
  add_array(cell_arrays, ncells, "cells", DataArray{name, STRING, 1, {}, {}, data});
}

/* ---------------------------------------------------------------------- */

void VTKWriter::set_active_scalars(const std::string &name)
{
  for (const auto &array : point_arrays)
    if (array.name == name) {
      scalars = name;
      return;
    }
  for (const auto &array : cell_arrays)
    if (array.name == name) {
      scalars = name;
      return;
    }
  throw VTKWriterException(fmt::format("VTK data array {} was not added to this writer", name));
}

/* ---------------------------------------------------------------------- */

void VTKWriter::write(const std::string &filename)
{
  SafeFilePtr fp(fopen(filename.c_str(), "wb"));
  if (!fp) throw VTKWriterException(fmt::format("Cannot open VTK file {}", filename));
  write(fp);
}

void VTKWriter::write(FILE *fp)
{
  if (dataset == NONE) throw VTKWriterException("VTK writer has no dataset selected");
  if (flavor == LEGACY)
    write_legacy(fp);
  else
    write_xml(fp);
}

/* ----------------------------------------------------------------------
   legacy format
------------------------------------------------------------------------- */

void VTKWriter::write_legacy_cells(FILE *fp, const char *keyword)
{
  // the legacy offsets start at zero, so there is one more of them than cells

  std::vector<std::int64_t> offsets(ncells + 1);
  const int per_cell = (celltype == VTK_HEXAHEDRON) ? 8 : 1;
  for (int i = 0; i <= ncells; ++i) offsets[i] = static_cast<std::int64_t>(i) * per_cell;

  std::vector<std::int64_t> connectivity(static_cast<std::size_t>(ncells) * per_cell);
  for (std::size_t i = 0; i < connectivity.size(); ++i)
    connectivity[i] = static_cast<std::int64_t>(i);

  utils::print(fp, "{} {} {}\n", keyword, ncells + 1, connectivity.size());
  utils::print(fp, "OFFSETS vtktypeint64\n");
  write_legacy_values(fp, binary, offsets, fmt_int64);
  utils::print(fp, "CONNECTIVITY vtktypeint64\n");
  write_legacy_values(fp, binary, connectivity, fmt_int64);
}

/* ---------------------------------------------------------------------- */

void VTKWriter::write_legacy_array_data(FILE *fp, const DataArray &array)
{
  switch (array.type) {
    case INT:
      write_legacy_values(fp, binary, array.ivalues, fmt_int);
      break;
    case DOUBLE:
      write_legacy_values(fp, binary, array.dvalues, fmt_double);
      break;
    case STRING:
      if (binary) {
        std::string raw;
        for (const auto &s : array.svalues) {
          if (s.size() > MAX_BINARY_STRING)
            throw VTKWriterException(
                fmt::format("Cannot write string \"{}\" to a binary legacy VTK file, it is longer "
                            "than {} characters",
                            s, MAX_BINARY_STRING));
          raw += static_cast<char>(0xc0 | s.size());
          raw += s;
        }
        fwrite(raw.data(), 1, raw.size(), fp);
        fputc('\n', fp);
      } else {
        for (const auto &s : array.svalues) utils::print(fp, "{}\n", s);
      }
      break;
  }
}

/* ---------------------------------------------------------------------- */

void VTKWriter::write_legacy_arrays(FILE *fp, const std::vector<DataArray> &arrays,
                                    const char *keyword, int nitems)
{
  if (arrays.empty()) return;
  utils::print(fp, "{} {}\n", keyword, nitems);

  // the active scalars get their own section, everything else is field data

  int nfield = 0;
  for (const auto &array : arrays) {
    if (array.name == scalars) {
      const char *type = (array.type == INT) ? "int" : "double";
      if (array.type == STRING)
        throw VTKWriterException("VTK string arrays cannot be used as active scalars");
      utils::print(fp, "SCALARS {} {}", array.name, type);
      if (array.ncomp > 1) utils::print(fp, " {}", array.ncomp);
      utils::print(fp, "\nLOOKUP_TABLE default\n");
      write_legacy_array_data(fp, array);
    } else {
      ++nfield;
    }
  }
  if (nfield == 0) return;

  utils::print(fp, "FIELD FieldData {}\n", nfield);
  for (const auto &array : arrays) {
    if (array.name == scalars) continue;
    const char *type = "double";
    if (array.type == INT) type = "int";
    if (array.type == STRING) type = "string";
    utils::print(fp, "{} {} {} {}\n", array.name, array.ncomp, nitems, type);
    write_legacy_array_data(fp, array);
  }
}

/* ---------------------------------------------------------------------- */

void VTKWriter::write_legacy(FILE *fp)
{
  utils::print(fp, "# vtk DataFile Version 5.1\n{}\n{}\n", title, binary ? "BINARY" : "ASCII");

  switch (dataset) {
    case POLYDATA:
      utils::print(fp, "DATASET POLYDATA\nPOINTS {} double\n", npoints);
      write_legacy_values(fp, binary, points, fmt_double);
      write_legacy_cells(fp, "VERTICES");
      break;

    case UNSTRUCTURED: {
      utils::print(fp, "DATASET UNSTRUCTURED_GRID\nPOINTS {} double\n", npoints);
      write_legacy_values(fp, binary, points, fmt_double);
      write_legacy_cells(fp, "CELLS");
      std::vector<int> types(ncells, celltype);
      utils::print(fp, "CELL_TYPES {}\n", ncells);
      write_legacy_values(fp, binary, types, fmt_int);
      break;
    }

    case RECTILINEAR:
      utils::print(fp, "DATASET RECTILINEAR_GRID\nDIMENSIONS {} {} {}\n", dims[0], dims[1],
                   dims[2]);
      utils::print(fp, "X_COORDINATES {} double\n", dims[0]);
      write_legacy_values(fp, binary, xcoord, fmt_double);
      utils::print(fp, "Y_COORDINATES {} double\n", dims[1]);
      write_legacy_values(fp, binary, ycoord, fmt_double);
      utils::print(fp, "Z_COORDINATES {} double\n", dims[2]);
      write_legacy_values(fp, binary, zcoord, fmt_double);
      break;

    case IMAGE:
      utils::print(fp, "DATASET STRUCTURED_POINTS\nDIMENSIONS {} {} {}\n", dims[0], dims[1],
                   dims[2]);
      utils::print(fp, "SPACING {} {} {}\n", fmt_double(spacing[0]), fmt_double(spacing[1]),
                   fmt_double(spacing[2]));
      utils::print(fp, "ORIGIN {} {} {}\n", fmt_double(origin[0]), fmt_double(origin[1]),
                   fmt_double(origin[2]));
      break;

    case NONE:
      break;
  }

  write_legacy_arrays(fp, point_arrays, "POINT_DATA", npoints);
  write_legacy_arrays(fp, cell_arrays, "CELL_DATA", ncells);
}

/* ----------------------------------------------------------------------
   XML format
------------------------------------------------------------------------- */

void VTKWriter::write_xml_data_array(FILE *fp, const char *type, const std::string &name, int ncomp,
                                     const std::string &payload, int indent)
{
  const std::string pad(indent, ' ');
  utils::print(fp, R"({}<DataArray type="{}" Name="{}")", pad, type, name);
  if (ncomp > 1) utils::print(fp, R"( NumberOfComponents="{}")", ncomp);
  utils::print(fp,
               R"( format="{}">)"
               "\n",
               binary ? "binary" : "ascii");
  fputs(payload.c_str(), fp);
  utils::print(fp, "{}</DataArray>\n", pad);
}

/* ---------------------------------------------------------------------- */

void VTKWriter::write_xml_array(FILE *fp, const DataArray &array, int indent)
{
  const std::string pad(indent, ' ');

  switch (array.type) {
    case INT:
      write_xml_data_array(fp, "Int32", array.name, array.ncomp,
                           xml_values(binary, array.ivalues, fmt_int, indent + 2), indent);
      break;

    case DOUBLE:
      write_xml_data_array(fp, "Float64", array.name, array.ncomp,
                           xml_values(binary, array.dvalues, fmt_double, indent + 2), indent);
      break;

    case STRING: {
      // string arrays use <Array> and store the characters of all strings,
      // each terminated by a zero byte

      std::string payload;
      if (binary) {
        std::string raw;
        for (const auto &s : array.svalues) {
          raw += s;
          raw += '\0';
        }
        payload = std::string(indent + 2, ' ') + encode_bytes(raw) + "\n";
      } else {
        std::vector<int> codes;
        for (const auto &s : array.svalues) {
          for (const char c : s) codes.push_back(static_cast<unsigned char>(c));
          codes.push_back(0);
        }
        payload = xml_values(false, codes, fmt_int, indent + 2);
      }
      utils::print(fp,
                   R"({}<Array type="String" Name="{}" format="{}">)"
                   "\n",
                   pad, array.name, binary ? "binary" : "ascii");
      fputs(payload.c_str(), fp);
      utils::print(fp, "{}</Array>\n", pad);
      break;
    }
  }
}

/* ---------------------------------------------------------------------- */

void VTKWriter::write_xml_arrays(FILE *fp, const std::vector<DataArray> &arrays, const char *tag,
                                 int indent)
{
  const std::string pad(indent, ' ');
  bool has_scalars = false;
  for (const auto &array : arrays)
    if (array.name == scalars) has_scalars = true;

  if (has_scalars)
    utils::print(fp,
                 R"({}<{} Scalars="{}">)"
                 "\n",
                 pad, tag, scalars);
  else
    utils::print(fp, "{}<{}>\n", pad, tag);

  for (const auto &array : arrays) write_xml_array(fp, array, indent + 2);
  utils::print(fp, "{}</{}>\n", pad, tag);
}

/* ---------------------------------------------------------------------- */

void VTKWriter::write_xml_cells(FILE *fp, const char *tag, int indent)
{
  const std::string pad(indent, ' ');
  const int per_cell = (celltype == VTK_HEXAHEDRON) ? 8 : 1;

  // XML offsets are end offsets, one per cell

  std::vector<std::int64_t> offsets(ncells);
  for (int i = 0; i < ncells; ++i) offsets[i] = static_cast<std::int64_t>(i + 1) * per_cell;

  std::vector<std::int64_t> connectivity(static_cast<std::size_t>(ncells) * per_cell);
  for (std::size_t i = 0; i < connectivity.size(); ++i)
    connectivity[i] = static_cast<std::int64_t>(i);

  utils::print(fp, "{}<{}>\n", pad, tag);
  write_xml_data_array(fp, "Int64", "connectivity", 1,
                       xml_values(binary, connectivity, fmt_int64, indent + 4), indent + 2);
  write_xml_data_array(fp, "Int64", "offsets", 1,
                       xml_values(binary, offsets, fmt_int64, indent + 4), indent + 2);

  if (strcmp(tag, "Cells") == 0) {
    std::vector<std::uint8_t> types(ncells, static_cast<std::uint8_t>(celltype));
    write_xml_data_array(fp, "UInt8", "types", 1, xml_values(binary, types, fmt_uint8, indent + 4),
                         indent + 2);
  }
  utils::print(fp, "{}</{}>\n", pad, tag);
}

/* ---------------------------------------------------------------------- */

void VTKWriter::write_xml(FILE *fp)
{
  const char *gridtype = "PolyData";
  if (dataset == UNSTRUCTURED) gridtype = "UnstructuredGrid";
  if (dataset == RECTILINEAR) gridtype = "RectilinearGrid";
  if (dataset == IMAGE) gridtype = "ImageData";

  const char *byte_order = little_endian() ? "LittleEndian" : "BigEndian";

  utils::print(fp, "<?xml version=\"1.0\"?>\n");
  utils::print(fp, R"(<VTKFile type="{}" version="0.1" byte_order="{}" header_type="UInt32")",
               gridtype, byte_order);
#ifdef LAMMPS_ZLIB
  if (binary) utils::print(fp, R"( compressor="vtkZLibDataCompressor")");
#endif
  utils::print(fp, ">\n");

  // the grid element carries the extent for the two structured dataset types

  if (dataset == RECTILINEAR || dataset == IMAGE) {
    utils::print(fp, R"(  <{} WholeExtent="0 {} 0 {} 0 {}")", gridtype, dims[0] - 1, dims[1] - 1,
                 dims[2] - 1);
    if (dataset == IMAGE)
      utils::print(fp, R"( Origin="{} {} {}" Spacing="{} {} {}")", fmt_double(origin[0]),
                   fmt_double(origin[1]), fmt_double(origin[2]), fmt_double(spacing[0]),
                   fmt_double(spacing[1]), fmt_double(spacing[2]));
    utils::print(fp, ">\n");
    utils::print(fp,
                 R"(    <Piece Extent="0 {} 0 {} 0 {}">)"
                 "\n",
                 dims[0] - 1, dims[1] - 1, dims[2] - 1);
  } else {
    utils::print(fp, "  <{}>\n", gridtype);
    if (dataset == POLYDATA)
      utils::print(fp,
                   R"(    <Piece NumberOfPoints="{}" NumberOfVerts="{}" NumberOfLines="0" )"
                   R"(NumberOfStrips="0" NumberOfPolys="0">)"
                   "\n",
                   npoints, ncells);
    else
      utils::print(fp,
                   R"(    <Piece NumberOfPoints="{}" NumberOfCells="{}">)"
                   "\n",
                   npoints, ncells);
  }

  write_xml_arrays(fp, point_arrays, "PointData", 6);
  write_xml_arrays(fp, cell_arrays, "CellData", 6);

  switch (dataset) {
    case POLYDATA:
      utils::print(fp, "      <Points>\n");
      write_xml_data_array(fp, "Float64", "Points", 3, xml_values(binary, points, fmt_double, 10),
                           8);
      utils::print(fp, "      </Points>\n");
      write_xml_cells(fp, "Verts", 6);
      break;

    case UNSTRUCTURED:
      utils::print(fp, "      <Points>\n");
      write_xml_data_array(fp, "Float64", "Points", 3, xml_values(binary, points, fmt_double, 10),
                           8);
      utils::print(fp, "      </Points>\n");
      write_xml_cells(fp, "Cells", 6);
      break;

    case RECTILINEAR:
      utils::print(fp, "      <Coordinates>\n");
      write_xml_data_array(fp, "Float64", "x", 1, xml_values(binary, xcoord, fmt_double, 10), 8);
      write_xml_data_array(fp, "Float64", "y", 1, xml_values(binary, ycoord, fmt_double, 10), 8);
      write_xml_data_array(fp, "Float64", "z", 1, xml_values(binary, zcoord, fmt_double, 10), 8);
      utils::print(fp, "      </Coordinates>\n");
      break;

    case IMAGE:
    case NONE:
      break;
  }

  utils::print(fp, "    </Piece>\n  </{}>\n</VTKFile>\n", gridtype);
}
