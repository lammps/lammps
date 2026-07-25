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

// Tests for the VTKWriter class.  These check the structure of the produced
// files and decode the binary payloads, so that the layouts the VTK library
// expects are pinned down without needing that library to be installed.

#include "vtk_writer.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#ifdef LAMMPS_ZLIB
#include <zlib.h>
#endif

using namespace LAMMPS_NS;
using testing::HasSubstr;
using testing::Not;

namespace {

// read a whole file, including any embedded null bytes

std::string slurp(const std::string &filename)
{
    std::ifstream in(filename, std::ios::binary);
    std::ostringstream buf;
    buf << in.rdbuf();
    return buf.str();
}

std::string base64_decode(const std::string &text)
{
    auto value = [](char c) -> int {
        if (c >= 'A' && c <= 'Z') return c - 'A';
        if (c >= 'a' && c <= 'z') return c - 'a' + 26;
        if (c >= '0' && c <= '9') return c - '0' + 52;
        if (c == '+') return 62;
        if (c == '/') return 63;
        return -1;
    };

    std::string out;
    int bits = 0, nbits = 0;
    for (const char c : text) {
        const int v = value(c);
        if (v < 0) continue; // skip padding and whitespace
        bits = (bits << 6) | v;
        nbits += 6;
        if (nbits >= 8) {
            nbits -= 8;
            out += static_cast<char>((bits >> nbits) & 0xff);
        }
    }
    return out;
}

// pull the text between the first ">" after "marker" and the following "<"

std::string payload_after(const std::string &file, const std::string &marker)
{
    auto pos = file.find(marker);
    if (pos == std::string::npos) return "";
    pos = file.find('>', pos);
    if (pos == std::string::npos) return "";
    const auto end = file.find('<', pos);
    return file.substr(pos + 1, end - pos - 1);
}

std::vector<double> parse_numbers(const std::string &text)
{
    std::istringstream in(text);
    std::vector<double> out;
    double v;
    while (in >> v)
        out.push_back(v);
    return out;
}

// decode the content of an inline binary XML data array back to raw bytes.
// without compression the byte count header and the data share one base64
// stream; with compression the header is encoded separately.

std::string decode_binary_payload(const std::string &raw_payload, std::uint32_t expected_bytes)
{
    // drop the indentation and line breaks so that character offsets into the
    // encoded stream are meaningful
    std::string payload;
    for (const char c : raw_payload)
        if (!isspace(static_cast<unsigned char>(c))) payload += c;

    const std::string first = base64_decode(payload);
    EXPECT_GE(first.size(), 4u);
    if (first.size() < 4) return "";

    std::uint32_t header[4] = {0, 0, 0, 0};
    memcpy(header, first.data(), std::min<std::size_t>(16, first.size()));

    // uncompressed: a single stream of [byte count][data]
    if (header[0] == expected_bytes && first.size() == expected_bytes + 4) return first.substr(4);

    // compressed: [nblocks][block size][last partial block][compressed size]
    EXPECT_EQ(header[0], 1u) << "expected exactly one compressed block";
    EXPECT_EQ(header[1], expected_bytes) << "uncompressed size in block header";
    EXPECT_EQ(header[2], 0u) << "no partial last block for a single block";

    // the header occupies the first 24 base64 characters, the data follows
    const std::string data = base64_decode(payload.substr(24));
    EXPECT_EQ(data.size(), header[3]) << "compressed size in block header";

#ifdef LAMMPS_ZLIB
    std::string out(expected_bytes, '\0');
    uLongf len = expected_bytes;
    const int rv =
        uncompress(reinterpret_cast<Bytef *>(&out[0]), &len,
                   reinterpret_cast<const Bytef *>(data.data()), static_cast<uLong>(data.size()));
    EXPECT_EQ(rv, Z_OK);
    EXPECT_EQ(len, expected_bytes);
    return out;
#else
    return "";
#endif
}

double be_double(const char *p)
{
    char buf[8];
    for (int i = 0; i < 8; ++i)
        buf[i] = p[7 - i];
    double v;
    memcpy(&v, buf, 8);
    return v;
}

std::int32_t be_int32(const char *p)
{
    char buf[4];
    for (int i = 0; i < 4; ++i)
        buf[i] = p[3 - i];
    std::int32_t v;
    memcpy(&v, buf, 4);
    return v;
}

// a small reusable data set: 4 points with a scalar, a vector and names

const std::vector<double> COORDS = {0.0, 0.0, 0.0, 1.5, 0.0, 0.0, 1.5, 2.5, 0.0, 0.0, 2.5, -3.25};
const std::vector<int> IDS       = {1, 2, 3, 4};
const std::vector<double> VECS   = {0.5, -1.5, 2.5, 1.0, -2.0, 3.0, 1.5, -2.5, 3.5, 2.0, -3.0, 4.0};
const std::vector<std::string> NAMES = {"C", "O", "Na", "C"};

class VTKWriterTest : public ::testing::Test {
protected:
    std::vector<std::string> tempfiles;

    std::string tempfile(const std::string &name)
    {
        tempfiles.push_back(name);
        return name;
    }

    void TearDown() override
    {
        for (const auto &f : tempfiles)
            remove(f.c_str());
    }
};

} // namespace

TEST_F(VTKWriterTest, LegacyAsciiPolyData)
{
    const auto file = tempfile("vtkwriter_poly.vtk");
    VTKWriter writer(VTKWriter::LEGACY, false);
    writer.set_title("test title");
    writer.set_polydata(COORDS);
    writer.add_point_array("id", 1, IDS);
    writer.add_point_array("v", 3, VECS);
    writer.write(file);

    const auto text = slurp(file);
    EXPECT_THAT(text, HasSubstr("# vtk DataFile Version 5.1\ntest title\nASCII\n"));
    EXPECT_THAT(text, HasSubstr("DATASET POLYDATA\nPOINTS 4 double\n"));

    // one vertex cell per point, offsets start at zero and have one extra entry
    EXPECT_THAT(text, HasSubstr("VERTICES 5 4\nOFFSETS vtktypeint64\n0 1 2 3 4 \n"));
    EXPECT_THAT(text, HasSubstr("CONNECTIVITY vtktypeint64\n0 1 2 3 \n"));

    // arrays that are not the active scalars become field data
    EXPECT_THAT(text, HasSubstr("POINT_DATA 4\nFIELD FieldData 2\n"));
    EXPECT_THAT(text, HasSubstr("id 1 4 int\n1 2 3 4 \n"));
    EXPECT_THAT(text, HasSubstr("v 3 4 double\n"));
    EXPECT_THAT(text, Not(HasSubstr("SCALARS")));

    EXPECT_EQ(parse_numbers(text.substr(text.find("POINTS 4 double\n") + 16)).size() >= 12, true);
}

TEST_F(VTKWriterTest, LegacyActiveScalars)
{
    const auto file = tempfile("vtkwriter_scalars.vtk");
    VTKWriter writer(VTKWriter::LEGACY, false);
    writer.set_polydata(COORDS);
    writer.add_point_array("intensity", 1, std::vector<double>{1.0, 2.0, 3.0, 4.0});
    writer.add_point_array("id", 1, IDS);
    writer.set_active_scalars("intensity");
    writer.write(file);

    const auto text = slurp(file);
    EXPECT_THAT(text, HasSubstr("SCALARS intensity double\nLOOKUP_TABLE default\n1 2 3 4 \n"));

    // the remaining array is still written as field data
    EXPECT_THAT(text, HasSubstr("FIELD FieldData 1\n"));
    EXPECT_THAT(text, HasSubstr("id 1 4 int\n"));
}

TEST_F(VTKWriterTest, LegacyBinaryIsBigEndian)
{
    const auto file = tempfile("vtkwriter_poly_bin.vtk");
    VTKWriter writer(VTKWriter::LEGACY, true);
    writer.set_polydata(COORDS);
    writer.add_point_array("id", 1, IDS);
    writer.write(file);

    const auto data = slurp(file);
    EXPECT_THAT(data, HasSubstr("BINARY\n"));

    auto pos = data.find("POINTS 4 double\n");
    ASSERT_NE(pos, std::string::npos);
    pos += strlen("POINTS 4 double\n");
    for (std::size_t i = 0; i < COORDS.size(); ++i)
        EXPECT_DOUBLE_EQ(be_double(data.data() + pos + 8 * i), COORDS[i]);

    pos = data.find("id 1 4 int\n");
    ASSERT_NE(pos, std::string::npos);
    pos += strlen("id 1 4 int\n");
    for (std::size_t i = 0; i < IDS.size(); ++i)
        EXPECT_EQ(be_int32(data.data() + pos + 4 * i), IDS[i]);
}

TEST_F(VTKWriterTest, LegacyStrings)
{
    const auto ascii = tempfile("vtkwriter_str.vtk");
    VTKWriter writer(VTKWriter::LEGACY, false);
    writer.set_polydata(COORDS);
    writer.add_point_array("element", NAMES);
    writer.write(ascii);

    // ASCII strings are written one per line
    EXPECT_THAT(slurp(ascii), HasSubstr("element 1 4 string\nC\nO\nNa\nC\n"));

    // binary strings carry a single byte 0xc0|length prefix
    const auto binary = tempfile("vtkwriter_str_bin.vtk");
    VTKWriter bwriter(VTKWriter::LEGACY, true);
    bwriter.set_polydata(COORDS);
    bwriter.add_point_array("element", NAMES);
    bwriter.write(binary);

    std::string expect("element 1 4 string\n");
    for (const auto &name : NAMES) {
        expect += static_cast<char>(0xc0 | name.size());
        expect += name;
    }
    expect += "\n";
    EXPECT_THAT(slurp(binary), HasSubstr(expect));
}

TEST_F(VTKWriterTest, LegacyBinaryRejectsLongStrings)
{
    const auto file = tempfile("vtkwriter_longstr.vtk");
    VTKWriter writer(VTKWriter::LEGACY, true);
    writer.set_polydata(COORDS);
    writer.add_point_array("element", {std::string(64, 'x'), "b", "c", "d"});
    EXPECT_THROW(writer.write(file), VTKWriterException);
}

TEST_F(VTKWriterTest, XmlAsciiPolyData)
{
    const auto file = tempfile("vtkwriter_poly.vtp");
    VTKWriter writer(VTKWriter::XML, false);
    writer.set_polydata(COORDS);
    writer.add_point_array("id", 1, IDS);
    writer.add_point_array("v", 3, VECS);
    writer.write(file);

    const auto text = slurp(file);
    EXPECT_THAT(text, HasSubstr(R"(<VTKFile type="PolyData" version="0.1")"));
    EXPECT_THAT(text, HasSubstr(R"(header_type="UInt32")"));
    EXPECT_THAT(text, HasSubstr(R"(<Piece NumberOfPoints="4" NumberOfVerts="4")"));
    EXPECT_THAT(text, HasSubstr(R"(<DataArray type="Int32" Name="id" format="ascii">)"));
    EXPECT_THAT(text, HasSubstr(R"(<DataArray type="Float64" Name="v" NumberOfComponents="3")"));
    EXPECT_THAT(text,
                HasSubstr(R"(<DataArray type="Float64" Name="Points" NumberOfComponents="3")"));
    EXPECT_THAT(text, HasSubstr(R"(<DataArray type="Int64" Name="connectivity")"));
    EXPECT_THAT(text, HasSubstr(R"(<DataArray type="Int64" Name="offsets")"));

    EXPECT_EQ(parse_numbers(payload_after(text, R"(Name="id")")),
              std::vector<double>({1, 2, 3, 4}));
    EXPECT_EQ(parse_numbers(payload_after(text, R"(Name="Points")")), COORDS);

    // XML offsets are end offsets, one per cell
    EXPECT_EQ(parse_numbers(payload_after(text, R"(Name="offsets")")),
              std::vector<double>({1, 2, 3, 4}));
    EXPECT_EQ(parse_numbers(payload_after(text, R"(Name="connectivity")")),
              std::vector<double>({0, 1, 2, 3}));
}

TEST_F(VTKWriterTest, XmlUnstructuredGridHasCellTypes)
{
    const auto file = tempfile("vtkwriter_grid.vtu");
    VTKWriter writer(VTKWriter::XML, false);
    writer.set_unstructured_grid(COORDS);
    writer.write(file);

    const auto text = slurp(file);
    EXPECT_THAT(text, HasSubstr(R"(<VTKFile type="UnstructuredGrid")"));
    EXPECT_THAT(text, HasSubstr(R"(<Piece NumberOfPoints="4" NumberOfCells="4">)"));
    // VTK_VERTEX is cell type 1
    EXPECT_EQ(parse_numbers(payload_after(text, R"(Name="types")")),
              std::vector<double>({1, 1, 1, 1}));
}

TEST_F(VTKWriterTest, Hexahedron)
{
    const auto file      = tempfile("vtkwriter_hex.vtu");
    double corners[8][3] = {{0, 0, 0}, {5, 0, 0}, {6, 4, 0}, {1, 4, 0},
                            {2, 1, 7}, {7, 1, 7}, {8, 5, 7}, {3, 5, 7}};
    VTKWriter writer(VTKWriter::XML, false);
    writer.set_hexahedron(corners);
    writer.write(file);

    const auto text = slurp(file);
    EXPECT_THAT(text, HasSubstr(R"(<Piece NumberOfPoints="8" NumberOfCells="1">)"));
    // VTK_HEXAHEDRON is cell type 12
    EXPECT_EQ(parse_numbers(payload_after(text, R"(Name="types")")), std::vector<double>({12}));
    EXPECT_EQ(parse_numbers(payload_after(text, R"(Name="offsets")")), std::vector<double>({8}));
    EXPECT_EQ(parse_numbers(payload_after(text, R"(Name="connectivity")")),
              std::vector<double>({0, 1, 2, 3, 4, 5, 6, 7}));
}

TEST_F(VTKWriterTest, XmlStringArray)
{
    const auto file = tempfile("vtkwriter_str.vtp");
    VTKWriter writer(VTKWriter::XML, false);
    writer.set_polydata(COORDS);
    writer.add_point_array("element", NAMES);
    writer.write(file);

    const auto text = slurp(file);
    // string arrays use <Array>, not <DataArray>
    EXPECT_THAT(text, HasSubstr(R"(<Array type="String" Name="element" format="ascii">)"));
    // the ASCII form lists character codes with a zero terminating each string
    EXPECT_EQ(parse_numbers(payload_after(text, R"(Name="element")")),
              std::vector<double>({67, 0, 79, 0, 78, 97, 0, 67, 0}));
}

TEST_F(VTKWriterTest, XmlBinaryPayloads)
{
    const auto file = tempfile("vtkwriter_poly_bin.vtp");
    VTKWriter writer(VTKWriter::XML, true);
    writer.set_polydata(COORDS);
    writer.add_point_array("id", 1, IDS);
    writer.write(file);

    const auto text = slurp(file);
    EXPECT_THAT(text, HasSubstr(R"(format="binary")"));

    // the test is compiled with zlib support whenever the writer is, so the
    // payload can always be decoded back to the original values

    const auto ids = decode_binary_payload(payload_after(text, R"(Name="id")"),
                                           IDS.size() * sizeof(std::int32_t));
    ASSERT_EQ(ids.size(), IDS.size() * sizeof(std::int32_t));
    for (std::size_t i = 0; i < IDS.size(); ++i) {
        std::int32_t v;
        memcpy(&v, ids.data() + 4 * i, 4);
        EXPECT_EQ(v, IDS[i]);
    }

    const auto pts = decode_binary_payload(payload_after(text, R"(Name="Points")"),
                                           COORDS.size() * sizeof(double));
    ASSERT_EQ(pts.size(), COORDS.size() * sizeof(double));
    for (std::size_t i = 0; i < COORDS.size(); ++i) {
        double v;
        memcpy(&v, pts.data() + 8 * i, 8);
        EXPECT_DOUBLE_EQ(v, COORDS[i]);
    }
}

TEST_F(VTKWriterTest, RectilinearGrid)
{
    const std::vector<double> xc = {0.0, 1.0, 3.0}, yc = {0.0, 2.0}, zc = {0.0, 0.5};
    const std::vector<double> cells = {1.0, 2.0};

    const auto file = tempfile("vtkwriter_rect.vtr");
    VTKWriter writer(VTKWriter::XML, false);
    writer.set_rectilinear_grid(xc, yc, zc);
    writer.add_cell_array("value", 1, cells);
    writer.write(file);

    const auto text = slurp(file);
    EXPECT_THAT(text, HasSubstr(R"(<VTKFile type="RectilinearGrid")"));
    EXPECT_THAT(text, HasSubstr(R"(<RectilinearGrid WholeExtent="0 2 0 1 0 1">)"));
    EXPECT_THAT(text, HasSubstr(R"(<Piece Extent="0 2 0 1 0 1">)"));
    EXPECT_EQ(parse_numbers(payload_after(text, R"(Name="x")")), xc);
    EXPECT_EQ(parse_numbers(payload_after(text, R"(Name="value")")), cells);

    // the same grid in the legacy format
    const auto legacy = tempfile("vtkwriter_rect.vtk");
    VTKWriter lwriter(VTKWriter::LEGACY, false);
    lwriter.set_rectilinear_grid(xc, yc, zc);
    lwriter.add_cell_array("value", 1, cells);
    lwriter.write(legacy);

    const auto ltext = slurp(legacy);
    EXPECT_THAT(ltext, HasSubstr("DATASET RECTILINEAR_GRID\nDIMENSIONS 3 2 2\n"));
    EXPECT_THAT(ltext, HasSubstr("X_COORDINATES 3 double\n0 1 3 \n"));
    EXPECT_THAT(ltext, HasSubstr("CELL_DATA 2\n"));
}

TEST_F(VTKWriterTest, ImageData)
{
    const int dim[3]        = {3, 2, 2};
    const double origin[3]  = {-1.0, -2.0, -3.0};
    const double spacing[3] = {0.25, 0.5, 1.0};
    std::vector<double> vals(12);
    for (int i = 0; i < 12; ++i)
        vals[i] = 0.5 * i;

    const auto file = tempfile("vtkwriter_image.vti");
    VTKWriter writer(VTKWriter::XML, false);
    writer.set_image_data(dim, origin, spacing);
    writer.add_point_array("intensity", 1, vals);
    writer.set_active_scalars("intensity");
    writer.write(file);

    const auto text = slurp(file);
    EXPECT_THAT(text, HasSubstr(R"(<VTKFile type="ImageData")"));
    EXPECT_THAT(text, HasSubstr(R"(Origin="-1 -2 -3" Spacing="0.25 0.5 1")"));
    EXPECT_THAT(text, HasSubstr(R"(<PointData Scalars="intensity">)"));
    EXPECT_EQ(parse_numbers(payload_after(text, R"(Name="intensity")")), vals);

    // the legacy spelling of the same data set
    const auto legacy = tempfile("vtkwriter_image.vtk");
    VTKWriter lwriter(VTKWriter::LEGACY, false);
    lwriter.set_image_data(dim, origin, spacing);
    lwriter.add_point_array("intensity", 1, vals);
    lwriter.set_active_scalars("intensity");
    lwriter.write(legacy);

    const auto ltext = slurp(legacy);
    EXPECT_THAT(ltext, HasSubstr("DATASET STRUCTURED_POINTS\nDIMENSIONS 3 2 2\n"));
    EXPECT_THAT(ltext, HasSubstr("SPACING 0.25 0.5 1\nORIGIN -1 -2 -3\n"));
    EXPECT_THAT(ltext, HasSubstr("POINT_DATA 12\nSCALARS intensity double\n"));
}

TEST_F(VTKWriterTest, ErrorsOnInconsistentInput)
{
    // no dataset selected
    {
        VTKWriter writer(VTKWriter::XML, false);
        EXPECT_THROW(writer.add_point_array("id", 1, IDS), VTKWriterException);
        EXPECT_THROW(writer.write("vtkwriter_never.vtp"), VTKWriterException);
    }
    // two datasets
    {
        VTKWriter writer(VTKWriter::XML, false);
        writer.set_polydata(COORDS);
        EXPECT_THROW(writer.set_unstructured_grid(COORDS), VTKWriterException);
    }
    // coordinate list that is not a multiple of 3
    {
        VTKWriter writer(VTKWriter::XML, false);
        EXPECT_THROW(writer.set_polydata({0.0, 1.0}), VTKWriterException);
    }
    // array length does not match the number of points
    {
        VTKWriter writer(VTKWriter::XML, false);
        writer.set_polydata(COORDS);
        EXPECT_THROW(writer.add_point_array("id", 1, std::vector<int>({1, 2})), VTKWriterException);
        EXPECT_THROW(writer.add_point_array("v", 3, std::vector<double>({1.0, 2.0, 3.0})),
                     VTKWriterException);
    }
    // marking an array that was never added
    {
        VTKWriter writer(VTKWriter::XML, false);
        writer.set_polydata(COORDS);
        EXPECT_THROW(writer.set_active_scalars("nope"), VTKWriterException);
    }
    // a location that cannot be written
    {
        VTKWriter writer(VTKWriter::XML, false);
        writer.set_polydata(COORDS);
        EXPECT_THROW(writer.write("no_such_directory/out.vtp"), VTKWriterException);
    }
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleMock(&argc, argv);
    return RUN_ALL_TESTS();
}
