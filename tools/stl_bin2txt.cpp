/* -----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------ */

// Convert a binary STL file to ASCII format
// Contributing author: Axel Kohlmeyer, Temple U, akohlmey at gmail.com
//
// Specs of the format taken from: https://en.wikipedia.org/wiki/STL_(file_format)

#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstring>

int main(int argc, char **argv)
{
    FILE *in, *out;
    char title[80];
    float normal[3], vert1[3], vert2[3], vert3[3];
    uint32_t ntriangles;
    size_t count;
    uint16_t attributes;

    if (argc != 3) {
        printf("Usage: %s <input file>  <output file>\n", argv[0]);
        return 1;
    }

    in = fopen(argv[1], "rb");
    if (!in) {
        printf("Error opening input file %s: %s\n", argv[1], strerror(errno));
        return 2;
    }
    out = fopen(argv[2], "w");
    if (!out) {
        printf("Error opening output file %s: %s\n", argv[1], strerror(errno));
        return 3;
    }

    /* read header */
    count = fread(title, sizeof(char), sizeof(title), in);
    if (count != sizeof(title)) {
        printf("Error reading binary STL header: %s\n", strerror(errno));
        return 4;
    }
    title[79] = '\0'; // ensure null termination of title string
    count = strlen(title);
    if (count == 0) snprintf(title, 80, "STL object from file %s", argv[1]);

    /* read triangle count */
    count = fread(&ntriangles, sizeof(uint32_t), 1, in);
    if (count != 1) {
        printf("Error reading binary STL triangle count: %s\n", strerror(errno));
        return 5;
    }

    /* write header */
    printf("Converting: %s with %u triangles\n", title, ntriangles);
    fprintf(out, "solid %s\n", title);

    /* loop over triangles */
    for (uint32_t i = 0; i < ntriangles; ++i) {
        count = fread(normal, sizeof(float), 3, in);
        count += fread(vert1, sizeof(float), 3, in);
        count += fread(vert2, sizeof(float), 3, in);
        count += fread(vert3, sizeof(float), 3, in);
        if (count != 12) {
            printf("Error reading binary STL vertices: %s\n", strerror(errno));
            return 6;
        }
        count = fread(&attributes, sizeof(uint16_t), 1, in);
        if (count != 1) {
            printf("Error reading binary STL facet attributes: %s\n", strerror(errno));
            return 7;
        }
        fprintf(out, "  facet normal %e %e %e\n", normal[0], normal[1], normal[2]);
        fputs("    outer loop\n", out);
        fprintf(out, "    vertex %e %e %e\n", vert1[0], vert1[1], vert1[2]);
        fprintf(out, "    vertex %e %e %e\n", vert2[0], vert2[1], vert2[2]);
        fprintf(out, "    vertex %e %e %e\n", vert3[0], vert3[1], vert3[2]);
        fputs("    endloop\n  endfacet\n", out);
        if (ferror(out)) {
            printf("Error writing text STL facet: %s\n", strerror(errno));
            return 7;
        }
    }
    fprintf(out, "endsolid %s\n", title);
    fclose(in);
    fclose(out);
    return 0;
}
