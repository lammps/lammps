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

#ifndef STDCAPTURE_H
#define STDCAPTURE_H

#include <string>

class StdCapture {
public:
    StdCapture();
    virtual ~StdCapture();

    void BeginCapture();
    bool EndCapture();
    std::string GetCapture();
    std::string GetChunk();

private:
    enum PIPES { READ, WRITE };
    int m_pipe[2];
    int m_oldStdOut;
    bool m_capturing;
    std::string m_captured;

    static constexpr int bufSize = 1025;
    char buf[bufSize];
};

#endif
// Local Variables:
// c-basic-offset: 4
// End:
