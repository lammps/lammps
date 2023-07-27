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

#ifndef LAMMPSRUNNER_H
#define LAMMPSRUNNER_H

#include <QThread>

class LammpsRunner : public QThread {
    Q_OBJECT

public:
    LammpsRunner(QObject *parent = nullptr);
    ~LammpsRunner() = default;

    // execute LAMMPS run in runner thread
    void run() override;

    // transfer info to worker thread
    void setup_run(void *_handle, const char *_input, void *_plugin = nullptr)
    {
        handle = _handle;
        plugin = _plugin;
        input  = _input;
    }

signals:
    void resultReady();

private:
    void *handle, *plugin;
    const char *input;
};

#endif
// Local Variables:
// c-basic-offset: 4
// End:
