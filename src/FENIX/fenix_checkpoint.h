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

#ifndef LMP_FENIX_CHECKPOINT_H
#define LMP_FENIX_CHECKPOINT_H

#include "error.h"
#include "write_restart.h"
#include "FENIX/local_serializer.h"

#include <string>

namespace LAMMPS_NS {

class Fenix;

class FenixCheckpoint : public WriteRestart {
 public:
  FenixCheckpoint(class LAMMPS *);

  void command(int, char**) override {
    error->all(FLERR, "fenix_checkpoint is not an invokable command!");
  };

  void write(const std::string &) final;

  void checkpoint();
  void recover();

  ~FenixCheckpoint();

 protected:
  friend Fenix;

  void create_group();

  void store_data();
  void store_metadata();

  void commit();

  void restore_data();
  void restore_metadata();

 private:
  void safe_create_member(int, void*, int, MPI_Datatype);
  void safe_store(int, int);
  void safe_storev(int, int);

  void safe_restore(int, void*, int);
  int safe_restore(int, char*&);

  int group_id;
  int data_member;
  int metadata_member;

  int group_policy;
  int* group_policy_args;

  LocalSerializer serializer {lmp};
};

}

#endif
