/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef TANGENTIAL_CONTACT_MODELS_H_
#define TANGENTIAL_CONTACT_MODELS_H_

#include "contact.h";
#include "sub_model.h"

namespace Contact {

class TangentialModel:SubModel
{
public:
  TangentialModel() {};
  virtual ~TangentialModel() {};
  virtual double calculate_forces() = 0;
  virtual void coeffs_to_local();
  virtual void mix_coeffs(TangentialModel*, TangentialModel*);
  int rescale_flag;
private:
  int beyond_contact;
  int allow_limit_damping;
};

/* ---------------------------------------------------------------------- */

class TangentialLinearNoHistory:TangentialModel
{
public:
  TangentialLinearNoHistory()
  virtual void coeffs_to_local();
  virtual void mix_coeffs(TangentialModel*, TangentialModel*);
  double calculate_forces();
private:
  double xt, damp, mu;
};

/* ---------------------------------------------------------------------- */

class TangentialLinearHistory:TangentialModel
{
public:
  TangentialLinearHistory()
  virtual void coeffs_to_local();
  virtual void mix_coeffs(TangentialModel*, TangentialModel*);
  double calculate_forces();
private:
  double kt, xt, damp, mu;
};

/* ---------------------------------------------------------------------- */

class TangentialMindlin:TangentialModel
{
public:
  TangentialMindlin();
  void coeffs_to_local();
  void coeffs_to_local(TangentialModel*, TangentialModel*);
  double calculate_forces();
private:
  double k_norm, damp, Emod, poiss, coh;
  int mindlin_rescale, mindlin_force;
};

/* ---------------------------------------------------------------------- */

class TangentialMindlinForce:TangentialMindlin
{
public:
  TangentialMindlinForce();
};

/* ---------------------------------------------------------------------- */

class TangentialMindlinRescale:TangentialMindlin
{
public:
  TangentialMindlinRescale();
};

/* ---------------------------------------------------------------------- */

class TangentialMindlinRescaleForce:TangentialMindlinRescale
{
public:
  TangentialMindlinForceRescale();
};
}

#endif /*TANGENTIAL_CONTACT_MODELS_H_ */

