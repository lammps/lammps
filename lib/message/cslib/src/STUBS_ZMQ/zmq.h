/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// ZMQ constants and dummy functions

#ifndef ZMQ_DUMMY_H
#define ZMQ_DUMMY_H

namespace CSLIB_NS {

#define ZMQ_REQ 0
#define ZMQ_REP 0

static void *zmq_ctx_new() {return NULL;}
static void *zmq_connect(void *, char *) {return NULL;}
static int zmq_bind(void *, char *) {return 0;}
static void *zmq_socket(void *,int) {return NULL;}
static void zmq_close(void *) {}
static void zmq_ctx_destroy(void *) {}
static void zmq_send(void *, void *, int, int) {}
static void zmq_recv(void *, void *, int, int) {}

};

#endif
