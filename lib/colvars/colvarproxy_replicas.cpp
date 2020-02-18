// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include "colvarmodule.h"
#include "colvarproxy.h"


colvarproxy_replicas::colvarproxy_replicas() {}


colvarproxy_replicas::~colvarproxy_replicas() {}


int colvarproxy_replicas::replica_enabled()
{
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy_replicas::replica_index()
{
  return 0;
}


int colvarproxy_replicas::num_replicas()
{
  return 1;
}


void colvarproxy_replicas::replica_comm_barrier() {}


int colvarproxy_replicas::replica_comm_recv(char* /* msg_data */,
                                            int /* buf_len */,
                                            int /* src_rep */)
{
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy_replicas::replica_comm_send(char* /* msg_data */,
                                            int /* msg_len */,
                                            int /* dest_rep */)
{
  return COLVARS_NOT_IMPLEMENTED;
}


