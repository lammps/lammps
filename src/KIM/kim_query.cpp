/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Axel Kohlmeyer (Temple U),
                         Ryan S. Elliott (UMN)
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the Free
   Software Foundation; either version 2 of the License, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
   more details.

   You should have received a copy of the GNU General Public License along with
   this program; if not, see <https://www.gnu.org/licenses>.

   Linking LAMMPS statically or dynamically with other modules is making a
   combined work based on LAMMPS. Thus, the terms and conditions of the GNU
   General Public License cover the whole combination.

   In addition, as a special exception, the copyright holders of LAMMPS give
   you permission to combine LAMMPS with free software programs or libraries
   that are released under the GNU LGPL and with code included in the standard
   release of the "kim-api" under the CDDL (or modified versions of such code,
   with unchanged license). You may copy and distribute such a system following
   the terms of the GNU GPL for LAMMPS and the licenses of the other code
   concerned, provided that you include the source code of that other code
   when and as the GNU GPL requires distribution of source code.

   Note that people who make modified versions of LAMMPS are not obligated to
   grant this special exception for their modified versions; it is their choice
   whether to do so. The GNU General Public License gives permission to release
   a modified version without this exception; this exception also makes it
   possible to release a modified version which carries forward this exception.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Designed for use with the kim-api-2.1.0 (and newer) package
------------------------------------------------------------------------- */

#include "kim_query.h"
#include <mpi.h>
#include <cstring>
#include <string>
#include <sstream>
#include "comm.h"
#include "error.h"
#include "input.h"
#include "modify.h"
#include "variable.h"
#include "fix_store_kim.h"

#if defined(LMP_KIM_CURL)
#include <sys/types.h>
#include <curl/curl.h>
#endif

using namespace LAMMPS_NS;

#if defined(LMP_KIM_CURL)

struct WriteBuf {
  char *dataptr;
  size_t sizeleft;
};

static char *do_query(char *, char *, int, char **, int, MPI_Comm);
static size_t write_callback(void *, size_t, size_t, void *);

#endif

/* ---------------------------------------------------------------------- */

void KimQuery::command(int narg, char **arg)
{
  char *varname, *function, *value;

  if (narg < 2) error->all(FLERR,"Illegal kim_query command");

  // check if we had a kim_init command by finding fix STORE/KIM
  // retrieve model name.
  char * model_name;

  int ifix = modify->find_fix("KIM_MODEL_STORE");
  if (ifix >= 0) {
    FixStoreKIM *fix_store = (FixStoreKIM *) modify->fix[ifix];
    model_name = (char *)fix_store->getptr("model_name");
  } else error->all(FLERR,"Must use 'kim_init' before 'kim_query'");


  varname = arg[0];
  bool split = false;
  if (0 == strcmp("split",arg[1])) {
    if (narg == 2) error->all(FLERR,"Illegal kim_query command");
    split = true;
    arg++;
    narg--;
  }
  function = arg[1];
  for (int i = 2; i < narg; ++i)
  {
    if (0 == strncmp("model=",arg[i], 6)) {
      error->all(FLERR,"Illegal 'model' key in kim_query command");
    }
  }


#if defined(LMP_KIM_CURL)

  value = do_query(function, model_name, narg-2, arg+2, comm->me, world);

  // check for valid result
  // on error the content of "value" is a '\0' byte
  // as the first element, and then the error message
  // that was returned by the web server

  char errmsg[1024];
  if (0 == strlen(value)) {
    sprintf(errmsg,"OpenKIM query failed: %s",value+1);
    error->all(FLERR,errmsg);
  } else if (0 == strcmp(value,"EMPTY")) {
    sprintf(errmsg,"OpenKIM query returned no results");
    error->all(FLERR,errmsg);
  }

  kim_query_log_delimiter("begin");
  char **varcmd = new char*[3];
  if (split) {
    int counter = 1;
    std::stringstream ss(value);
    std::string token;
    varcmd[1] = (char *) "string";

    while(std::getline(ss, token, ',')) {
      token.erase(0,token.find_first_not_of(" \n\r\t"));  // ltrim
      token.erase(token.find_last_not_of(" \n\r\t") + 1);  // rtrim
      std::stringstream splitname;
      splitname << varname << "_" << counter++;
      varcmd[0] = const_cast<char *>(splitname.str().c_str());
      varcmd[2] = const_cast<char *>(token.c_str());
      input->variable->set(3,varcmd);
      echo_var_assign(splitname.str(), varcmd[2]);
    }
  }
  else {
    varcmd[0] = varname;
    varcmd[1] = (char *) "string";
    varcmd[2] = value;
    input->variable->set(3,varcmd);

    echo_var_assign(varname, value);
  }
  kim_query_log_delimiter("end");

  delete[] varcmd;
  delete[] value;
#else
  error->all(FLERR,"Cannot use 'kim_query' command when KIM package "
             "is compiled without support for libcurl");
#endif
}

#if defined(LMP_KIM_CURL)

// copy data to the user provided data structure, optionally in increments

size_t write_callback(void *data, size_t size, size_t nmemb, void *userp)
{
  struct WriteBuf *buf = (struct WriteBuf *)userp;
  size_t buffer_size = size*nmemb;

  // copy chunks into the buffer for as long as there is space left
  if (buf->sizeleft) {
    size_t copy_this_much = buf->sizeleft;
    if (copy_this_much > buffer_size)
      copy_this_much = buffer_size;
    memcpy(buf->dataptr, data, copy_this_much);

    buf->dataptr += copy_this_much;
    buf->sizeleft -= copy_this_much;
    return copy_this_much;
  }
  return 0; // done
}

char *do_query(char *qfunction, char * model_name, int narg, char **arg,
               int rank, MPI_Comm comm)
{
  char value[512], *retval;

  // run the web query from rank 0 only

  if (rank == 0) {
    CURL *handle;
    CURLcode res;

    // set up and clear receive buffer

    struct WriteBuf buf;
    buf.dataptr = value;
    buf.sizeleft = 511;
    memset(value,0,512);

    // create curl web query instance
    curl_global_init(CURL_GLOBAL_DEFAULT);
    handle = curl_easy_init();

    if (handle) {
      std::string url("https://query.openkim.org/api/");
      url += qfunction;

      std::string query(arg[0]);
      query += "&model=[\"";
      query += model_name;
      query += "\"]";
      for (int i=1; i < narg; ++i) {
        query += '&';
        query += arg[i];
      }

#if LMP_DEBUG_CURL
      curl_easy_setopt(handle, CURLOPT_VERBOSE, 1L);
#endif

#if defined(LMP_NO_SSL_CHECK)
      // disable verifying SSL certificate and host name
      curl_easy_setopt(handle, CURLOPT_SSL_VERIFYPEER, 0L);
      curl_easy_setopt(handle, CURLOPT_SSL_VERIFYHOST, 0L);
#endif

      curl_easy_setopt(handle, CURLOPT_URL, url.c_str());
      curl_easy_setopt(handle, CURLOPT_FOLLOWLOCATION, 1L);
      curl_easy_setopt(handle, CURLOPT_POSTFIELDS, query.c_str());

      curl_easy_setopt(handle, CURLOPT_WRITEFUNCTION,write_callback);
      curl_easy_setopt(handle, CURLOPT_WRITEDATA,&buf);

      // perform OpenKIM query and check for errors
      res = curl_easy_perform(handle);
      if (res != CURLE_OK) {
        // on error we return an "empty" string but add error message after it
        value[0]= '\0';
        strcpy(value+1,curl_easy_strerror(res));
      }
      curl_easy_cleanup(handle);
    }
    curl_global_cleanup();
  }
  MPI_Bcast(value, 512, MPI_CHAR, 0, comm);

  // we must make a proper copy of the query, as the stack allocation
  // for "value" will go out of scope. a valid query has a '[' as
  // the first character. skip over it (and the last character ']', too)
  // an error messages starts with a '\0' character. copy that and
  // the remaining string, as that is the error message.

  if (value[0] == '[') {
    int len = strlen(value)-1;
    if (value[len] == ']') {
      retval = new char[len];
      value[len] = '\0';
      if (0 == strcmp(value+1, "")) {
        strcpy(retval,"EMPTY");
      }
      else
        strcpy(retval,value+1);
    } else {
      retval = new char[len+2];
      retval[0] = '\0';
      strcpy(retval+1,value);
    }
  } else if (value[0] == '\0') {
    int len = strlen(value+1)+2;
    retval = new char[len];
    retval[0] = '\0';
    strcpy(retval+1,value+1);
  } else {
    // unknown response type. we should not get here.
    // we return an "empty" string but add error message after it
    int len = strlen(value)+2;
    retval = new char[len];
    retval[0] = '\0';
    strcpy(retval+1,value);
  }
  return retval;
}
#endif

/* ---------------------------------------------------------------------- */

void KimQuery::kim_query_log_delimiter(std::string const begin_end) const
{
  if (comm->me == 0) {
    std::string mesg;
    if (begin_end == "begin")
      mesg =
          "#=== BEGIN kim-query =========================================\n";
    else if (begin_end == "end")
      mesg =
          "#=== END kim-query ===========================================\n\n";

    input->write_echo(mesg.c_str());
  }
}

/* ---------------------------------------------------------------------- */

void KimQuery::echo_var_assign(std::string const & name,
                               std::string const & value) const
{
  if (comm->me == 0) {
    std::string mesg;
    mesg += "variable " + name + " string " + value + "\n";
    input->write_echo(mesg.c_str());
  }
}
