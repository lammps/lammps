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
   Contributing authors:  Axel Kohlmeyer (Temple U),
------------------------------------------------------------------------- */

#include "geturl.h"

#include "comm.h"
#include "error.h"

#if defined(LAMMPS_CURL)
#include <curl/curl.h>
#endif

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

void GetURL::command(int narg, char **arg)
{
#if !defined(LAMMPS_CURL)
  error->all(FLERR, "LAMMPS has not been compiled with libcurl support");
#else
  if (narg < 1) utils::missing_cmd_args(FLERR, "geturl", error);
  int verify = 1;
  int overwrite = 1;
  int verbose = 0;

  // process arguments

  std::string url = arg[0];

  // sanity check

  if ((url.find(':') == std::string::npos) || (url.find('/') == std::string::npos))
    error->all(FLERR, "URL '{}' is not a supported URL", url);

  std::string output = url.substr(url.find_last_of('/') + 1);
  if (output.empty()) error->all(FLERR, "URL '{}' must end in a file string", url);

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "output") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "geturl output", error);
      output = arg[iarg + 1];
      ++iarg;
    } else if (strcmp(arg[iarg], "overwrite") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "geturl overwrite", error);
      overwrite = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      ++iarg;
    } else if (strcmp(arg[iarg], "verify") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "geturl verify", error);
      verify = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      ++iarg;
    } else if (strcmp(arg[iarg], "verbose") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "geturl verbose", error);
      verbose = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      ++iarg;
    } else {
      error->all(FLERR, "Unknown geturl keyword: {}", arg[iarg]);
    }
    ++iarg;
  }

  // only download files from rank 0

  if (comm->me != 0) return;

  if (!overwrite && platform::file_is_readable(output)) return;

  // open output file for writing

  FILE *out = fopen(output.c_str(), "wb");
  if (!out)
    error->all(FLERR, "Cannot open output file {} for writing: {}", output, utils::getsyserror());

  // initialize curl and perform download

  CURL *curl;
  curl_global_init(CURL_GLOBAL_DEFAULT);
  curl = curl_easy_init();
  if (curl) {
    (void) curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    (void) curl_easy_setopt(curl, CURLOPT_WRITEDATA, (void *) out);
    (void) curl_easy_setopt(curl, CURLOPT_FILETIME, 1L);
    (void) curl_easy_setopt(curl, CURLOPT_FAILONERROR, 1L);
    if (verbose && screen) {
      (void) curl_easy_setopt(curl, CURLOPT_VERBOSE, 1L);
      (void) curl_easy_setopt(curl, CURLOPT_STDERR, (void *) screen);
    }
    if (!verify) {
      (void) curl_easy_setopt(curl, CURLOPT_SSL_VERIFYPEER, 0L);
      (void) curl_easy_setopt(curl, CURLOPT_SSL_VERIFYHOST, 0L);
    }
    auto res = curl_easy_perform(curl);
    if (res != CURLE_OK) {
      long response = 0L;
      curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &response);
      error->one(FLERR, "Download of {} failed with: {} {}", output, curl_easy_strerror(res),
                 response);
    }
    curl_easy_cleanup(curl);
  }
  curl_global_cleanup();
  fclose(out);
#endif
}
