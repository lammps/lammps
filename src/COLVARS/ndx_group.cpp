// clang-format off
// -*- c++ -*-

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
   Contributing author:  Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "ndx_group.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "group.h"
#include "tokenizer.h"

using namespace LAMMPS_NS;
#define BUFLEN 4096
#define DELTA 16384

// read file until next section "name" or any next section if name == ""

static std::string find_section(FILE *fp, const std::string &name)
{
  char linebuf[BUFLEN];

  std::string pattern = "^\\s*\\[\\s+\\S+\\s+\\]\\s*$";
  if (!name.empty())
    pattern = fmt::format("^\\s*\\[\\s+{}\\s+\\]\\s*$",name);

  fgets(linebuf,BUFLEN,fp);
  while (!feof(fp)) {
    if (utils::strmatch(linebuf,pattern))
      return Tokenizer(linebuf).as_vector()[1];
    fgets(linebuf,BUFLEN,fp);
  }
  return "";
}

static std::vector<tagint> read_section(FILE *fp, std::string &name)
{
  char linebuf[BUFLEN];
  std::vector<tagint> tagbuf;
  std::string pattern = "^\\s*\\[\\s+\\S+\\s+\\]\\s*$";

  while (fgets(linebuf,BUFLEN,fp)) {
    // start of new section. we are done, update "name"
    if (utils::strmatch(linebuf,pattern)) {
      name = Tokenizer(linebuf).as_vector()[1];
      return tagbuf;
    }
    ValueTokenizer values(linebuf);
    while (values.has_next())
      tagbuf.push_back(values.next_tagint());
  }
  // set empty name to indicate end of file
  name = "";
  return tagbuf;
}

/* ---------------------------------------------------------------------- */

void Ndx2Group::command(int narg, char **arg)
{
  int len;
  bigint num;
  FILE *fp;
  std::string name, next;

  if (narg < 1) error->all(FLERR,"Illegal ndx2group command");
  if (atom->tag_enable == 0)
    error->all(FLERR,"Must have atom IDs for ndx2group command");
  if (atom->map_style == Atom::MAP_NONE)
    error->all(FLERR,"Must have an atom map for ndx2group command");
  if (comm->me == 0) {
    fp = fopen(arg[0], "r");
    if (fp == nullptr)
      error->one(FLERR,"Cannot open index file for reading: {}",
                                   utils::getsyserror());
    utils::logmesg(lmp,"Reading groups from index file {}:\n",arg[0]);
  }

  if (narg == 1) {              // restore all groups

    if (comm->me == 0) {
      name = find_section(fp,"");
      while (!name.empty()) {

        // skip over group "all", which is called "System" in gromacs
        if (name == "System") {
          name = find_section(fp,"");
          continue;
        }

        utils::logmesg(lmp," Processing group '{}'\n",name);
        len = name.size()+1;
        MPI_Bcast(&len,1,MPI_INT,0,world);
        if (len > 1) {
          MPI_Bcast((void *)name.c_str(),len,MPI_CHAR,0,world);

          // read tags for atoms in group and broadcast
          std::vector<tagint> tags = read_section(fp,next);
          num = tags.size();
          MPI_Bcast(&num,1,MPI_LMP_BIGINT,0,world);
          MPI_Bcast((void *)tags.data(),num,MPI_LMP_TAGINT,0,world);
          create(name,tags);
          name = next;
        }
      }
      len = -1;
      MPI_Bcast(&len,1,MPI_INT,0,world);

    } else {

      while (true) {
        MPI_Bcast(&len,1,MPI_INT,0,world);
        if (len < 0) break;
        if (len > 1) {
          char *buf = new char[len];
          MPI_Bcast(buf,len,MPI_CHAR,0,world);
          MPI_Bcast(&num,1,MPI_LMP_BIGINT,0,world);
          tagint *tbuf = new tagint[num];
          MPI_Bcast(tbuf,num,MPI_LMP_TAGINT,0,world);
          create(buf,std::vector<tagint>(tbuf,tbuf+num));
          delete[] buf;
          delete[] tbuf;
        }
      }
    }

  } else {            // restore selected groups

    for (int idx=1; idx < narg; ++idx) {
      if (comm->me == 0) {

        // find named section, search from beginning of file
        rewind(fp);
        name = find_section(fp,arg[idx]);
        utils::logmesg(lmp," {} group '{}'\n", name.size()
                       ? "Processing" : "Skipping",arg[idx]);
        len = name.size()+1;
        MPI_Bcast(&len,1,MPI_INT,0,world);
        if (len > 1) {
          MPI_Bcast((void *)name.c_str(),len,MPI_CHAR,0,world);

          // read tags for atoms in group and broadcast
          std::vector<tagint> tags = read_section(fp,next);
          num = tags.size();
          MPI_Bcast(&num,1,MPI_LMP_BIGINT,0,world);
          MPI_Bcast((void *)tags.data(),num,MPI_LMP_TAGINT,0,world);
          create(name,tags);
          name = next;
        }
      } else {
        MPI_Bcast(&len,1,MPI_INT,0,world);
        if (len > 1) {
          char *buf = new char[len];
          MPI_Bcast(buf,len,MPI_CHAR,0,world);
          MPI_Bcast(&num,1,MPI_LMP_BIGINT,0,world);
          tagint *tbuf = new tagint[num];
          MPI_Bcast(tbuf,num,MPI_LMP_TAGINT,0,world);
          create(buf,std::vector<tagint>(tbuf,tbuf+num));
          delete[] buf;
          delete[] tbuf;
        }
      }
    }
  }
  if (comm->me == 0) fclose(fp);
}

/* ---------------------------------------------------------------------- */

void Ndx2Group::create(const std::string &name, const std::vector<tagint> &tags)
{
  // wipe out all members if the group exists. gid==0 is group "all"
  int gid = group->find(name);
  if (gid > 0) group->assign(name + " clear");

  // map from global to local
  const int nlocal = atom->nlocal;
  int *flags = (int *)calloc(nlocal,sizeof(int));
  for (bigint i=0; i < (int)tags.size(); ++i) {
    const int id = atom->map(tags[i]);
    if (id < nlocal && id >= 0) flags[id] = 1;
  }
  group->create(name,flags);
  free(flags);
}
