
/* ----------------------------------------------------------------------
   gather the named atom-based entity for all atoms
     return it in user-allocated data
   data will be ordered by atom ID
     requirement for consecutive atom IDs (1 to N)
   see gather_atoms_concat() to return data for all atoms, unordered
   see gather_atoms_subset() to return data for only a subset of atoms
   name = desired quantity, e.g. x or charge
   type = 0 for integer values, 1 for double values
   count = # of per-atom values, e.g. 1 for type or charge, 3 for x or f
     use count = 3 with "image" if want single image flag unpacked into xyz
   return atom-based values in 1d data, ordered by count, then by atom ID
     e.g. x[0][0],x[0][1],x[0][2],x[1][0],x[1][1],x[1][2],x[2][0],...
     data must be pre-allocated by caller to correct length
     correct length = count*Natoms, as queried by get_natoms()
   method:
     alloc and zero count*Natom length vector
     loop over Nlocal to fill vector with my values
     Allreduce to sum vector into data across all procs
------------------------------------------------------------------------- */

#if defined(LAMMPS_BIGBIG)
void lammps_gather_atoms(void *ptr, char * /*name */,
                         int /*type*/, int /*count*/, void * /*data*/)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  BEGIN_CAPTURE
  {
    lmp->error->all(FLERR,"Library function lammps_gather_atoms() "
                    "is not compatible with -DLAMMPS_BIGBIG");
  }
  END_CAPTURE
}
#else
void lammps_gather_atoms(void *ptr, char *name,
                         int type, int count, void *data)
{
  LAMMPS *lmp = (LAMMPS *) ptr;
  BEGIN_CAPTURE
  {
    int i,j,offset;

    // error if tags are not defined or not consecutive
    // NOTE: test that name = image or ids is not a 64-bit int in code?

    int flag = 0;
    if (lmp->atom->tag_enable == 0 || lmp->atom->tag_consecutive() == 0)
      flag = 1;
    if (lmp->atom->natoms > MAXSMALLINT) flag = 1;
    if (flag) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR,"Library error in lammps_gather_atoms");
      return;
    }

    int natoms = static_cast<int> (lmp->atom->natoms);

    void *vptr = lammps_extract_atom(lmp,name);

    if (vptr == NULL) {
      lmp->error->warning(FLERR,"lammps_gather_atoms: unknown property name");
      return;
    }

    // copy = Natom length vector of per-atom values
    // use atom ID to insert each atom's values into copy
    // MPI_Allreduce with MPI_SUM to merge into data, ordered by atom ID

    if (type == 0) {
      int *vector = NULL;
      int **array = NULL;
      const int imgunpack = (count == 3) && (strcmp(name,"image") == 0);

      if ((count == 1) || imgunpack) vector = (int *) vptr;
      else array = (int **) vptr;

      int *copy;
      lmp->memory->create(copy,count*natoms,"lib/gather:copy");
      for (i = 0; i < count*natoms; i++) copy[i] = 0;

      tagint *tag = lmp->atom->tag;
      int nlocal = lmp->atom->nlocal;

      if (count == 1) {
        for (i = 0; i < nlocal; i++)
          copy[tag[i]-1] = vector[i];

      } else if (imgunpack) {
        for (i = 0; i < nlocal; i++) {
          offset = count*(tag[i]-1);
          const int image = vector[i];
          copy[offset++] = (image & IMGMASK) - IMGMAX;
          copy[offset++] = ((image >> IMGBITS) & IMGMASK) - IMGMAX;
          copy[offset++] = ((image >> IMG2BITS) & IMGMASK) - IMGMAX;
        }

      } else {
        for (i = 0; i < nlocal; i++) {
          offset = count*(tag[i]-1);
          for (j = 0; j < count; j++)
            copy[offset++] = array[i][j];
        }
      }

      MPI_Allreduce(copy,data,count*natoms,MPI_INT,MPI_SUM,lmp->world);
      lmp->memory->destroy(copy);

    } else {
      double *vector = NULL;
      double **array = NULL;
      if (count == 1) vector = (double *) vptr;
      else array = (double **) vptr;

      double *copy;
      lmp->memory->create(copy,count*natoms,"lib/gather:copy");
      for (i = 0; i < count*natoms; i++) copy[i] = 0.0;

      tagint *tag = lmp->atom->tag;
      int nlocal = lmp->atom->nlocal;

      if (count == 1) {
        for (i = 0; i < nlocal; i++)
          copy[tag[i]-1] = vector[i];

      } else {
        for (i = 0; i < nlocal; i++) {
          offset = count*(tag[i]-1);
          for (j = 0; j < count; j++)
            copy[offset++] = array[i][j];
        }
      }

      MPI_Allreduce(copy,data,count*natoms,MPI_DOUBLE,MPI_SUM,lmp->world);
      lmp->memory->destroy(copy);
    }
  }
  END_CAPTURE
}
#endif

/* ----------------------------------------------------------------------
   gather the named atom-based entity for all atoms
     return it in user-allocated data
   data will be a concatenation of chunks of each proc's atoms,
     in whatever order the atoms are on each proc
     no requirement for consecutive atom IDs (1 to N)
     can do a gather_atoms_concat for "id" if need to know atom IDs
   see gather_atoms() to return data ordered by consecutive atom IDs
   see gather_atoms_subset() to return data for only a subset of atoms
   name = desired quantity, e.g. x or charge
   type = 0 for integer values, 1 for double values
   count = # of per-atom values, e.g. 1 for type or charge, 3 for x or f
     use count = 3 with "image" if want single image flag unpacked into xyz
   return atom-based values in 1d data, ordered by count, then by atom
     e.g. x[0][0],x[0][1],x[0][2],x[1][0],x[1][1],x[1][2],x[2][0],...
     data must be pre-allocated by caller to correct length
     correct length = count*Natoms, as queried by get_natoms()
   method:
     Allgather Nlocal atoms from each proc into data
------------------------------------------------------------------------- */

#if defined(LAMMPS_BIGBIG)
void lammps_gather_atoms_concat(void *ptr, char * /*name */,
                                int /*type*/, int /*count*/, void * /*data*/)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  BEGIN_CAPTURE
  {
    lmp->error->all(FLERR,"Library function lammps_gather_atoms_concat() "
                    "is not compatible with -DLAMMPS_BIGBIG");
  }
  END_CAPTURE
}
#else
void lammps_gather_atoms_concat(void *ptr, char *name,
                                int type, int count, void *data)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  BEGIN_CAPTURE
  {
    int i,offset;

    // error if tags are not defined
    // NOTE: test that name = image or ids is not a 64-bit int in code?

    int flag = 0;
    if (lmp->atom->tag_enable == 0) flag = 1;
    if (lmp->atom->natoms > MAXSMALLINT) flag = 1;
    if (flag) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR,"Library error in lammps_gather_atoms");
      return;
    }

    int natoms = static_cast<int> (lmp->atom->natoms);

    void *vptr = lammps_extract_atom(lmp,name);

    if (vptr == NULL) {
      lmp->error->warning(FLERR,"lammps_gather_atoms: unknown property name");
      return;
    }

    // perform MPI_Allgatherv on each proc's chunk of Nlocal atoms

    int nprocs = lmp->comm->nprocs;

    int *recvcounts,*displs;
    lmp->memory->create(recvcounts,nprocs,"lib/gather:recvcounts");
    lmp->memory->create(displs,nprocs,"lib/gather:displs");

    if (type == 0) {
      int *vector = NULL;
      int **array = NULL;
      const int imgunpack = (count == 3) && (strcmp(name,"image") == 0);

      if ((count == 1) || imgunpack) vector = (int *) vptr;
      else array = (int **) vptr;

      int *copy;
      lmp->memory->create(copy,count*natoms,"lib/gather:copy");
      for (i = 0; i < count*natoms; i++) copy[i] = 0;

      int nlocal = lmp->atom->nlocal;

      if (count == 1) {
        MPI_Allgather(&nlocal,1,MPI_INT,recvcounts,1,MPI_INT,lmp->world);
        displs[0] = 0;
        for (i = 1; i < nprocs; i++)
          displs[i] = displs[i-1] + recvcounts[i-1];
        MPI_Allgatherv(vector,nlocal,MPI_INT,data,recvcounts,displs,
                       MPI_INT,lmp->world);

      } else if (imgunpack) {
        int *copy;
        lmp->memory->create(copy,count*nlocal,"lib/gather:copy");
        offset = 0;
        for (i = 0; i < nlocal; i++) {
          const int image = vector[i];
          copy[offset++] = (image & IMGMASK) - IMGMAX;
          copy[offset++] = ((image >> IMGBITS) & IMGMASK) - IMGMAX;
          copy[offset++] = ((image >> IMG2BITS) & IMGMASK) - IMGMAX;
        }
        int n = count*nlocal;
        MPI_Allgather(&n,1,MPI_INT,recvcounts,1,MPI_INT,lmp->world);
        displs[0] = 0;
        for (i = 1; i < nprocs; i++)
          displs[i] = displs[i-1] + recvcounts[i-1];
        MPI_Allgatherv(copy,count*nlocal,MPI_INT,
                       data,recvcounts,displs,MPI_INT,lmp->world);
        lmp->memory->destroy(copy);

      } else {
        int n = count*nlocal;
        MPI_Allgather(&n,1,MPI_INT,recvcounts,1,MPI_INT,lmp->world);
        displs[0] = 0;
        for (i = 1; i < nprocs; i++)
          displs[i] = displs[i-1] + recvcounts[i-1];
        MPI_Allgatherv(&array[0][0],count*nlocal,MPI_INT,
                       data,recvcounts,displs,MPI_INT,lmp->world);
      }

    } else {
      double *vector = NULL;
      double **array = NULL;
      if (count == 1) vector = (double *) vptr;
      else array = (double **) vptr;

      int nlocal = lmp->atom->nlocal;

      if (count == 1) {
        MPI_Allgather(&nlocal,1,MPI_INT,recvcounts,1,MPI_INT,lmp->world);
        displs[0] = 0;
        for (i = 1; i < nprocs; i++)
          displs[i] = displs[i-1] + recvcounts[i-1];
        MPI_Allgatherv(vector,nlocal,MPI_DOUBLE,data,recvcounts,displs,
                       MPI_DOUBLE,lmp->world);

      } else {
        int n = count*nlocal;
        MPI_Allgather(&n,1,MPI_INT,recvcounts,1,MPI_INT,lmp->world);
        displs[0] = 0;
        for (i = 1; i < nprocs; i++)
          displs[i] = displs[i-1] + recvcounts[i-1];
        MPI_Allgatherv(&array[0][0],count*nlocal,MPI_DOUBLE,
                       data,recvcounts,displs,MPI_DOUBLE,lmp->world);
      }
    }

    lmp->memory->destroy(recvcounts);
    lmp->memory->destroy(displs);
  }
  END_CAPTURE
}
#endif

/* ----------------------------------------------------------------------
   gather the named atom-based entity for a subset of atoms
     return it in user-allocated data
   data will be ordered by requested atom IDs
     no requirement for consecutive atom IDs (1 to N)
   see gather_atoms() to return data for all atoms, ordered by consecutive IDs
   see gather_atoms_concat() to return data for all atoms, unordered
   name = desired quantity, e.g. x or charge
   type = 0 for integer values, 1 for double values
   count = # of per-atom values, e.g. 1 for type or charge, 3 for x or f
     use count = 3 with "image" if want single image flag unpacked into xyz
   ndata = # of atoms to return data for (could be all atoms)
   ids = list of Ndata atom IDs to return data for
   return atom-based values in 1d data, ordered by count, then by atom
     e.g. x[0][0],x[0][1],x[0][2],x[1][0],x[1][1],x[1][2],x[2][0],...
     data must be pre-allocated by caller to correct length
     correct length = count*Ndata
   method:
     alloc and zero count*Ndata length vector
     loop over Ndata to fill vector with my values
     Allreduce to sum vector into data across all procs
------------------------------------------------------------------------- */

#if defined(LAMMPS_BIGBIG)
void lammps_gather_atoms_subset(void *ptr, char * /*name */,
                                int /*type*/, int /*count*/,
                                int /*ndata*/, int * /*ids*/, void * /*data*/)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  BEGIN_CAPTURE
  {
    lmp->error->all(FLERR,"Library function lammps_gather_atoms_subset() "
                    "is not compatible with -DLAMMPS_BIGBIG");
  }
  END_CAPTURE
}
#else
void lammps_gather_atoms_subset(void *ptr, char *name,
                                int type, int count,
                                int ndata, int *ids, void *data)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  BEGIN_CAPTURE
  {
    int i,j,m,offset;
    tagint id;

    // error if tags are not defined
    // NOTE: test that name = image or ids is not a 64-bit int in code?

    int flag = 0;
    if (lmp->atom->tag_enable == 0) flag = 1;
    if (lmp->atom->natoms > MAXSMALLINT) flag = 1;
    if (flag) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR,"Library error in lammps_gather_atoms_subset");
      return;
    }

    void *vptr = lammps_extract_atom(lmp,name);
    if (vptr == NULL) {
      lmp->error->warning(FLERR,"lammps_gather_atoms_subset: "
                          "unknown property name");
      return;
    }

    // copy = Ndata length vector of per-atom values
    // use atom ID to insert each atom's values into copy
    // MPI_Allreduce with MPI_SUM to merge into data

    if (type == 0) {
      int *vector = NULL;
      int **array = NULL;
      const int imgunpack = (count == 3) && (strcmp(name,"image") == 0);

      if ((count == 1) || imgunpack) vector = (int *) vptr;
      else array = (int **) vptr;

      int *copy;
      lmp->memory->create(copy,count*ndata,"lib/gather:copy");
      for (i = 0; i < count*ndata; i++) copy[i] = 0;

      int nlocal = lmp->atom->nlocal;

      if (count == 1) {
        for (i = 0; i < ndata; i++) {
          id = ids[i];
          if ((m = lmp->atom->map(id)) >= 0 && m < nlocal)
            copy[i] = vector[m];
        }

      } else if (imgunpack) {
        for (i = 0; i < ndata; i++) {
          id = ids[i];
          if ((m = lmp->atom->map(id)) >= 0 && m < nlocal) {
            offset = count*i;
            const int image = vector[m];
            copy[offset++] = (image & IMGMASK) - IMGMAX;
            copy[offset++] = ((image >> IMGBITS) & IMGMASK) - IMGMAX;
            copy[offset++] = ((image >> IMG2BITS) & IMGMASK) - IMGMAX;
          }
        }

      } else {
        for (i = 0; i < ndata; i++) {
          id = ids[i];
          if ((m = lmp->atom->map(id)) >= 0 && m < nlocal) {
            offset = count*i;
            for (j = 0; j < count; j++)
              copy[offset++] = array[m][j];
          }
        }
      }

      MPI_Allreduce(copy,data,count*ndata,MPI_INT,MPI_SUM,lmp->world);
      lmp->memory->destroy(copy);

    } else {
      double *vector = NULL;
      double **array = NULL;
      if (count == 1) vector = (double *) vptr;
      else array = (double **) vptr;

      double *copy;
      lmp->memory->create(copy,count*ndata,"lib/gather:copy");
      for (i = 0; i < count*ndata; i++) copy[i] = 0.0;

      int nlocal = lmp->atom->nlocal;

      if (count == 1) {
        for (i = 0; i < ndata; i++) {
          id = ids[i];
          if ((m = lmp->atom->map(id)) >= 0 && m < nlocal)
            copy[i] = vector[m];
        }

      } else {
        for (i = 0; i < ndata; i++) {
          id = ids[i];
          if ((m = lmp->atom->map(id)) >= 0 && m < nlocal) {
            offset = count*i;
            for (j = 0; j < count; j++)
              copy[offset++] = array[m][j];
          }
        }
      }

      MPI_Allreduce(copy,data,count*ndata,MPI_DOUBLE,MPI_SUM,lmp->world);
      lmp->memory->destroy(copy);
    }
  }
  END_CAPTURE
}
#endif

/* ----------------------------------------------------------------------
   scatter the named atom-based entity in data to all atoms
   data is ordered by atom ID
     requirement for consecutive atom IDs (1 to N)
   see scatter_atoms_subset() to scatter data for some (or all) atoms, unordered
   name = desired quantity, e.g. x or charge
   type = 0 for integer values, 1 for double values
   count = # of per-atom values, e.g. 1 for type or charge, 3 for x or f
     use count = 3 with "image" for xyz to be packed into single image flag
   data = atom-based values in 1d data, ordered by count, then by atom ID
     e.g. x[0][0],x[0][1],x[0][2],x[1][0],x[1][1],x[1][2],x[2][0],...
     data must be correct length = count*Natoms, as queried by get_natoms()
   method:
     loop over Natoms, if I own atom ID, set its values from data
------------------------------------------------------------------------- */

#if defined(LAMMPS_BIGBIG)
void lammps_scatter_atoms(void *ptr, char * /*name */,
                          int /*type*/, int /*count*/, void * /*data*/)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  BEGIN_CAPTURE
  {
    lmp->error->all(FLERR,"Library function lammps_scatter_atoms() "
                    "is not compatible with -DLAMMPS_BIGBIG");
  }
  END_CAPTURE
}
#else
void lammps_scatter_atoms(void *ptr, char *name,
                          int type, int count, void *data)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  BEGIN_CAPTURE
  {
    int i,j,m,offset;

    // error if tags are not defined or not consecutive or no atom map
    // NOTE: test that name = image or ids is not a 64-bit int in code?

    int flag = 0;
    if (lmp->atom->tag_enable == 0 || lmp->atom->tag_consecutive() == 0)
      flag = 1;
    if (lmp->atom->natoms > MAXSMALLINT) flag = 1;
    if (lmp->atom->map_style == 0) flag = 1;
    if (flag) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR,"Library error in lammps_scatter_atoms");
      return;
    }

    int natoms = static_cast<int> (lmp->atom->natoms);

    void *vptr = lammps_extract_atom(lmp,name);

    if(vptr == NULL) {
        lmp->error->warning(FLERR,
                            "lammps_scatter_atoms: unknown property name");
        return;
    }

    // copy = Natom length vector of per-atom values
    // use atom ID to insert each atom's values into copy
    // MPI_Allreduce with MPI_SUM to merge into data, ordered by atom ID

    if (type == 0) {
      int *vector = NULL;
      int **array = NULL;
      const int imgpack = (count == 3) && (strcmp(name,"image") == 0);

      if ((count == 1) || imgpack) vector = (int *) vptr;
      else array = (int **) vptr;
      int *dptr = (int *) data;

      if (count == 1) {
        for (i = 0; i < natoms; i++)
          if ((m = lmp->atom->map(i+1)) >= 0)
            vector[m] = dptr[i];

      } else if (imgpack) {
        for (i = 0; i < natoms; i++)
          if ((m = lmp->atom->map(i+1)) >= 0) {
            offset = count*i;
            int image = dptr[offset++] + IMGMAX;
            image += (dptr[offset++] + IMGMAX) << IMGBITS;
            image += (dptr[offset++] + IMGMAX) << IMG2BITS;
            vector[m] = image;
          }

      } else {
        for (i = 0; i < natoms; i++)
          if ((m = lmp->atom->map(i+1)) >= 0) {
            offset = count*i;
            for (j = 0; j < count; j++)
              array[m][j] = dptr[offset++];
          }
      }

    } else {
      double *vector = NULL;
      double **array = NULL;
      if (count == 1) vector = (double *) vptr;
      else array = (double **) vptr;
      double *dptr = (double *) data;

      if (count == 1) {
        for (i = 0; i < natoms; i++)
          if ((m = lmp->atom->map(i+1)) >= 0)
            vector[m] = dptr[i];

      } else {
        for (i = 0; i < natoms; i++) {
          if ((m = lmp->atom->map(i+1)) >= 0) {
            offset = count*i;
            for (j = 0; j < count; j++)
              array[m][j] = dptr[offset++];
          }
        }
      }
    }
  }
  END_CAPTURE
}
#endif

/* ----------------------------------------------------------------------
   scatter the named atom-based entity in data to a subset of atoms
   data is ordered by provided atom IDs
     no requirement for consecutive atom IDs (1 to N)
   see scatter_atoms() to scatter data for all atoms, ordered by consecutive IDs
   name = desired quantity, e.g. x or charge
   type = 0 for integer values, 1 for double values
   count = # of per-atom values, e.g. 1 for type or charge, 3 for x or f
     use count = 3 with "image" for xyz to be packed into single image flag
   ndata = # of atoms in ids and data (could be all atoms)
   ids = list of Ndata atom IDs to scatter data to
   data = atom-based values in 1d data, ordered by count, then by atom ID
     e.g. x[0][0],x[0][1],x[0][2],x[1][0],x[1][1],x[1][2],x[2][0],...
     data must be correct length = count*Ndata
   method:
     loop over Ndata, if I own atom ID, set its values from data
------------------------------------------------------------------------- */

#if defined(LAMMPS_BIGBIG)
void lammps_scatter_atoms_subset(void *ptr, char * /*name */,
                                int /*type*/, int /*count*/,
                                int /*ndata*/, int * /*ids*/, void * /*data*/)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  BEGIN_CAPTURE
  {
    lmp->error->all(FLERR,"Library function lammps_scatter_atoms_subset() "
                    "is not compatible with -DLAMMPS_BIGBIG");
  }
  END_CAPTURE
}
#else
void lammps_scatter_atoms_subset(void *ptr, char *name,
                                 int type, int count,
                                 int ndata, int *ids, void *data)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  BEGIN_CAPTURE
  {
    int i,j,m,offset;
    tagint id;

    // error if tags are not defined or no atom map
    // NOTE: test that name = image or ids is not a 64-bit int in code?

    int flag = 0;
    if (lmp->atom->tag_enable == 0) flag = 1;
    if (lmp->atom->natoms > MAXSMALLINT) flag = 1;
    if (lmp->atom->map_style == 0) flag = 1;
    if (flag) {
      if (lmp->comm->me == 0)
        lmp->error->warning(FLERR,"Library error in lammps_scatter_atoms_subset");
      return;
    }

    void *vptr = lammps_extract_atom(lmp,name);

    if(vptr == NULL) {
        lmp->error->warning(FLERR,"lammps_scatter_atoms_subset: "
                                  "unknown property name");
        return;
    }

    // copy = Natom length vector of per-atom values
    // use atom ID to insert each atom's values into copy
    // MPI_Allreduce with MPI_SUM to merge into data, ordered by atom ID

    if (type == 0) {
      int *vector = NULL;
      int **array = NULL;
      const int imgpack = (count == 3) && (strcmp(name,"image") == 0);

      if ((count == 1) || imgpack) vector = (int *) vptr;
      else array = (int **) vptr;
      int *dptr = (int *) data;

      if (count == 1) {
        for (i = 0; i < ndata; i++) {
          id = ids[i];
          if ((m = lmp->atom->map(id)) >= 0)
            vector[m] = dptr[i];
        }

      } else if (imgpack) {
        for (i = 0; i < ndata; i++) {
          id = ids[i];
          if ((m = lmp->atom->map(id)) >= 0) {
            offset = count*i;
            int image = dptr[offset++] + IMGMAX;
            image += (dptr[offset++] + IMGMAX) << IMGBITS;
            image += (dptr[offset++] + IMGMAX) << IMG2BITS;
            vector[m] = image;
          }
        }

      } else {
        for (i = 0; i < ndata; i++) {
          id = ids[i];
          if ((m = lmp->atom->map(id)) >= 0) {
            offset = count*i;
            for (j = 0; j < count; j++)
              array[m][j] = dptr[offset++];
          }
        }
      }

    } else {
      double *vector = NULL;
      double **array = NULL;
      if (count == 1) vector = (double *) vptr;
      else array = (double **) vptr;
      double *dptr = (double *) data;

      if (count == 1) {
        for (i = 0; i < ndata; i++) {
          id = ids[i];
          if ((m = lmp->atom->map(id)) >= 0)
            vector[m] = dptr[i];
        }

      } else {
        for (i = 0; i < ndata; i++) {
          id = ids[i];
          if ((m = lmp->atom->map(id)) >= 0) {
            offset = count*i;
            for (j = 0; j < count; j++)
              array[m][j] = dptr[offset++];
          }
        }
      }
    }
  }
  END_CAPTURE
}
#endif
