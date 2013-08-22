///////////////////////////////////////////////////////////////////////////////
//
//  Copyright 2010, Alexander Stukowski
//  All rights reserved. See README.txt for more information.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef __DXA_DISLOCATION_EXTRACTOR_H
#define __DXA_DISLOCATION_EXTRACTOR_H

// Include the LAMMPS headers that we need.
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "domain.h"
#include "update.h"
#include "neigh_list.h"
using namespace LAMMPS_NS;

// This is must be defined when we build the LAMMPS module (and not the DXA standalone tool).
#define DXA_LAMMPS_MODULE

// Enables additional sanity checks.
#define DEBUG_DISLOCATIONS

// In LAMMPS, we always use double precision floating-point numbers.
#define USE_DOUBLE_PRECISION_FP

#include "../src/dxa/DXATracing.h"
#include "../src/util/FILEStream.h"

/**
 *
 */
class DXADislocationExtractor : private DXATracing, Pointers
{
public:

  /// Constructor.
  ///
  /// Expects a pointer to the global LAMMPS object.
  ///
  /// Exact mode flag:
  ///
  /// Enables serial processing, i.e. the master processor performs the complete analysis.
  /// This allows for proper handling of periodic boundary conditions.
  DXADislocationExtractor(LAMMPS* lmp, bool _exactMode = false) : DXATracing(msgLoggerStream, verboseLoggerStream), Pointers(lmp)
  {
    exactMode = _exactMode;

    // Connect internal logging streams to global LAMMPS streams.
    if(comm->me == 0) {
      msgLoggerStream.rdbuf()->setFiles(screen, logfile);
      verboseLoggerStream.rdbuf()->setFiles(screen, logfile);
    }
    this->processor = comm->me;
  }

  /// Extracts all dislocations from the current simulation state.
  ///
  /// First copies atoms from LAMMPS array to internal memory.
  /// Also adopts the nearest neighbor lists from LAMMPS.
  ///
  /// Parameters:
  ///
  ///      neighborList   A full neighbor list provided by LAMMPS, which contains neighbors of ghost atoms.
  ///      cnaCutoff      Cutoff radius to be used for the common neighbor analysis.
  ///
  void extractDislocations(NeighList* neighborList, double cnaCutoff, int lineSmoothingLevel = DEFAULT_LINE_SMOOTHING_LEVEL, int lineCoarseningLevel = DEFAULT_LINE_COARSENING_LEVEL)
  {
    // Check parameters.
    if(cnaCutoff <= 0.0) lammps->error->all(FLERR,"Common neighbor analysis cutoff radius must be positive.");

    // Release any memory allocated during last call.
    cleanup();

    // Setup simulation cell.
    if(exactMode) {
      pbc[0] = domain->periodicity[0] != 0;
      pbc[1] = domain->periodicity[1] != 0;
      pbc[2] = domain->periodicity[2] != 0;
    }
    else {
      pbc[0] = pbc[1] = pbc[2] = false;
    }
    simulationCell(0,0) = domain->h[0];
    simulationCell(1,1) = domain->h[1];
    simulationCell(2,2) = domain->h[2];
    simulationCell(1,2) = domain->h[3];
    simulationCell(0,2) = domain->h[4];
    simulationCell(0,1) = domain->h[5];
    simulationCell(1,0) = simulationCell(2,0) = simulationCell(2,1) = 0;
    simulationCellOrigin.X = domain->boxlo[0];
    simulationCellOrigin.Y = domain->boxlo[1];
    simulationCellOrigin.Z = domain->boxlo[2];
    setCNACutoff(cnaCutoff);
    setupSimulationCell(cnaCutoff);
    timestep = update->ntimestep;

    // Transfer LAMMPS atoms to internal array.
    if(exactMode)
      transferAllLAMMPSAtoms();
    else
      transferLocalLAMMPSAtoms(neighborList);

    if(exactMode == false || processor == 0) {

      // Perform common neighbor analysis to identify crystalline atoms.
      performCNA();

      // Order the neighbor lists of crystalline atoms.
      orderCrystallineAtoms();

      // Cluster crystalline atoms.
      clusterAtoms();

#ifdef PROCESSED_ATOMS
      // Enable this code block to get the processed atoms.
      {
        ostringstream ss;
        ss << "proc_" << comm->me << ".atoms.dump";
        ofstream stream(ss.str().c_str());
        writeAtomsDumpFile(stream);
      }
#endif

      // Create the nodes of the interface mesh.
      createInterfaceMeshNodes();

      // Create the facets of the interface mesh.
      createInterfaceMeshFacets();

#ifdef DEBUG_DISLOCATIONS
      // Check the generated mesh.
      validateInterfaceMesh();
#endif

      // Generate and advance Burgers circuits on the interface mesh.
      traceDislocationSegments();

      // Smooth dislocation lines.
      smoothDislocationSegments(lineSmoothingLevel, lineCoarseningLevel);
    }

    if(exactMode) {
      // Wrap segments at simulation cell boundaries.
      wrapDislocationSegments();
      // Distribute extracted dislocation segments to back all processors.
      broadcastDislocations();
    }

    // Clip dislocation lines at processor domain.
    Point3 subOrigin;
    Matrix3 subCell;
        if(domain->triclinic) {
      subOrigin = simulationCellOrigin + simulationCell * Vector3(domain->sublo_lamda);
      subCell.setColumn(0, simulationCell.column(0) * (domain->subhi_lamda[0] - domain->sublo_lamda[0]));
      subCell.setColumn(1, simulationCell.column(1) * (domain->subhi_lamda[1] - domain->sublo_lamda[1]));
      subCell.setColumn(2, simulationCell.column(2) * (domain->subhi_lamda[2] - domain->sublo_lamda[2]));
        } else {
      subOrigin = Point3(domain->sublo);
      subCell.setColumn(0, Vector3(domain->subhi[0] - domain->sublo[0], 0, 0));
      subCell.setColumn(1, Vector3(0, domain->subhi[1] - domain->sublo[1], 0));
      subCell.setColumn(2, Vector3(0, 0, domain->subhi[2] - domain->sublo[2]));
        }
    clipDislocationLines(subOrigin, subCell);

#ifdef OUTPUT_SEGMENTS
    // Enable this code block to see the extracted dislocation segments on each processor.
    {
      ostringstream ss;
      ss << "proc_" << comm->me << ".dislocations.vtk";
      ofstream stream(ss.str().c_str());
      writeDislocationsVTKFile(stream);
    }
#endif
  }

  /// Returns the list of extracted dislocation segments.
  const vector<DislocationSegment*>& getSegments() const { return this->segments; }

protected:

  /// Copies local atoms and neighbor lists from LAMMPS to internal array.
  void transferLocalLAMMPSAtoms(NeighList* neighborList)
  {
    // Allocate internal atoms array.
    int numTotalAtoms = atom->nlocal + atom->nghost;
    inputAtoms.resize(numTotalAtoms);
    numLocalInputAtoms = numTotalAtoms;

    // Transfer local and ghost atoms to internal array.
    for(int i = 0; i < numTotalAtoms; i++) {
      InputAtom& a = inputAtoms[i];
      a.tag = i + 1;
      a.pos.X = atom->x[i][0];
      a.pos.Y = atom->x[i][1];
      a.pos.Z = atom->x[i][2];
      a.flags = 0;
      a.cluster = NULL;
      a.numNeighbors = 0;
      a.setFlag(ATOM_IS_LOCAL_ATOM);
    }

    // Adopt nearest neighbor lists from LAMMPS.
    if(neighborList->ghostflag) {
      DISLOCATIONS_ASSERT(neighborList->inum + neighborList->gnum == numTotalAtoms);
      int tt = 0;
      for(int ii = 0; ii < numTotalAtoms; ii++) {
        int i = neighborList->ilist[ii];
        DISLOCATIONS_ASSERT(i < numTotalAtoms);
        InputAtom& ai = inputAtoms[i];
        int* jlist = neighborList->firstneigh[i];
        int jnum = neighborList->numneigh[i];
        if(ii >= neighborList->inum)
          tt += jnum;
        for(int jj = 0; jj < jnum; jj++) {
          int j = jlist[jj];
          j &= NEIGHMASK;
          DISLOCATIONS_ASSERT(j >= 0 && j < numTotalAtoms);
          if(j < i) continue;
          InputAtom& aj = inputAtoms[j];
          if(DistanceSquared(aj.pos, ai.pos) >= cnaCutoff * cnaCutoff) continue;
          if(ai.numNeighbors == MAX_ATOM_NEIGHBORS)
            raiseError("Maximum number of nearest neighbors exceeded. Atom %i has more than %i nearest neighbors (built-in maximum number).", ai.tag, MAX_ATOM_NEIGHBORS);
          ai.addNeighbor(&aj);
          if(aj.numNeighbors == MAX_ATOM_NEIGHBORS)
            raiseError("Maximum number of nearest neighbors exceeded. Atom %i has more than %i nearest neighbors (built-in maximum number).", aj.tag, MAX_ATOM_NEIGHBORS);
          aj.addNeighbor(&ai);
        }
      }
    }
    else {
      for(int ii = 0; ii < neighborList->inum; ii++) {
        int i = neighborList->ilist[ii];
        InputAtom& ai = inputAtoms[i];
        int* jlist = neighborList->firstneigh[i];
        int jnum = neighborList->numneigh[i];
        for(int jj = 0; jj < jnum; jj++) {
          int j = jlist[jj];
          j &= NEIGHMASK;
          DISLOCATIONS_ASSERT(j >= 0 && j < numTotalAtoms);
          if(j < i) continue;
          InputAtom& aj = inputAtoms[j];
          if(DistanceSquared(aj.pos, ai.pos) < cnaCutoff * cnaCutoff) {
            if(ai.numNeighbors == MAX_ATOM_NEIGHBORS)
              raiseError("Maximum number of nearest neighbors exceeded. Atom %i has more than %i nearest neighbors (built-in maximum number).", ai.tag, MAX_ATOM_NEIGHBORS);
            ai.addNeighbor(&aj);
            if(aj.numNeighbors == MAX_ATOM_NEIGHBORS)
              raiseError("Maximum number of nearest neighbors exceeded. Atom %i has more than %i nearest neighbors (built-in maximum number).", aj.tag, MAX_ATOM_NEIGHBORS);
            aj.addNeighbor(&ai);
          }
          if(j >= atom->nlocal) {
            for(int kk = 0; kk < jnum; kk++) {
              int k = jlist[kk];
              k &= NEIGHMASK;
              DISLOCATIONS_ASSERT(k < (int)inputAtoms.size());
              if(k == j) continue;
              if(k < atom->nlocal) continue;
              InputAtom& ak = inputAtoms[k];
              if(DistanceSquared(aj.pos, ak.pos) >= cnaCutoff * cnaCutoff) continue;
              if(aj.hasNeighbor(&ak)) continue;
              if(aj.numNeighbors == MAX_ATOM_NEIGHBORS)
                raiseError("Maximum number of nearest neighbors exceeded. Atom %i has more than %i nearest neighbors (built-in maximum number).", aj.tag, MAX_ATOM_NEIGHBORS);
              aj.addNeighbor(&ak);
            }
          }
        }
      }
    }
  }

  /// Copies atoms from LAMMPS to internal array and gathers all atoms on the master processor.
  void transferAllLAMMPSAtoms()
  {
    // Determine maximum number of local atoms a single processor has.
    int nlocalatoms_max;
    MPI_Reduce(&atom->nlocal, &nlocalatoms_max, 1, MPI_INT, MPI_MAX, 0, world);

      if(processor == 0) {
      // Allocate internal atoms array.
      inputAtoms.resize(atom->natoms);
      numLocalInputAtoms = atom->natoms;
      vector<InputAtom>::iterator currentAtom = inputAtoms.begin();

      // Copy local atoms of master processor into internal array.
      for(int i = 0; i < atom->nlocal; i++, ++currentAtom) {
        currentAtom->tag = currentAtom - inputAtoms.begin() + 1;
        currentAtom->pos.X = atom->x[i][0];
        currentAtom->pos.Y = atom->x[i][1];
        currentAtom->pos.Z = atom->x[i][2];
        currentAtom->flags = 0;
        currentAtom->cluster = NULL;
        currentAtom->numNeighbors = 0;
        currentAtom->setFlag(ATOM_IS_LOCAL_ATOM);
      }

      if(comm->nprocs > 1) {
        // Allocate receive buffer.
        vector<double> buffer(nlocalatoms_max * 3);

        // Receive atoms from other processors.
        for(int iproc = 1; iproc < comm->nprocs; iproc++) {
          MPI_Status status;
          MPI_Recv(buffer.empty() ? NULL : &buffer.front(), nlocalatoms_max * 3, MPI_DOUBLE, iproc, 0, world, &status);
          int ndoubles;
          MPI_Get_count(&status, MPI_DOUBLE, &ndoubles);
          int nReceived = ndoubles / 3;

          vector<double>::const_iterator data = buffer.begin();
          for(int i = 0; i < nReceived; i++, ++currentAtom) {
            currentAtom->tag = currentAtom - inputAtoms.begin() + 1;
            currentAtom->pos.X = *data++;
            currentAtom->pos.Y = *data++;
            currentAtom->pos.Z = *data++;
            currentAtom->flags = 0;
            currentAtom->cluster = NULL;
            currentAtom->numNeighbors = 0;
            currentAtom->setFlag(ATOM_IS_LOCAL_ATOM);
          }
        }
      }
        DISLOCATIONS_ASSERT(currentAtom == inputAtoms.end());
    }
    else {
      // Allocate and fill send buffer.
      vector<double> buffer(atom->nlocal * 3);
      vector<double>::iterator data = buffer.begin();
      for(int i = 0; i < atom->nlocal; i++) {
        *data++ = atom->x[i][0];
        *data++ = atom->x[i][1];
        *data++ = atom->x[i][2];
      }
      // Send local atom coordinates to master proc.
      MPI_Send(buffer.empty() ? NULL : &buffer.front(), buffer.size(), MPI_DOUBLE, 0, 0, world);
    }

      // Make sure all input atoms are wrapped at periodic boundary conditions.
      wrapInputAtoms(NULL_VECTOR);

      // Build nearest neighbor lists.
      buildNearestNeighborLists();
  }

  /// This structure is used to communicate dislocation segments between processors.
  struct DislocationSegmentComm {
    LatticeVector burgersVector;
    Vector3 burgersVectorWorld;
    int numPoints;
  };

  /// After the master processor has extracted all dislocation segments, this function broadcasts the dislocations
  /// back to all other processor.
  void broadcastDislocations()
  {
    if(comm->nprocs == 1) return;  // Nothing to do in serial mode.

    // Broadcast number of segments.
    int numSegments = segments.size();
    MPI_Bcast(&numSegments, 1, MPI_INT, 0, world);

    // Allocate send/receive buffer for dislocation segments.
    vector<DislocationSegmentComm> segmentBuffer(numSegments);

    // The total number of line points (sum of all segments).
    int numLinePoints = 0;

    if(processor == 0) {
      // Fill segment send buffer.
      vector<DislocationSegmentComm>::iterator sendItem = segmentBuffer.begin();
      for(vector<DislocationSegment*>::const_iterator segment = segments.begin(); segment != segments.end(); ++segment, ++sendItem) {
        sendItem->burgersVector = (*segment)->burgersVector;
        sendItem->burgersVectorWorld = (*segment)->burgersVectorWorld;
        sendItem->numPoints = (*segment)->line.size();
        numLinePoints += (*segment)->line.size();
      }
    }
    // Broadcast segments.
    MPI_Bcast(segmentBuffer.empty() ? NULL : &segmentBuffer.front(), segmentBuffer.size() * sizeof(segmentBuffer[0]), MPI_CHAR, 0, world);

    if(processor != 0) {
      // Extract segments from receive buffer.
      segments.reserve(segmentBuffer.size());
      for(vector<DislocationSegmentComm>::const_iterator receiveItem = segmentBuffer.begin(); receiveItem != segmentBuffer.end(); ++receiveItem) {
        DislocationSegment* newSegment = segmentPool.construct(receiveItem->burgersVector, receiveItem->burgersVectorWorld);
        newSegment->index = segments.size() + 1;
        segments.push_back(newSegment);
        numLinePoints += receiveItem->numPoints;
      }
    }

    // Allocate send/receive buffer for dislocation points.
    vector<Point3> pointBuffer(numLinePoints);

    if(processor == 0) {
      // Fill point send buffer.
      vector<Point3>::iterator sendItem = pointBuffer.begin();
      for(vector<DislocationSegment*>::const_iterator segment = segments.begin(); segment != segments.end(); ++segment) {
        for(deque<Point3>::const_iterator p = (*segment)->line.begin(); p != (*segment)->line.end(); ++p, ++sendItem)
          *sendItem = *p;
      }
      DISLOCATIONS_ASSERT(sendItem == pointBuffer.end());
    }
    // Broadcast segments.
    MPI_Bcast(pointBuffer.empty() ? NULL : &pointBuffer.front(), pointBuffer.size() * sizeof(pointBuffer[0]), MPI_CHAR, 0, world);

    if(processor != 0) {
      // Extract points from receive buffer.
      vector<DislocationSegment*>::const_iterator segment = segments.begin();
      vector<Point3>::const_iterator p = pointBuffer.begin();
      for(vector<DislocationSegmentComm>::const_iterator receiveItem = segmentBuffer.begin(); receiveItem != segmentBuffer.end(); ++receiveItem, ++segment) {
        (*segment)->line.assign(p, p + receiveItem->numPoints);
        p += receiveItem->numPoints;
      }
      DISLOCATIONS_ASSERT(p == pointBuffer.end());
    }
  }

protected:

  /// Pointer to the main LAMMPS object.
  LAMMPS* lammps;

  /// Enables serial processing, which provides exact handling of periodic boundary conditions.
  bool exactMode;

  /// The output stream to which log messages are sent.
  stdio_osyncstream msgLoggerStream;

  /// The output stream to which verbose log messages are sent.
  stdio_osyncstream verboseLoggerStream;;
};

#endif // __DXA_DISLOCATION_EXTRACTOR_H
