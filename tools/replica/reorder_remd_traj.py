#!/usr/bin/env python

"""
LAMMPS Replica Exchange Molecular Dynamics (REMD) trajectories are arranged by
replica, i.e., each trajectory is a continuous replica that records all the
ups and downs in temperature. However, often the requirement is trajectories
that are continuous in temperature, which is achieved by this tool.

Author:
Tanmoy Sanyal, Shell lab, Chemical Engineering, UC Santa Barbara
Email: tanmoy dot 7989 at gmail dot com

Usage
-----
To get detailed information about the arguments, flags, etc use:
python reorder_remd_traj.py -h or
python reorder_remd_traj.py --help

Features of this script
-----------------------
a) reorder LAMMPS REMD trajectories by temperature keeping only desired frames.
Note: this only handles LAMMPS format trajectories (i.e. .lammpstrj format)
Trajectories can be gzipped or bz2-compressed. The trajectories are assumed to
be named as <prefix>.%d.lammpstrj[.gz or .bz2]

b) (optionally) calculate configurational weights for each frame at each
temperature if potential energies are supplied. But this if for the canonical
(NVT) ensemble only.

Dependencies
------------
mpi4py
pymbar (for getting configurational weights)
tqdm (for printing pretty progress bars)
StringIO (or io if in Python 3.x)

"""



import os, numpy as np, argparse, time, pickle
from scipy.special import logsumexp
from mpi4py import MPI

from tqdm import tqdm
import gzip, bz2
try:
    # python-2
    from StringIO import StringIO as IOBuffer
except ImportError:
    # python-3
    from io import BytesIO as IOBuffer



#### INITIALISE MPI ####
# (note that all output on screen will be printed only on the ROOT proc)
ROOT = 0
comm = MPI.COMM_WORLD
me = comm.rank # my proc id
nproc = comm.size


#### HELPER FUNCTIONS ####
def _get_nearest_temp(temps, query_temp):
    """
    Helper function to get the nearest temp in a list
    from a given query_temp

    :param temps: list of temps.

    :param query_temp: query temp

    Returns:
    idx: index of nearest temp in the list

    out_temp: nearest temp from the list
    """

    if isinstance(temps, list): temps = np.array(temps)
    return temps[np.argmin(np.abs(temps-query_temp))]


def readwrite(trajfn, mode):
    """
    Helper function for input/output LAMMPS traj files.
    Trajectories may be plain text, .gz or .bz2 compressed.

    :param trajfn: name of LAMMPS traj

    :param mode: "r" ("w") and "rb" ("wb") depending on read or write

    Returns: file pointer
    """

    if trajfn.endswith(".gz"):
        of = gzip.open(trajfn, mode)
        #return gzip.GzipFile(trajfn, mode)
    elif trajfn.endswith(".bz2"):
        of = bz2.open(trajfn, mode)
        #return bz2.BZ2File(trajfn, mode)
    else:
        of = open(trajfn, mode)
    return of


def get_replica_frames(logfn, temps, nswap, writefreq):
    """
    Get a list of frames from each replica that is
    at a particular temp. Do this for all temps.

    :param logfn: master LAMMPS log file that contains the temp
                  swap history of all replicas

    :param temps: list of all temps used in the REMD simulation.

    :param nswap: swap frequency of the REMD simulation

    :param writefreq: traj dump frequency in LAMMPS

    Returns: master_frametuple_dict:
             dict containing a tuple (replica #, frame #) for each temp.
    """

    n_rep = len(temps)
    swap_history = np.loadtxt(logfn, skiprows = 3)
    master_frametuple_dict = dict( (n, []) for n in range(n_rep) )

    # walk through the replicas
    print("Getting frames from all replicas at temperature:")
    for n in range(n_rep):
        print("%3.2f K" % temps[n])
        rep_inds = [np.where(x[1:] == n)[0][0] for x in swap_history]

        # case-1: when frames are dumped faster than temp. swaps
        if writefreq <= nswap:
            for ii, i in enumerate(rep_inds[:-1]):
                start = int(ii * nswap / writefreq)
                stop = int( (ii+1) * nswap / writefreq)
                [master_frametuple_dict[n].append( (i,x) ) \
                                        for x in range(start, stop)]

        # case-2: when temps. are swapped faster than dumping frames
        else:
            nskip = int(writefreq / nswap)
            [master_frametuple_dict[n].append( (i,ii) ) \
            for ii, i in enumerate(rep_inds[0::nskip])]

    return master_frametuple_dict


def get_byte_index(rep_inds, byteindfns, intrajfns):
    """
    Get byte indices from (un-ordered) trajectories.

    :param rep_inds: indices of replicas to process on this proc

    :param byteindsfns: list of filenames that will contain the byte indices

    :param intrajfns: list of (unordered) input traj filenames
    """
    for n in rep_inds:
        # check if the byte indices for this traj has aleady been computed
        if os.path.isfile(byteindfns[n]): continue

        # extract bytes
        fobj = readwrite(intrajfns[n], "rb")
        byteinds = [ [0,0] ]

        # place file pointer at first line
        nframe = 0
        first_line = fobj.readline()
        cur_pos = fobj.tell()

        # status printed only for replica read on root proc
        # this assumes that each proc takes roughly the same time
        if me == ROOT:
            pb = tqdm(desc = "Reading replicas", leave = True,
                  position = ROOT + 2*me,
                  unit = "B/replica", unit_scale = True,
                  unit_divisor = 1024)

        # start crawling through the bytes
        while True:
            next_line = fobj.readline()
            if len(next_line) == 0: break
            # this will only work with lammpstrj traj format.
            # this condition essentially checks periodic recurrences
            # of the token TIMESTEP. Each time it is found,
            # we have crawled through a frame (snapshot)
            if next_line == first_line:
                nframe += 1
                byteinds.append( [nframe, cur_pos] )
                if me == ROOT: pb.update()
            cur_pos = fobj.tell()
            if me == ROOT: pb.update(0)
        if me == ROOT: pb.close()

        # take care of the EOF
        cur_pos = fobj.tell()
        byteinds.append( [nframe+1, cur_pos] ) # dummy index for the EOF

        # write to file
        np.savetxt(byteindfns[n], np.array(byteinds), fmt = "%d")

        # close the trajfile object
        fobj.close()

        return


def write_reordered_traj(temp_inds, byte_inds, outtemps, temps,
                         frametuple_dict, nprod, writefreq,
                         outtrajfns, infobjs):
    """
    Reorders trajectories by temp. and writes them to disk

    :param temp_inds: list index of temps (in the list of all temps) for which
                      reordered trajs will be produced on this proc.

    :param byte_inds: dict containing the (previously stored) byte indices
                      for each replica file (key = replica number)

    :param outtemps: list of all temps for which to produce reordered trajs.

    :param temps: list of all temps used in the REMD simulation.

    :param outtrajfns: list of filenames for output (ordered) trajs.

    :param frametuple_dict: dict containing a tuple (replica #, frame #)
                            for each temp.

    :param nprod: number of production timesteps.
                  Last (nprod / writefreq) frames
                  from the end will be written to disk.

    :param writefreq: traj dump frequency in LAMMPS

    :param infobjs: list of file pointers to input (unordered) trajs.
    """

    nframes = int(nprod / writefreq)

    for n in temp_inds:
        # open string-buffer and file
        buf = IOBuffer()
        of = readwrite(outtrajfns[n], "wb")

        # get frames
        abs_temp_ind = np.argmin( abs(temps - outtemps[n]) )
        frametuple = frametuple_dict[abs_temp_ind][-nframes:]

        # write frames to buffer
        if me == ROOT:
            pb = tqdm(frametuple,
                  desc = ("Buffering trajectories for writing"),
                  leave = True, position = ROOT + 2*me,
                  unit = 'frame/replica', unit_scale = True)

            iterable = pb
        else:
            iterable = frametuple

        for i, (rep, frame) in enumerate(iterable):
            infobj = infobjs[rep]
            start_ptr = int(byte_inds[rep][frame,1])
            stop_ptr = int(byte_inds[rep][frame+1,1])
            byte_len = stop_ptr - start_ptr
            infobj.seek(start_ptr)
            buf.write(infobj.read(byte_len))
        if me == ROOT: pb.close()

        # write buffer to disk
        if me == ROOT: print("Writing buffer to file")
        of.write(buf.getvalue())
        of.close()
        buf.close()

    for i in infobjs: i.close()

    return


def get_canonical_logw(enefn, frametuple_dict, temps, nprod, writefreq,
                       kB):
    """
    Gets configurational log-weights (logw) for each frame and at each temp.
    from the REMD simulation. ONLY WRITTEN FOR THE CANONICAL (NVT) ensemble.

    This weights can be used to calculate the
    ensemble averaged value of any simulation observable X at a given temp. T :
    <X> (T) = \sum_{k=1, ntemps} \sum_{n=1, nframes} w[idx][k,n] X[k,n]
    where nframes is the number of frames to use from each *reordered* traj

    :param enefn: ascii file (readable by numpy.loadtxt) containing an array
                  u[r,n] of *total* potential energy for the n-th frame for
                  the r-th replica.

    :param frametuple_dict: dict containing a tuple (replica #, frame #)
                            for each temp.

    :param temps: array of temps. used in the REMD simulation

    :param nprod: number of production timesteps. Last (nprod / writefreq)
                  frames from the end will be written to disk.

    :param writefreq: traj dump frequency in LAMMPS

    :param kB : Boltzmann constant to set the energy scale.
                Default is in kcal/mol

    Returns: logw: dict, logw[l][k,n] gives the log weights from the
                   n-th frame of the k-th temp. *ordered* trajectory
                   to reweight to the l-th temp.

    """

    try:
        import pymbar
    except ImportError:
        print("""
              Configurational log-weight calculation requires pymbar.
              Here are some options to install it:
              conda install -c omnia pymbar
              pip install --user pymbar
              sudo pip install pymbar

              To install the dev. version directly from github, use:
              pip install pip install git+https://github.com/choderalab/pymbar.git
              """)

    u_rn = np.loadtxt(enefn)
    ntemps = u_rn.shape[0] # number of temps.
    nframes = int(nprod / writefreq) # number of frames at each temp.

    # reorder the temps
    u_kn = np.zeros([ntemps, nframes], float)
    for k in range(ntemps):
        frame_tuple = frametuple_dict[k][-nframes:]
        for i, (rep, frame) in enumerate(frame_tuple):
            u_kn[k, i] = u_rn[rep, frame]

    # prep input for pymbar
    #1) array of frames at each temp.
    nframes_k = nframes * np.ones(ntemps, np.uint8)

    #2) inverse temps. for chosen energy scale
    beta_k = 1.0 / (kB * temps)

    #3) get reduced energies (*ONLY FOR THE CANONICAL ENSEMBLE*)
    u_kln = np.zeros([ntemps, ntemps, nframes], float)
    for k in range(ntemps):
        u_kln[k] = np.outer(beta_k, u_kn[k])

    # run pymbar and extract the free energies
    print("\nRunning pymbar...")
    mbar = pymbar.mbar.MBAR(u_kln, nframes_k, verbose = True)
    f_k = mbar.f_k # (1 x k array)

    # calculate the log-weights
    print("\nExtracting log-weights...")
    log_nframes = np.log(nframes)
    logw = dict( (k, np.zeros([ntemps, nframes], float)) for k in range(ntemps) )
    # get log-weights to reweight to this temp.
    for k in range(ntemps):
        for n in range(nframes):
            num = -beta_k[k] * u_kn[k,n]
            denom = f_k - beta_k[k] * u_kn[k,n]
            for l in range(ntemps):
                logw[l][k,n] = num - logsumexp(denom) - log_nframes

    return logw



#### MAIN WORKFLOW ####
if __name__ == "__main__":
    # accept user inputs
    parser = argparse.ArgumentParser(description = __doc__,
             formatter_class = argparse.RawDescriptionHelpFormatter)

    parser.add_argument("prefix",
                        help = "Prefix of REMD LAMMPS trajectories.\
                        Supply full path. Trajectories assumed to be named as \
                        <prefix>.%%d.lammpstrj. \
                        Can be in compressed (.gz or .bz2) format. \
                        This is a required argument")

    parser.add_argument("-logfn", "--logfn", default = "log.lammps",
                        help = "LAMMPS log file that contains swap history \
                        of temperatures among replicas. \
                        Default = 'lammps.log'")

    parser.add_argument("-tfn", "--tempfn", default = "temps.txt",
                        help = "ascii file (readable by numpy.loadtxt) with \
                        the temperatures used in the REMD simulation.")

    parser.add_argument("-ns", "--nswap", type = int,
                        help = "Swap frequency used in LAMMPS temper command")

    parser.add_argument("-nw", "--nwrite", type = int, default = 1,
                        help = "Trajectory writing frequency used \
                        in LAMMPS dump command")

    parser.add_argument("-np", "--nprod", type = int, default = 0,
                        help = "Number of timesteps to save in the reordered\
                        trajectories.\
                        This should be in units of the LAMMPS timestep")

    parser.add_argument("-logw", "--logw", action = 'store_true',
                        help = "Supplying this flag \
                        calculates *canonical* (NVT ensemble) log weights")

    parser.add_argument("-e", "--enefn",
                        help = "File that has n_replica x n_frames array\
                        of total potential energies")

    parser.add_argument("-kB", "--boltzmann_const",
                        type = float, default = 0.001987,
                        help = "Boltzmann constant in appropriate units. \
                        Default is kcal/mol")

    parser.add_argument("-ot", "--out_temps", nargs = '+', type = np.float64,
                        help = "Reorder trajectories at these temperatures.\n \
                        Default is all temperatures used in the simulation")

    parser.add_argument("-od", "--outdir", default = ".",
                        help = "All output will be saved to this directory")

    # parse inputs
    args = parser.parse_args()
    traj_prefix = os.path.abspath(args.prefix)
    logfn = os.path.abspath(args.logfn)
    tempfn = os.path.abspath(args.tempfn)

    nswap = args.nswap
    writefreq = args.nwrite
    nprod = args.nprod

    enefn = args.enefn
    if not enefn is None: enefn = os.path.abspath(enefn)
    get_logw = args.logw
    kB = args.boltzmann_const

    out_temps = args.out_temps
    outdir = os.path.abspath(args.outdir)
    if not os.path.isdir(outdir):
        if me == ROOT: os.mkdir(outdir)

    # check that all input files are present (only on the ROOT proc)
    if me == ROOT:
        if not os.path.isfile(tempfn):
            raise IOError("Temperature file %s not found." % tempfn)
        elif not os.path.isfile(logfn):
            raise IOError("LAMMPS log file %s not found." % logfn)
        elif get_logw and not os.path.isfile(enefn):
            raise IOError("Canonical log-weight calculation requested but\
                          energy file %s not found" % enefn)

    # get (unordered) trajectories
    temps = np.loadtxt(tempfn)
    ntemps = len(temps)
    intrajfns = ["%s.%d.lammpstrj" % (traj_prefix, k) for k in range(ntemps)]
    # check if the trajs. (or their zipped versions are present)
    for i in range(ntemps):
        this_intrajfn = intrajfns[i]
        x = this_intrajfn + ".gz"
        if os.path.isfile(this_intrajfn): continue
        elif os.path.isfile(this_intrajfn + ".gz"):
            intrajfns[i] = this_intrajfn + ".gz"
        elif os.path.isfile(this_intrajfn + ".bz2"):
            intrajfns[i] = this_intrajfn + ".bz2"
        else:
            if me == ROOT:
                raise IOError("Trajectory for replica # %d missing" % i)

    # set output filenames
    outprefix = os.path.join(outdir, traj_prefix.split('/')[-1])
    outtrajfns = ["%s.%3.2f.lammpstrj.gz" % \
                 (outprefix, _get_nearest_temp(temps, t)) \
                  for t in out_temps]
    byteindfns = [os.path.join(outdir, ".byteind_%d.gz" % k) \
                  for k in range(ntemps)]
    frametuplefn = outprefix + '.frametuple.pickle'
    if get_logw:
        logwfn = outprefix + ".logw.pickle"


    # get a list of all frames at a particular temp visited by each replica
    # this is fast so run only on ROOT proc.
    master_frametuple_dict = {}
    if me == ROOT:
        master_frametuple_dict = get_replica_frames(logfn = logfn,
                                                    temps = temps,
                                                    nswap = nswap,
                                                    writefreq = writefreq)
        # save to a pickle from the ROOT proc
        with open(frametuplefn, 'wb') as of:
            pickle.dump(master_frametuple_dict, of)

    # broadcast to all procs
    master_frametuple_dict = comm.bcast(master_frametuple_dict, root = ROOT)

    # define a chunk of replicas  to process on each proc
    CHUNKSIZE_1 = int(ntemps/nproc)
    if me < nproc - 1:
        my_rep_inds = range( (me*CHUNKSIZE_1), (me+1)*CHUNKSIZE_1 )
    else:
        my_rep_inds = range( (me*CHUNKSIZE_1), ntemps )

    # get byte indices from replica (un-ordered) trajs. in parallel
    get_byte_index(rep_inds = my_rep_inds,
                   byteindfns = byteindfns,
                   intrajfns = intrajfns)

    # block until all procs have finished
    comm.barrier()

    # open all replica files for reading
    infobjs = [readwrite(i, "rb") for i in intrajfns]

    # open all byteindex files
    byte_inds = dict( (i, np.loadtxt(fn)) for i, fn in enumerate(byteindfns) )

    # define a chunk of output trajs. to process for each proc.
    # # of reordered trajs. to write may be less than the total # of replicas
    # which is usually equal to the requested nproc. If that is indeed the case,
    # retire excess procs
    n_out_temps = len(out_temps)
    CHUNKSIZE_2 = int(n_out_temps / nproc)
    if CHUNKSIZE_2 == 0:
        nproc_active = n_out_temps
        CHUNKSIZE_2 = 1
        if me == ROOT:
            print("\nReleasing %d excess procs" % (nproc - nproc_active))
    else:
        nproc_active = nproc
    if me < nproc_active-1:
        my_temp_inds = range( (me*CHUNKSIZE_2), (me+1)*CHUNKSIZE_1 )
    else:
        my_temp_inds = range( (me*CHUNKSIZE_2), n_out_temps)

    # retire the excess procs
    # dont' forget to close any open file objects
    if me >= nproc_active:
        for fobj in infobjs: fobj.close()
        exit()

    # write reordered trajectories to disk from active procs in parallel
    write_reordered_traj(temp_inds = my_temp_inds,
                         byte_inds = byte_inds,
                         outtemps = out_temps, temps = temps,
                         frametuple_dict = master_frametuple_dict,
                         nprod = nprod, writefreq = writefreq,
                         outtrajfns = outtrajfns,
                         infobjs = infobjs)

    # calculate canonical log-weights if requested
    # usually this is very fast so retire all but the ROOT proc
    if not get_logw: exit()
    if not me == ROOT: exit()

    logw = get_canonical_logw(enefn = enefn, temps = temps,
                              frametuple_dict = master_frametuple_dict,
                              nprod = nprod, writefreq = writefreq,
                              kB = kB)


    # save the logweights to a pickle
    with open(logwfn, 'wb') as of:
        pickle.dump(logw, of)


