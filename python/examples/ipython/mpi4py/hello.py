from mpi4py import MPI

comm=MPI.COMM_WORLD
print("Hello from rank %d of %d" % (comm.rank, comm.size))
