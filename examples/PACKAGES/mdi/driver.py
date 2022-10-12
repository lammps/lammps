import sys
import mdi

use_mpi4py = False
try:
    from mpi4py import MPI
    use_mpi4py = True
except:
    pass

# Initialize the MDI Library
mdi.MDI_Init(sys.argv[2])

# Connect to the engine
comm = mdi.MDI_Accept_communicator()

# Determine the name of the engine
mdi.MDI_Send_Command("<NAME", comm)
name = mdi.MDI_Recv(mdi.MDI_NAME_LENGTH, mdi.MDI_CHAR, comm)

print("Engine name: " + str(name))

# Send the "EXIT" command to the engine
mdi.MDI_Send_Command("EXIT", comm)
