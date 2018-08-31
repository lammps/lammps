# sample launch script

# message on 1 proc each

mpirun -np 1 lmp_mpi -log log.message.g++.1 < in.message

mpirun -np 1 lmp_mpi -v mode file -log log.message.client.file.g++.1 < in.message.client & 
mpirun -np 1 lmp_mpi -v mode file -log log.message.server.file.g++.1 < in.message.server

mpirun -np 1 lmp_mpi -v mode zmq -log log.message.client.zmq.g++.1 < in.message.client & 
mpirun -np 1 lmp_mpi -v mode zmq -log log.message.server.zmq.g++.1 < in.message.server

mpirun -np 1 lmp_mpi -v mode mpitwo -log log.message.client.mpitwo.g++.1 < in.message.client & 
mpirun -np 1 lmp_mpi -v mode mpitwo -log log.message.server.mpitwo.g++.1 < in.message.server

mpirun -np 1 lmp_mpi -m 0 -in in.message.client -v mode mpione -log log.message.client.mpione.g++.1 : -np 1 lmp_mpi -m 1 -in in.message.server -v mode mpione -log log.message.server.mpione.g++.1 

# message on 2/4 procs each

mpirun -np 4 lmp_mpi -log log.message.g++.4 < in.message

mpirun -np 2 lmp_mpi -v mode file -log log.message.client.file.g++.2 < in.message.client & 
mpirun -np 4 lmp_mpi -v mode file -log log.message.server.file.g++.4 < in.message.server

mpirun -np 2 lmp_mpi -v mode zmq -log log.message.client.zmq.g++.2 < in.message.client & 
mpirun -np 4 lmp_mpi -v mode zmq -log log.message.server.zmq.g++.4 < in.message.server

mpirun -np 2 lmp_mpi -v mode mpitwo -log log.message.client.mpitwo.g++.2 < in.message.client & 
mpirun -np 4 lmp_mpi -v mode mpitwo -log log.message.server.mpitwo.g++.4 < in.message.server

mpirun -np 2 lmp_mpi -m 0 -in in.message.client -v mode mpione -log log.message.client.mpione.g++.2 : -np 4 lmp_mpi -m 1 -in in.message.server -v mode mpione -log log.message.server.mpione.g++.4

# message.tilt on 1 proc each

mpirun -np 1 lmp_mpi -log log.message.tilt.g++.1 < in.message.tilt

mpirun -np 1 lmp_mpi -v mode zmq -log log.message.tilt.client.zmq.g++.1 < in.message.tilt.client & 
mpirun -np 1 lmp_mpi -v mode zmq -log log.message.tilt.server.zmq.g++.1 < in.message.tilt.server

mpirun -np 1 lmp_mpi -v mode mpitwo -log log.message.tilt.client.mpitwo.g++.1 < in.message.tilt.client & 
mpirun -np 1 lmp_mpi -v mode mpitwo -log log.message.tilt.server.mpitwo.g++.1 < in.message.tilt.server

mpirun -np 1 lmp_mpi -m 0 -in in.message.tilt.client -v mode mpione -log log.message.tilt.client.mpione.g++.1 : -np 1 lmp_mpi -m 1 -in in.message.tilt.server -v mode mpione -log log.message.tilt.server.mpione.g++.1 

# message.tilt on 2/4 procs each

mpirun -np 1 lmp_mpi -log log.message.tilt.g++.4 < in.message.tilt

mpirun -np 2 lmp_mpi -v mode zmq -log log.message.tilt.client.zmq.g++.2 < in.message.tilt.client & 
mpirun -np 4 lmp_mpi -v mode zmq -log log.message.tilt.server.zmq.g++.4 < in.message.tilt.server

mpirun -np 2 lmp_mpi -v mode mpitwo -log log.message.tilt.client.mpitwo.g++.2 < in.message.tilt.client & 
mpirun -np 4 lmp_mpi -v mode mpitwo -log log.message.tilt.server.mpitwo.g++.4 < in.message.tilt.server

mpirun -np 2 lmp_mpi -m 0 -in in.message.tilt.client -v mode mpione -log log.message.tilt.client.mpione.g++.2 : -np 4 lmp_mpi -m 1 -in in.message.tilt.server -v mode mpione -log log.message.tilt.server.mpione.g++.4
