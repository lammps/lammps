from fedora:latest

RUN yum update -y
RUN yum install -y make gcc cmake git openmpi openmpi-devel gcc-gfortran libpng-devel libjpeg-devel openmpi-devel python3-devel fftw3-devel zlib-devel gsl-devel voro++-devel gcc-c++ ocl-icd-devel opencl-headers ccache

RUN useradd -m lammps
USER lammps
ENV PATH=/usr/lib64/ccache:/usr/lib64/openmpi/bin${PATH:+:}${PATH}
ENV LD_LIBRARY_PATH=/usr/lib64/openmpi/lib${LD_LIBRARY_PATH:+:}${LD_LIBRARY_PATH}
WORKDIR /home/lammps
RUN mkdir -p .ccache

COPY lammps/ /home/lammps/lammps
RUN rm -rf /home/lammps/.ccache
COPY ccache/ /home/lammps/.ccache
USER root
RUN chown -R lammps:lammps /home/lammps/lammps /home/lammps/.ccache
USER lammps

WORKDIR /home/lammps/lammps
RUN mkdir build
WORKDIR build
RUN ccache -z
RUN cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=/usr \
  -DENABLE_ALL=ON \
  -DLAMMPS_XDR=ON  \
  -DENABLE_PYTHON=ON \
  -DENABLE_VORONOI=ON \
  -DENABLE_TESTING=ON \
  -DENABLE_GPU=ON -DGPU_API=OpenCL \
  -DFFT=FFTW3 \
  -DENABLE_MPI=ON \
  -DENABLE_MPIIO=ON \
  ../cmake
RUN make VERBOSE=1 -j2
RUN ccache -s
RUN make install DESTDIR=${PWD}/install && rm -rf ${PWD}/install/usr && rmdir ${PWD}/install
USER root
RUN make install
USER lammps
