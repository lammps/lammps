def compile_on_env(image_name, compiler) {
    def envImage = docker.image(image_name)

    docker.withRegistry('https://registry.hub.docker.com', 'docker-registry-login') {
        // ensure image is current
        envImage.pull()

        // use workaround (see https://issues.jenkins-ci.org/browse/JENKINS-34276)
        docker.image(envImage.imageName()).inside {
            switch(compiler) {
                case 'gcc':
                    env.CC = 'gcc'
                    env.CXX = 'g++'
                    env.OMPI_CC = 'gcc'
                    env.OMPI_CXX = 'g++'
                    break

                case 'clang':
                    env.CC = 'clang'
                    env.CXX = 'clang++'
                    env.OMPI_CC = 'clang'
                    env.OMPI_CXX = 'clang++'
                    break
            }

            sh 'ccache -C'
            sh 'ccache -M 5G'

            // clean up project directory
            sh '''
            date
            make -C src clean-all
            make -C src yes-all
            make -C src no-lib
            date
            '''

            // build libraries
            sh '''
            date
            make -C lib/colvars clean-all
            make -j 8 -C lib/colvars -f Makefile.g++ CXX="${COMP}"
            make -j 8 -C lib/poems -f Makefile.g++ CXX="${COMP}"
            make -j 8 -C lib/voronoi -f Makefile.g++ CXX="${COMP}"
            make -j 8 -C lib/awpmd -f Makefile.mpicc CC="${COMP}"
            make -j 8 -C lib/meam -f Makefile.gfortran CC=gcc F90=gfortran
            make -j 8 -C lib/h5md
            date
            '''

            // enable modules
            sh '''
            date
            make -C src yes-user-smd yes-user-molfile yes-compress yes-python
            make -C src yes-poems yes-voronoi yes-user-colvars yes-user-awpmd yes-meam
            make -C src yes-user-h5md
            date
            '''

            // add additonal modules if MPI is used
            if(env.MACH == "mpi") {
                sh 'make -C src yes-mpiio yes-user-lb'
            }

            switch(env.MACH) {
                case 'serial':
                case 'mpi':
                    sh '''
                    date
                    mpicxx -v
                    make -j 8 -C src ${MACH} MPICMD="${MPICMD}" CC="${COMP}" LINK="${COMP}" LMP_INC="${LMP_INC}" JPG_LIB="${JPG_LIB}" TAG="${TAG}-$CC" LMPFLAGS="${LMPFLAGS}"
                    date
                    '''
                    sh '''
                    make -C src test-${MACH} MPICMD="${MPICMD}" TAG="${TAG}-$CC" LMPFLAGS="${LMPFLAGS}"
                    '''
                    break

                case 'shlib':
                    sh '''
                    date
                    make -j 8 -C src mode=shlib serial MACH=serial MPICMD="${MPICMD}" CC="${COMP}" LINK="${COMP}" LMP_INC="${LMP_INC}" JPG_LIB="${JPG_LIB}" TAG="${TAG}-$CC" LMPFLAGS="${LMPFLAGS}"
                    date
                    '''
                    break
            }

            sh 'ccache -s'
        }
    }
}

node {
    stage 'Checkout'
    checkout scm

    step([$class: 'GitHubCommitStatusSetter', contextSource: [$class: 'ManuallyEnteredCommitContextSource', context: 'continuous-integration/jenkins'], statusResultSource: [$class: 'ConditionalStatusResultSource', results: [[$class: 'AnyBuildResult', message: '', state: 'PENDING']]]])

    env.CCACHE_DIR= pwd() + '/.ccache'

    stage 'Serial Binary'

    env.COMP     = 'g++'
    env.MACH     = 'serial'
    env.LMPFLAGS = '-sf off'
    env.LMP_INC  = '-I../../src/STUBS -DFFT_KISSFFT -DLAMMPS_GZIP -DLAMMPS_PNG -DLAMMPS_JPEG'
    env.JPG_LIB  = '-L../../src/STUBS/ -lmpi_stubs -ljpeg -lpng -lz'

    compile_on_env('rbberger/lammps-testing:ubuntu_latest', 'gcc')

    stage 'Shared Library'

    env.COMP     = 'g++'
    env.MACH     = 'shlib'
    env.LMPFLAGS = '-sf off'
    env.LMP_INC  = '-I../../src/STUBS -DFFT_KISSFFT -DLAMMPS_GZIP -DLAMMPS_PNG -DLAMMPS_JPEG'
    env.JPG_LIB  = '-L../../src/STUBS/ -lmpi_stubs -ljpeg -lpng -lz'

    compile_on_env('rbberger/lammps-testing:ubuntu_latest', 'gcc')

    stage 'OpenMPI binary'

    env.COMP     = 'mpicxx'
    env.MACH     = 'mpi'
    env.MPICMD   = 'mpirun -np 4'
    env.LMPFLAGS = '-sf off'
    env.LMP_INC  = '-DFFT_KISSFFT -DLAMMPS_GZIP -DLAMMPS_PNG -DLAMMPS_JPEG -DLAMMPS_SMALLSMALL'
    env.JPG_LIB  = '-ljpeg -lpng -lz'

    compile_on_env('rbberger/lammps-testing:ubuntu_latest', 'gcc')

    step([$class: 'WarningsPublisher', canComputeNew: false, consoleParsers: [[parserName: 'GNU Make + GNU C Compiler (gcc)']], defaultEncoding: '', excludePattern: '', healthy: '', includePattern: '', messagesPattern: '', unHealthy: ''])

    step([$class: 'GitHubCommitStatusSetter', contextSource: [$class: 'ManuallyEnteredCommitContextSource', context: 'continuous-integration/jenkins']])
}
