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

            //sh 'ccache -C'
            sh 'ccache -M 5G'

            // clean up project directory
            sh '''
            make -C src clean-all
            make -C src yes-all
            make -C src no-lib
            '''

            // build libraries
            sh '''
            make -C lib/colvars clean-all
            make -j 8 -C lib/colvars -f Makefile.g++ CXX="${COMP}"
            make -j 8 -C lib/poems -f Makefile.g++ CXX="${COMP}"
            make -j 8 -C lib/voronoi -f Makefile.g++ CXX="${COMP}"
            make -j 8 -C lib/awpmd -f Makefile.mpicc CC="${COMP}"
            make -j 8 -C lib/meam -f Makefile.gfortran CC=gcc F90=gfortran
            make -j 8 -C lib/h5md
            '''

            // enable modules
            sh '''
            make -C src yes-user-smd yes-user-molfile yes-compress yes-python
            make -C src yes-poems yes-voronoi yes-user-colvars yes-user-awpmd yes-meam
            make -C src yes-user-h5md
            '''

            // add additonal modules if MPI is used
            if(env.MACH == "mpi") {
                sh 'make -C src yes-mpiio yes-user-lb'
            }

            switch(env.MACH) {
                case 'serial':
                case 'mpi':
                    sh '''
                    make -j 8 -C src ${MACH} MPICMD="${MPICMD}" CC="${COMP}" LINK="${COMP}" LMP_INC="${LMP_INC}" JPG_LIB="${JPG_LIB}" TAG="${TAG}-$CC" LMPFLAGS="${LMPFLAGS}"
                    '''
                    //sh '''
                    //make -C src test-${MACH} MPICMD="${MPICMD}" TAG="${TAG}-$CC" LMPFLAGS="${LMPFLAGS}"
                    //'''
                    break

                case 'shlib':
                    sh '''
                    make -j 8 -C src mode=shlib serial MACH=serial MPICMD="${MPICMD}" CC="${COMP}" LINK="${COMP}" LMP_INC="${LMP_INC}" JPG_LIB="${JPG_LIB}" TAG="${TAG}-$CC" LMPFLAGS="${LMPFLAGS}"
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

    stash includes: 'src/**', name: 'source'
    stash includes: 'lib/**', name: 'libraries'

    env.CCACHE_DIR= pwd() + '/.ccache'

    parallel (
        "Serial Binary" : {
            node {
                unstash 'source'
                unstash 'libraries'

                step([$class: 'GitHubCommitStatusSetter', contextSource: [$class: 'ManuallyEnteredCommitContextSource', context: 'continuous-integration/jenkins/serial'], statusResultSource: [$class: 'ConditionalStatusResultSource', results: [[$class: 'AnyBuildResult', message: 'running serial compilation', state: 'PENDING']]]])

                env.COMP     = 'g++'
                env.MACH     = 'serial'
                env.LMPFLAGS = '-sf off'
                env.LMP_INC  = '-I../../src/STUBS -DFFT_KISSFFT -DLAMMPS_GZIP -DLAMMPS_PNG -DLAMMPS_JPEG'
                env.JPG_LIB  = '-L../../src/STUBS/ -lmpi_stubs -ljpeg -lpng -lz'

                try {
                    compile_on_env('rbberger/lammps-testing:ubuntu_latest', 'gcc')
                    step([$class: 'GitHubCommitStatusSetter', contextSource: [$class: 'ManuallyEnteredCommitContextSource', context: 'continuous-integration/jenkins/serial'], statusResultSource: [$class: 'ConditionalStatusResultSource', results: [[$class: 'AnyBuildResult', message: 'serial build successful!', state: 'SUCCESS']]]])
                } catch (err) {
                    echo "Caught: ${err}"
                    currentBuild.result = 'FAILURE'
                    step([$class: 'GitHubCommitStatusSetter', contextSource: [$class: 'ManuallyEnteredCommitContextSource', context: 'continuous-integration/jenkins/serial'], statusResultSource: [$class: 'ConditionalStatusResultSource', results: [[$class: 'AnyBuildResult', message: 'serial build failed!', state: 'FAILURE']]]])
                }
            }
        },
        "Shared Library" : {
            node {
                unstash 'source'
                unstash 'libraries'

                step([$class: 'GitHubCommitStatusSetter', contextSource: [$class: 'ManuallyEnteredCommitContextSource', context: 'continuous-integration/jenkins/shlib'], statusResultSource: [$class: 'ConditionalStatusResultSource', results: [[$class: 'AnyBuildResult', message: 'running shlib compilation', state: 'PENDING']]]])

                env.COMP     = 'g++'
                env.MACH     = 'shlib'
                env.LMPFLAGS = '-sf off'
                env.LMP_INC  = '-I../../src/STUBS -DFFT_KISSFFT -DLAMMPS_GZIP -DLAMMPS_PNG -DLAMMPS_JPEG'
                env.JPG_LIB  = '-L../../src/STUBS/ -lmpi_stubs -ljpeg -lpng -lz'

                try {
                    compile_on_env('rbberger/lammps-testing:ubuntu_latest', 'gcc')
                    step([$class: 'GitHubCommitStatusSetter', contextSource: [$class: 'ManuallyEnteredCommitContextSource', context: 'continuous-integration/jenkins/shlib'], statusResultSource: [$class: 'ConditionalStatusResultSource', results: [[$class: 'AnyBuildResult', message: 'shlib build successful!', state: 'SUCCESS']]]])
                } catch (err) {
                    echo "Caught: ${err}"
                    currentBuild.result = 'FAILURE'
                    step([$class: 'GitHubCommitStatusSetter', contextSource: [$class: 'ManuallyEnteredCommitContextSource', context: 'continuous-integration/jenkins/shlib'], statusResultSource: [$class: 'ConditionalStatusResultSource', results: [[$class: 'AnyBuildResult', message: 'shlib build failed!', state: 'FAILURE']]]])
                }
            }

        },
        "OpenMPI binary" : {
            node {
                step([$class: 'GitHubCommitStatusSetter', contextSource: [$class: 'ManuallyEnteredCommitContextSource', context: 'continuous-integration/jenkins/openmpi'], statusResultSource: [$class: 'ConditionalStatusResultSource', results: [[$class: 'AnyBuildResult', message: 'running shlib compilation', state: 'PENDING']]]])

                unstash 'source'
                unstash 'libraries'

                env.COMP     = 'mpicxx'
                env.MACH     = 'mpi'
                env.MPICMD   = 'mpirun -np 4'
                env.LMPFLAGS = '-sf off'
                env.LMP_INC  = '-DFFT_KISSFFT -DLAMMPS_GZIP -DLAMMPS_PNG -DLAMMPS_JPEG -DLAMMPS_SMALLSMALL'
                env.JPG_LIB  = '-ljpeg -lpng -lz'

                try {
                    compile_on_env('rbberger/lammps-testing:ubuntu_latest', 'gcc')
                    step([$class: 'GitHubCommitStatusSetter', contextSource: [$class: 'ManuallyEnteredCommitContextSource', context: 'continuous-integration/jenkins/openmpi'], statusResultSource: [$class: 'ConditionalStatusResultSource', results: [[$class: 'AnyBuildResult', message: 'openmpi build successful!', state: 'SUCCESS']]]])
                } catch (err) {
                    echo "Caught: ${err}"
                    currentBuild.result = 'FAILURE'
                    step([$class: 'GitHubCommitStatusSetter', contextSource: [$class: 'ManuallyEnteredCommitContextSource', context: 'continuous-integration/jenkins/openmpi'], statusResultSource: [$class: 'ConditionalStatusResultSource', results: [[$class: 'AnyBuildResult', message: 'openmpi build failed!', state: 'FAILURE']]]])
                }
            }
        }
    )

    step([$class: 'WarningsPublisher', canComputeNew: false, consoleParsers: [[parserName: 'GNU Make + GNU C Compiler (gcc)']], defaultEncoding: '', excludePattern: '', healthy: '', includePattern: '', messagesPattern: '', unHealthy: ''])
}
