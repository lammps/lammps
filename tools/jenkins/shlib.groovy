node {
    def build_name = 'jenkins/shlib'

    step([$class: 'GitHubCommitStatusSetter', contextSource: [$class: 'ManuallyEnteredCommitContextSource', context: build_name], statusResultSource: [$class: 'ConditionalStatusResultSource', results: [[$class: 'AnyBuildResult', message: 'building...', state: 'PENDING']]]])

    stage 'Checkout'
    checkout scm

    env.CCACHE_DIR= pwd() + '/.ccache'
    env.COMP     = 'g++'
    env.MACH     = 'shlib'
    env.LMPFLAGS = '-sf off'
    env.LMP_INC  = '-I../../src/STUBS -DFFT_KISSFFT -DLAMMPS_GZIP -DLAMMPS_PNG -DLAMMPS_JPEG'
    env.JPG_LIB  = '-L../../src/STUBS/ -lmpi_stubs -ljpeg -lpng -lz'

    env.CC = 'gcc'
    env.CXX = 'g++'
    env.OMPI_CC = 'gcc'
    env.OMPI_CXX = 'g++'

    stage 'Setting up build environment'

    def envImage = docker.image('rbberger/lammps-testing:ubuntu_latest')

    try {
        docker.withRegistry('https://registry.hub.docker.com', 'docker-registry-login') {
            // ensure image is current
            envImage.pull()

            // use workaround (see https://issues.jenkins-ci.org/browse/JENKINS-34276)
            docker.image(envImage.imageName()).inside {
                sh 'ccache -C'
                sh 'ccache -M 5G'

                // clean up project directory
                sh '''
                make -C src clean-all
                make -C src yes-all
                make -C src no-lib
                '''

                stage 'Building libraries'

                sh '''
                make -C lib/colvars clean-all
                make -j 8 -C lib/colvars -f Makefile.g++ CXX="${COMP}"
                make -j 8 -C lib/poems -f Makefile.g++ CXX="${COMP}"
                make -j 8 -C lib/voronoi -f Makefile.g++ CXX="${COMP}"
                make -j 8 -C lib/awpmd -f Makefile.mpicc CC="${COMP}"
                make -j 8 -C lib/meam -f Makefile.gfortran CC=gcc F90=gfortran
                make -j 8 -C lib/h5md
                '''

                stage 'Enabling modules'

                sh '''
                make -C src yes-user-smd yes-user-molfile yes-compress yes-python
                make -C src yes-poems yes-voronoi yes-user-colvars yes-user-awpmd yes-meam
                make -C src yes-user-h5md
                '''

                stage 'Compiling'
                sh 'make -j 8 -C src mode=shlib serial MACH=serial MPICMD="${MPICMD}" CC="${COMP}" LINK="${COMP}" LMP_INC="${LMP_INC}" JPG_LIB="${JPG_LIB}" TAG="${TAG}-$CC" LMPFLAGS="${LMPFLAGS}"'
            }
        }
        step([$class: 'GitHubCommitStatusSetter', contextSource: [$class: 'ManuallyEnteredCommitContextSource', context: build_name], statusResultSource: [$class: 'ConditionalStatusResultSource', results: [[$class: 'AnyBuildResult', message: 'build successful!', state: 'SUCCESS']]]])
    } catch (err) {
        echo "Caught: ${err}"
        currentBuild.result = 'FAILURE'
        step([$class: 'GitHubCommitStatusSetter', contextSource: [$class: 'ManuallyEnteredCommitContextSource', context: build_name], statusResultSource: [$class: 'ConditionalStatusResultSource', results: [[$class: 'AnyBuildResult', message: 'build failed!', state: 'FAILURE']]]])
    }

    step([$class: 'WarningsPublisher', canComputeNew: false, consoleParsers: [[parserName: 'GNU Make + GNU C Compiler (gcc)']], defaultEncoding: '', excludePattern: '', healthy: '', includePattern: '', messagesPattern: '', unHealthy: ''])
}
