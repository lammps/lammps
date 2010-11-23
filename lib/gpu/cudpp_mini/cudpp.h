// -------------------------------------------------------------
// CUDPP -- CUDA Data Parallel Primitives library
// -------------------------------------------------------------
// $Revision$
// $Date$
// ------------------------------------------------------------- 
// This source code is distributed under the terms of license.txt in
// the root directory of this source distribution.
// ------------------------------------------------------------- 
 
/**
 * @file
 * cudpp.h
 * 
 * @brief Main library header file.  Defines public interface.
 *
 * The CUDPP public interface is a C-only interface to enable 
 * linking with code written in other languages (e.g. C, C++, 
 * and Fortran).  While the internals of CUDPP are not limited 
 * to C (C++ features are used), the public interface is 
 * entirely C (thus it is declared "extern C").
 */

/**
 * \mainpage
 *
 * \section introduction Introduction
 * 
 * CUDPP is the CUDA Data Parallel Primitives Library. CUDPP is a
 * library of data-parallel algorithm primitives such as 
 * parallel-prefix-sum ("scan"), parallel sort and parallel reduction. 
 * Primitives such as these are important building blocks for a wide 
 * variety of data-parallel algorithms, including sorting, stream 
 * compaction, and building data structures such as trees and 
 * summed-area tables.
 *
 * \section overview Overview Presentation
 * 
 * A brief set of slides that describe the features, design principles,
 * applications and impact of CUDPP is available here:
 * <a href="http://cudpp.googlecode.com/svn/trunk/cudpp/doc/CUDPP_slides.pdf">CUDPP Presentation</a>.
 *
 * \section homepage Homepage
 * Homepage for CUDPP: http://code.google.com/p/cudpp
 * 
 * Announcements and discussion of CUDPP are hosted on the
 * <a href="http://groups.google.com/group/cudpp?hl=en">CUDPP Google Group</a>.
 * 
 * \section getting-started Getting Started with CUDPP
 *
 * You may want to start by browsing the \link publicInterface CUDPP Public 
 * Interface\endlink. For information on building CUDPP, see 
 * \ref building-cudpp "Building CUDPP".
 *
 * The "apps" subdirectory included with CUDPP has a few source code samples 
 * that use CUDPP:
 * - \ref example_simpleCUDPP "simpleCUDPP", a simple example of using 
 * cudppScan()
 * - satGL, an example of using cudppMultiScan() to generate a summed-area 
 * table (SAT) of a scene rendered in real time.  The SAT is then used to simulate 
 * depth of field blur.
 * - cudpp_testrig, a comprehensive test application for all the functionality 
 * of CUDPP
 *
 * We have also provided a code walkthrough of the 
 * \ref example_simpleCUDPP "simpleCUDPP" example.
 *
 * \section getting-help Getting Help and Reporting Problems
 *
 * To get help using CUDPP, please use the 
 * <a href="http://groups.google.com/group/cudpp?hl=en">CUDPP Google Group</a>.
 *
 * To report CUDPP bugs or request features, you may use either the above 
 * CUDPP Google Group, or you can file an issue directly using 
 * <a href="http://code.google.com/p/cudpp/issues/list">Google Code</a>.
 *
 * \section release-notes Release Notes
 *
 * For specific release details see the \ref changelog "Change Log".
 *
 * This release (1.1.1) is a bugfix release to CUDPP 1.1 that includes
 * fixes to support CUDA 3.0 and the new NVIDIA Fermi architecture, 
 * including GeForce 400 series and Tesla 20 series GPUs.  It also has
 * bug fixes for 64-bit OSes.
 * 
 * \section opSys Operating System Support
 * 
 * This release (1.1.1) has been thoroughly tested on the following OSes.
 * - Windows XP (32-bit) (CUDA 2.2, 3.0)
 * - Windows 7 (64-bit) (CUDA 3.0)
 * - Redhat Enterprise Linux 5 (64-bit) (CUDA 3.0)
 * - and Mac OS X 10.6 (Snow Leopard, 64-bit) (CUDA 3.0)
 *
 * We expect CUDPP to build and run correctly on other flavors of Linux 
 * and Windows, but these are not actively tested by the developers at 
 * this time.
 *
 * Notes: CUDPP is not compatible with CUDA 2.1.  A compiler bug in 2.1 
 * causes the compiler to crash.  Also, starting with CUDPP 1.1.1, we are 
 * no longer testing CUDA device emulation, because it is deprecated in 
 * CUDA 3.0 and will be removed from future CUDA versions.  
 *
 * \section cuda CUDA
 * CUDPP is implemented in
 * <a href="http://developer.nvidia.com/cuda">CUDA C/C++</a>. It requires the 
 * CUDA Toolkit version 2.2 or later.  Please see the NVIDIA 
 * <a href="http://developer.nvidia.com/cuda">CUDA</a> homepage to download 
 * CUDA as well as the CUDA Programming Guide and CUDA SDK, which includes many 
 * CUDA code examples.  Some of the samples in the CUDA SDK (including 
 * "marchingCubes", "lineOfSight", and radixSort) also use CUDPP.
 *
 * \section design-goals Design Goals
 * Design goals for CUDPP include:
 * 
 * - Performance. We aim to provide best-of-class performance for our
 *   primitives. We welcome suggestions and contributions that will improve 
 *   CUDPP performance. We also want to provide primitives that can be easily 
 *   benchmarked, and compared against other implementations on GPUs and other 
 *   processors.
 * - Modularity. We want our primitives to be easily included in other
 *   applications. To that end we have made the following design decisions:
 *   - CUDPP is provided as a library that can link against other applications. 
 *   - CUDPP calls run on the GPU on GPU data. Thus they can be used
 *     as standalone calls on the GPU (on GPU data initialized by the 
 *     calling application) and, more importantly, as GPU components in larger 
 *     CPU/GPU applications.
 *   - CUDPP is implemented as 4 layers:
 *     -# The \link publicInterface Public Interface\endlink is the external 
 *        library interface, which is the intended entry point for most 
 *        applications. The public interface calls into the 
 *        \link cudpp_app Application-Level API\endlink.
 *     -# The \link cudpp_app Application-Level API\endlink comprises functions
 *        callable from CPU code. These functions execute code jointly on the 
 *        CPU (host) and the GPU by calling into the 
 *        \link cudpp_kernel Kernel-Level API\endlink below them.
 *     -# The \link cudpp_kernel Kernel-Level API\endlink comprises functions
 *        that run entirely on the GPU across an entire grid of thread blocks.  
 *        These functions may call into the \link cudpp_cta CTA-Level API\endlink 
 *        below them.
 *     -# The \link cudpp_cta CTA-Level API\endlink comprises functions that run 
 *        entirely on the GPU within a single Cooperative Thread Array (CTA, 
 *        aka thread block). These are low-level functions that implement core 
 *        data-parallel algorithms, typically by processing data within shared 
 *        (CUDA \c __shared__) memory.
 *
 * Programmers may use any of the lower three CUDPP layers in their own 
 * programs by building the source directly into their application.  However, 
 * the typical usage of CUDPP is to link to the library and invoke functions in 
 * the CUDPP \link publicInterface Public Interface\endlink, as in the 
 * \ref example_simpleCUDPP "simpleCUDPP", satGL, and cudpp_testrig application 
 * examples included in the CUDPP distribution.
 *
 * In the future, if and when CUDA supports building device-level libraries, we 
 * hope to enhance CUDPP to ease the use of CUDPP internal algorithms at all 
 * levels.
 *
 * \subsection uses Use Cases
 * We expect the normal use of CUDPP will be in one of two ways:
 * -# Linking the CUDPP library against another application. 
 * -# Running our "test" application, cudpp_testrig, that exercises
 *   CUDPP functionality.
 *
 * \section references References
 * The following publications describe work incorporated in CUDPP.
 * 
 * - Mark Harris, Shubhabrata Sengupta, and John D. Owens. "Parallel Prefix Sum (Scan) with CUDA". In Hubert Nguyen, editor, <i>GPU Gems 3</i>, chapter 39, pages 851&ndash;876. Addison Wesley, August 2007. http://graphics.idav.ucdavis.edu/publications/print_pub?pub_id=916
 * - Shubhabrata Sengupta, Mark Harris, Yao Zhang, and John D. Owens. "Scan Primitives for GPU Computing". In <i>Graphics Hardware 2007</i>, pages 97&ndash;106, August 2007. http://graphics.idav.ucdavis.edu/publications/print_pub?pub_id=915
 * - Shubhabrata Sengupta, Mark Harris, and Michael Garland. "Efficient parallel scan algorithms for GPUs". NVIDIA Technical Report NVR-2008-003, December 2008. http://mgarland.org/papers.html#segscan-tr
 * - Nadathur Satish, Mark Harris, and Michael Garland. "Designing Efficient Sorting Algorithms for Manycore GPUs". In <i>Proceedings of the 23rd IEEE International Parallel & Distributed Processing Symposium</i>, May 2009. http://mgarland.org/papers.html#gpusort
 * - Stanley Tzeng, Li-Yi Wei. "Parallel White Noise Generation on a GPU via Cryptographic Hash". In <i>Proceedings of the 2008 Symposium on Interactive 3D Graphics and Games</i>, pages 79&ndash;87, February 2008. http://research.microsoft.com/apps/pubs/default.aspx?id=70502
 *
 * Many researchers are using CUDPP in their work, and there are many publications 
 * that have used it \ref cudpp_refs "(references)". If your work uses CUDPP, please 
 * let us know by sending us a reference (preferably in BibTeX format) to your work.
 * 
 * \section citing Citing CUDPP
 *
 * If you make use of CUDPP primitives in your work and want to cite
 * CUDPP (thanks!), we would prefer for you to cite the appropriate 
 * papers above, since they form the core of CUDPP. To be more specific, 
 * the GPU Gems paper describes (unsegmented) scan, multi-scan for 
 * summed-area tables, and stream compaction. The NVIDIA technical report 
 * describes the current scan and segmented scan algorithms used in the 
 * library, and the Graphics Hardware paper describes an earlier 
 * implementation of segmented scan, quicksort, and sparse matrix-vector 
 * multiply. The IPDPS paper describes the radix sort used in CUDPP, and 
 * the I3D paper describes the random number generation algorithm.
 *
 * \section credits Credits
 * \subsection developers CUDPP Developers
 * - <a href="http://www.markmark.net">Mark Harris</a>, NVIDIA Corporation
 * - <a href="http://www.ece.ucdavis.edu/~jowens/">John D. Owens</a>, University of California, Davis
 * - <a href="http://graphics.cs.ucdavis.edu/~shubho/">Shubho Sengupta</a>, University of California, Davis
 * - Stanley Tzeng,   University of California, Davis
 * - <a href="http://www.ece.ucdavis.edu/~yaozhang/">Yao Zhang</a>,       University of California, Davis
 * - <a href="http://www.ece.ucdavis.edu/~aaldavid/">Andrew Davidson</a>, University of California, Davis (formerly Louisiana State University)
 * 
 * \subsection contributors Other CUDPP Contributors
 * - <a href="http://www.eecs.berkeley.edu/~nrsatish/">Nadatur Satish</a>,  University of California, Berkeley
 *
 * \subsection acknowledgments Acknowledgments
 *
 * Thanks to Jim Ahrens, Timo Aila, Nathan Bell, Ian Buck, Guy Blelloch, 
 * Jeff Bolz, Michael Garland, Jeff Inman, Eric Lengyel, Samuli Laine, 
 * David Luebke, Pat McCormick, and Richard Vuduc for their contributions 
 * during the development of this library. 
 * 
 * CUDPP Developers from UC Davis thank their funding agencies:
 * - Department of Energy Early Career Principal Investigator Award
 *   DE-FG02-04ER25609
 * - SciDAC Institute for Ultrascale Visualization (http://www.iusv.org/)
 * - Los Alamos National Laboratory
 * - National Science Foundation (grant 0541448)
 * - Generous hardware donations from NVIDIA
 *
 * \section license-overview CUDPP Copyright and Software License
 * CUDPP is copyright The Regents of the University of California, Davis campus 
 * and NVIDIA Corporation.  The library, examples, and all source code are 
 * released under the BSD license, designed to encourage reuse of this software 
 * in other projects, both commercial and non-commercial.  For details, please 
 * see the \ref license page. 
 * 
 * Note that prior to release 1.1 of CUDPP, the license used was a modified
 * BSD license.  With release 1.1, this license was replaced with the pure BSD
 * license to facilitate the use of open source hosting of the code.
 */

/**
 * @page license CUDPP License
 *
 * \section licenseBSD CUDPP License
 *
 * CUDPP is released under the 
 * <a href="http://www.opensource.org/licenses/bsd-license.php">BSD license</a>.
 * 
 * @include license.txt
 *
 */

/** 
 * @page changelog CUDPP Change Log
 *
 * @include changelog.txt
 */

/** 
 * @page cudpp_refs Publications that use CUDPP
 *
 * @htmlinclude doc/bib/cudpp_refs.html
 */

/** 
 * @page cudpp_refs_bib Bibliography for publications that use CUDPP
 *
 * @htmlinclude doc/bib/cudpp_refs_bib.html
 */

/**
 * @page building-cudpp Building CUDPP
 *
 * CUDPP has currently been tested in Windows XP, Windows Vista, Mac OS X 
 * and Linux.  See \ref release-notes for release specific platform support.
 *
 * \section build-win32 Building CUDPP on Windows XP
 *
 * CUDPP can be built using either or MSVC 8 (2005) or MSVC 9 (2008).  To 
 * build, open cudpp/cudpp.sln. Then you can build the library 
 * using the "build" command as you would with any other workspace. There are 
 * four configurations: debug, release, emudebug, and emurelease.  The first 
 * two are self-explanatory.  The second two are built to use CUDA device 
 * emulation, meaning they will be run (slowly) on the CPU.
 *
 * \section build-linux Building CUDPP on Linux and Mac OS X
 *
 * CUDPP can be built using standard g++ and Make tools on Linux, by typing 
 * "make" in the "cudpp/" subdirectory.  Before building CUDPP, you should 
 * first build the CUDA Utility Library (libcutil) by typing "make; make dbg=1" 
 * in the "common/" subdirectory.  This will generate libcutil.a and 
 * libcutilD.a.  
 * 
 * The makefile for CUDPP and all sample applications take the optional 
 * arguments "emu=1" and "dbg=1".  The former builds CUDPP for device emulation,
 * and the latter for debugging. The two flags can be combined. "verbose=1"
 * can be used to see all compiler output.
 *
 * \section build-apps Building CUDPP Sample Applications
 * 
 * The sample applications in the "apps/" subdirectory can be built exactly 
 * like CUDPP is--either by opening the appropriate .sln/.vcproj file in MSVC 
 * in Windows, or using "make" in Linux.
 *
 * On some Linux installations you will get linker errors relating to "-lXi"
 * and "-lXmu".  To fix this, you will need to install libXi and libXmu.  On 
 * Debian and Ubuntu, for example, you can simply run 
 *  "sudo apt-get install libxi-dev", and 
 *  "sudo apt-get install libxmu-dev"
 * 
 */

#ifndef __CUDPP_H__
#define __CUDPP_H__

#include <stdlib.h> // for size_t

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief CUDPP Result codes returned by CUDPP API functions.
 */
enum CUDPPResult
{
    CUDPP_SUCCESS = 0,                 /**< No error. */
    CUDPP_ERROR_INVALID_HANDLE,        /**< Specified handle (for example, 
                                            to a plan) is invalid. **/
    CUDPP_ERROR_ILLEGAL_CONFIGURATION, /**< Specified configuration is
                                            illegal. For example, an
                                            invalid or illogical
                                            combination of options. */
    CUDPP_ERROR_UNKNOWN = 9999         /**< Unknown or untraceable error. */
};

/** 
 * @brief Options for configuring CUDPP algorithms.
 * 
 * @see CUDPPConfiguration, cudppPlan, CUDPPAlgorithm
 */
enum CUDPPOption
{
    CUDPP_OPTION_FORWARD   = 0x1,  /**< Algorithms operate forward:
                                    * from start to end of input
                                    * array */
    CUDPP_OPTION_BACKWARD  = 0x2,  /**< Algorithms operate backward:
                                    * from end to start of array */
    CUDPP_OPTION_EXCLUSIVE = 0x4,  /**< Exclusive (for scans) - scan
                                    * includes all elements up to (but
                                    * not including) the current
                                    * element */
    CUDPP_OPTION_INCLUSIVE = 0x8,  /**< Inclusive (for scans) - scan
                                    * includes all elements up to and
                                    * including the current element */
    CUDPP_OPTION_CTA_LOCAL = 0x10, /**< Algorithm performed only on
                                    * the CTAs (blocks) with no
                                    * communication between blocks.
                                    * @todo Currently ignored. */
    CUDPP_OPTION_KEYS_ONLY = 0x20, /**< No associated value to a key 
                                    * (for global radix sort) */
    CUDPP_OPTION_KEY_VALUE_PAIRS = 0x40, /**< Each key has an associated value */
};


/** 
 * @brief Datatypes supported by CUDPP algorithms.
 *
 * @see CUDPPConfiguration, cudppPlan
 */
enum CUDPPDatatype
{
    CUDPP_CHAR,     //!< Character type (C char)
    CUDPP_UCHAR,    //!< Unsigned character (byte) type (C unsigned char)
    CUDPP_INT,      //!< Integer type (C int)
    CUDPP_UINT,     //!< Unsigned integer type (C unsigned int)
    CUDPP_FLOAT     //!< Float type (C float)
};

/** 
 * @brief Operators supported by CUDPP algorithms (currently scan and
 * segmented scan).
 *
 * These are all binary associative operators.
 *
 * @see CUDPPConfiguration, cudppPlan
 */
enum CUDPPOperator
{
    CUDPP_ADD,      //!< Addition of two operands
    CUDPP_MULTIPLY, //!< Multiplication of two operands
    CUDPP_MIN,      //!< Minimum of two operands
    CUDPP_MAX       //!< Maximum of two operands
};

/**
* @brief Algorithms supported by CUDPP.  Used to create appropriate plans using
* cudppPlan.
* 
* @see CUDPPConfiguration, cudppPlan
*/
enum CUDPPAlgorithm
{
    CUDPP_SCAN,              //!< Scan or prefix-sum
    CUDPP_SEGMENTED_SCAN,    //!< Segmented scan
    CUDPP_COMPACT,           //!< Stream compact
    CUDPP_REDUCE,            //!< Parallel reduction (NOTE: currently unimplemented)
    CUDPP_SORT_RADIX,        //!< Radix sort
    CUDPP_SPMVMULT,          //!< Sparse matrix-dense vector multiplication
    CUDPP_RAND_MD5,          //!< PseudoRandom Number Generator using MD5 hash algorithm
    CUDPP_ALGORITHM_INVALID, //!< Placeholder at end of enum
};

/**
* @brief Configuration struct used to specify algorithm, datatype,
* operator, and options when creating a plan for CUDPP algorithms.
*
* @see cudppPlan
*/
struct CUDPPConfiguration
{
    CUDPPAlgorithm algorithm; //!< The algorithm to be used
    CUDPPOperator  op;        //!< The numerical operator to be applied
    CUDPPDatatype  datatype;  //!< The datatype of the input arrays
    unsigned int   options;   //!< Options to configure the algorithm
};

#define CUDPP_INVALID_HANDLE 0xC0DABAD1
typedef size_t CUDPPHandle;

/* To use CUDPP as a static library, #define CUDPP_STATIC_LIB before 
 * including cudpp.h
 */
#define CUDPP_STATIC_LIB
#ifndef CUDPP_DLL
    #ifdef _WIN32
        #ifdef CUDPP_STATIC_LIB
            #define CUDPP_DLL
        #else
        #ifdef BUILD_DLL
            #define CUDPP_DLL __declspec(dllexport)
        #else
            #define CUDPP_DLL __declspec(dllimport)
        #endif
        #endif
    #else
        #define CUDPP_DLL
    #endif
#endif

// Plan allocation (for scan, sort, and compact)

CUDPP_DLL
CUDPPResult cudppPlan(CUDPPHandle        *planHandle, 
                      CUDPPConfiguration config, 
                      size_t             n, 
                      size_t             rows, 
                      size_t             rowPitch);

CUDPP_DLL
CUDPPResult cudppDestroyPlan(CUDPPHandle plan);

// Scan and sort algorithms

CUDPP_DLL
CUDPPResult cudppScan(CUDPPHandle planHandle,
                      void        *d_out, 
                      const void  *d_in, 
                      size_t      numElements);

CUDPP_DLL
CUDPPResult cudppMultiScan(CUDPPHandle planHandle,
                           void        *d_out, 
                           const void  *d_in, 
                           size_t      numElements,
                           size_t      numRows);

CUDPP_DLL
CUDPPResult cudppSegmentedScan(CUDPPHandle        planHandle,
                               void               *d_out, 
                               const void         *d_idata,
                               const unsigned int *d_iflags,
                               size_t             numElements);

CUDPP_DLL
CUDPPResult cudppCompact(CUDPPHandle        planHandle,
                         void               *d_out, 
                         size_t             *d_numValidElements,
                         const void         *d_in, 
                         const unsigned int *d_isValid,
                         size_t             numElements);

CUDPP_DLL
CUDPPResult cudppSort(CUDPPHandle planHandle,
                      void        *d_keys,                                          
                      void        *d_values,                                                                       
                      int         keybits,
                      size_t      numElements);

// Sparse matrix allocation

CUDPP_DLL
CUDPPResult cudppSparseMatrix(CUDPPHandle        *sparseMatrixHandle, 
                              CUDPPConfiguration config, 
                              size_t             n, 
                              size_t             rows, 
                              const void         *A,
                              const unsigned int *h_rowIndices,
                              const unsigned int *h_indices);

CUDPP_DLL
CUDPPResult cudppDestroySparseMatrix(CUDPPHandle sparseMatrixHandle);

// Sparse matrix-vector algorithms

CUDPP_DLL
CUDPPResult cudppSparseMatrixVectorMultiply(CUDPPHandle sparseMatrixHandle,
                                            void        *d_y,
                                            const void  *d_x);

// random number generation algorithms
CUDPP_DLL
CUDPPResult cudppRand(CUDPPHandle planHandle,void * d_out, size_t numElements);

CUDPP_DLL
CUDPPResult cudppRandSeed(const CUDPPHandle planHandle, unsigned int seed);

#ifdef __cplusplus
}
#endif

#endif

// Leave this at the end of the file
// Local Variables:
// mode:c++
// c-file-style: "NVIDIA"
// End:
