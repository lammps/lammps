# metatensor requires C++17 due to Torch requiring C++17
if(CMAKE_CXX_STANDARD LESS 17)
  message(FATAL_ERROR "The ML-METATOMIC package requires the C++ standard to
be set to at least C++17")
endif()

if (BUILD_OMP AND APPLE)
    message(FATAL_ERROR
        "Can not enable both BUILD_OMP and PGK_ML-METATOMIC on Apple systems, "
        "since this results in two different versions of the OpenMP library (one "
        "from the system and one from Torch) being linked to the final "
        "executable, which then crashes"
    )
endif()

# Bring the `torch` target in scope to allow evaluation
# of cmake generator expression from `metatensor_torch`
find_package(Torch REQUIRED)

# The caffe2::mkl target contains MKL_INCLUDE_DIR in it's
# INTERFACE_INCLUDE_DIRECTORIES even if MKL was not found, causing a build
# failure with "Imported target "torch" includes non-existent path" down the
# line. This code removes the missing path from INTERFACE_INCLUDE_DIRECTORIES,
# allowing the build to continue further.
if (TARGET caffe2::mkl)
    get_target_property(CAFFE2_MKL_INCLUDE_DIRECTORIES caffe2::mkl INTERFACE_INCLUDE_DIRECTORIES)
    set(MKL_INCLUDE_DIR_NOTFOUND "")
    foreach(_include_dir_ ${CAFFE2_MKL_INCLUDE_DIRECTORIES})
        if ("${_include_dir_}" MATCHES "MKL_INCLUDE_DIR-NOTFOUND")
            set(MKL_INCLUDE_DIR_NOTFOUND "${_include_dir_}")
        endif()
    endforeach()

    if (NOT "${MKL_INCLUDE_DIR_NOTFOUND}" STREQUAL "")
        list(REMOVE_ITEM CAFFE2_MKL_INCLUDE_DIRECTORIES "${MKL_INCLUDE_DIR_NOTFOUND}")
    endif()
    set_target_properties(caffe2::mkl PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${CAFFE2_MKL_INCLUDE_DIRECTORIES}"
    )
endif()

################ definition of metatensor and metatomic targets ################

set(METATENSOR_CORE_VERSION "0.1.14")
set(METATENSOR_CORE_SHA1 "9e21c48d9059d8a37618958d9d253220dedf7562")

set(METATENSOR_TORCH_VERSION "0.7.6")
set(METATENSOR_TORCH_SHA1 "5668f5088a42507e9ca4a7b723b3baac0286035c")

set(METATOMIC_TORCH_VERSION "0.1.1")
set(METATOMIC_TORCH_SHA1 "12b830c23674339fc185ce6e94e5a445416199ff")

set(DOWNLOAD_METATENSOR_DEFAULT ON)
find_package(metatensor_torch QUIET ${METATENSOR_TORCH_VERSION})
if (metatensor_torch_FOUND)
    set(DOWNLOAD_METATENSOR_DEFAULT OFF)
endif()

set(DOWNLOAD_METATOMIC_DEFAULT ON)
find_package(metatomic_torch QUIET ${METATOMIC_TORCH_VERSION})
if (metatomic_torch_FOUND)
    set(DOWNLOAD_METATOMIC_DEFAULT OFF)
endif()


option(DOWNLOAD_METATENSOR "Download metatensor package instead of using an already installed one" ${DOWNLOAD_METATENSOR_DEFAULT})
option(DOWNLOAD_METATOMIC "Download metatomic package instead of using an already installed one" ${DOWNLOAD_METATOMIC_DEFAULT})

if (DOWNLOAD_METATENSOR)
    include(FetchContent)

    set(URL_BASE "https://github.com/metatensor/metatensor/releases/download")
    FetchContent_Declare(metatensor
        URL ${URL_BASE}/metatensor-core-v${METATENSOR_CORE_VERSION}/metatensor-core-cxx-${METATENSOR_CORE_VERSION}.tar.gz
        URL_HASH SHA1=${METATENSOR_CORE_SHA1}
    )

    message(STATUS "Fetching metatensor v${METATENSOR_CORE_VERSION} from github")
    FetchContent_MakeAvailable(metatensor)

    FetchContent_Declare(metatensor-torch
        URL ${URL_BASE}/metatensor-torch-v${METATENSOR_TORCH_VERSION}/metatensor-torch-cxx-${METATENSOR_TORCH_VERSION}.tar.gz
        URL_HASH SHA1=${METATENSOR_TORCH_SHA1}
    )

    message(STATUS "Fetching metatensor-torch v${METATENSOR_TORCH_VERSION} from github")
    FetchContent_MakeAvailable(metatensor-torch)
else()
    # make sure to fail the configuration if cmake can not find metatensor-torch
    find_package(metatensor_torch REQUIRED ${METATENSOR_TORCH_VERSION})
endif()

if (DOWNLOAD_METATOMIC)
    include(FetchContent)

    set(URL_BASE "https://github.com/metatensor/metatomic/releases/download")
    FetchContent_Declare(metatomic-torch
        URL ${URL_BASE}/metatomic-torch-v${METATOMIC_TORCH_VERSION}/metatomic-torch-cxx-${METATOMIC_TORCH_VERSION}.tar.gz
        URL_HASH SHA1=${METATOMIC_TORCH_SHA1}
    )

    message(STATUS "Fetching metatomic-torch v${METATOMIC_TORCH_VERSION} from github")
    FetchContent_MakeAvailable(metatomic-torch)
else()
    # make sure to fail the configuration if cmake can not find metatomic-torch
    find_package(metatomic_torch REQUIRED ${METATOMIC_TORCH_VERSION})
endif()


################ lammps target modifications ################

target_link_libraries(lammps PUBLIC metatomic_torch metatensor_torch)
