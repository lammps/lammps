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

set(METATENSOR_CORE_VERSION "0.1.17")
set(METATENSOR_CORE_SHA256 "42119e11908239915ccc187d7ca65449b461f1d4b5af4d6df1fb613d687da76a")

set(METATENSOR_TORCH_VERSION "0.8.1")
set(METATENSOR_TORCH_SHA256 "9da124e8e09dc1859700723a76ff29aef7a216b84a19d38746cc45bf45bc599b")

set(METATOMIC_TORCH_VERSION "0.1.5")
set(METATOMIC_TORCH_SHA256 "8ecd1587797fe1cf6b2162ddc10cc84c558fdfd55ab225bc5de4fe15ace8fc3d")

set(DOWNLOAD_METATENSOR_DEFAULT ON)
find_package(metatensor_torch ${METATENSOR_TORCH_VERSION} QUIET)
if (metatensor_torch_FOUND)
    set(DOWNLOAD_METATENSOR_DEFAULT OFF)
endif()

set(DOWNLOAD_METATOMIC_DEFAULT ON)
find_package(metatomic_torch ${METATOMIC_TORCH_VERSION} QUIET)
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
        URL_HASH SHA256=${METATENSOR_CORE_SHA256}
    )

    message(STATUS "Fetching metatensor v${METATENSOR_CORE_VERSION} from github")
    FetchContent_MakeAvailable(metatensor)

    FetchContent_Declare(metatensor-torch
        URL ${URL_BASE}/metatensor-torch-v${METATENSOR_TORCH_VERSION}/metatensor-torch-cxx-${METATENSOR_TORCH_VERSION}.tar.gz
        URL_HASH SHA256=${METATENSOR_TORCH_SHA256}
    )

    message(STATUS "Fetching metatensor-torch v${METATENSOR_TORCH_VERSION} from github")
    FetchContent_MakeAvailable(metatensor-torch)
else()
    # make sure to fail the configuration if cmake can not find metatensor-torch
    find_package(metatensor_torch ${METATENSOR_TORCH_VERSION} REQUIRED)
endif()

if (DOWNLOAD_METATOMIC)
    include(FetchContent)

    set(URL_BASE "https://github.com/metatensor/metatomic/releases/download")
    FetchContent_Declare(metatomic-torch
        URL ${URL_BASE}/metatomic-torch-v${METATOMIC_TORCH_VERSION}/metatomic-torch-cxx-${METATOMIC_TORCH_VERSION}.tar.gz
        URL_HASH SHA256=${METATOMIC_TORCH_SHA256}
    )

    message(STATUS "Fetching metatomic-torch v${METATOMIC_TORCH_VERSION} from github")
    FetchContent_MakeAvailable(metatomic-torch)
else()
    # make sure to fail the configuration if cmake can not find metatomic-torch
    find_package(metatomic_torch ${METATOMIC_TORCH_VERSION} REQUIRED)
endif()


################ lammps target modifications ################

target_link_libraries(lammps PUBLIC metatomic_torch metatensor_torch)
