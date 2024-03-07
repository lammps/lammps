include(FetchContent)

set(URL_BASE "https://github.com/lab-cosmo/metatensor/releases/download")

set(METATENSOR_CORE_VERSION "0.1.4")
FetchContent_Declare(metatensor
    URL ${URL_BASE}/metatensor-core-v${METATENSOR_CORE_VERSION}/metatensor-core-cxx-${METATENSOR_CORE_VERSION}.tar.gz
    URL_HASH SHA1=fd896b3f63911761e82ff56be612a8a19143ce6d
)

message(STATUS "Fetching metatensor v${METATENSOR_CORE_VERSION} from github")
FetchContent_MakeAvailable(metatensor)


set(METATENSOR_TORCH_VERSION "0.3.0")
FetchContent_Declare(metatensor-torch
    URL ${URL_BASE}/metatensor-torch-v${METATENSOR_TORCH_VERSION}/metatensor-torch-cxx-${METATENSOR_TORCH_VERSION}.tar.gz
    URL_HASH SHA1=a0c4891776766a3691106ba5c7131481291d11a1
)

message(STATUS "Fetching metatensor-torch v${METATENSOR_TORCH_VERSION} from github")
FetchContent_MakeAvailable(metatensor-torch)


################ lammps target modifications ################

target_link_libraries(lammps PRIVATE metatensor_torch)

if (BUILD_OMP AND APPLE)
    message(FATAL_ERROR
        "Can not enable BUILD_OMP on Apple systems, since this results in "
        "two different versions of libiomp5.dylib (one from theb system and "
        "one from Torch) being linked to the final executable, which then "
        "segfaults"
    )
endif()
