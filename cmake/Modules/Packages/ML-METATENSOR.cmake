include(FetchContent)

set(URL_BASE "https://github.com/lab-cosmo/metatensor/releases/download")

set(METATENSOR_CORE_VERSION "0.1.5")
FetchContent_Declare(metatensor
    URL ${URL_BASE}/metatensor-core-v${METATENSOR_CORE_VERSION}/metatensor-core-cxx-${METATENSOR_CORE_VERSION}.tar.gz
    URL_HASH SHA1=b895025056b6a479503356cb13f9b5011514fd72
)

message(STATUS "Fetching metatensor v${METATENSOR_CORE_VERSION} from github")
FetchContent_MakeAvailable(metatensor)


set(METATENSOR_TORCH_VERSION "0.4.0")
FetchContent_Declare(metatensor-torch
    URL ${URL_BASE}/metatensor-torch-v${METATENSOR_TORCH_VERSION}/metatensor-torch-cxx-${METATENSOR_TORCH_VERSION}.tar.gz
    URL_HASH SHA1=de8272d2b2e2ea74cef8da49d889e4029e2c9fe9
)

message(STATUS "Fetching metatensor-torch v${METATENSOR_TORCH_VERSION} from github")
FetchContent_MakeAvailable(metatensor-torch)


################ lammps target modifications ################

# Bring the `torch` target in scope to allow evaluation
# of cmake generator expression from `metatensor_torch`
find_package(Torch REQUIRED)

target_link_libraries(lammps PRIVATE metatensor_torch)

if (BUILD_OMP AND APPLE)
    message(FATAL_ERROR
        "Can not enable BUILD_OMP on Apple systems, since this results in "
        "two different versions of libiomp5.dylib (one from theb system and "
        "one from Torch) being linked to the final executable, which then "
        "segfaults"
    )
endif()
