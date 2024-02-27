include(FetchContent)
set(BUILD_METATENSOR_TORCH ON)
FetchContent_Declare(metatensor
    GIT_REPOSITORY https://github.com/lab-cosmo/metatensor
    GIT_TAG        ec74215ec48bb0dd902bf6f8fbb980d1d5e89117
    GIT_SHALLOW    ON
)

message(STATUS "Fetching metatensor from github")
FetchContent_MakeAvailable(metatensor)

################ lammps target modifications ################

# find_package(metatensor_torch)
target_link_libraries(lammps PRIVATE metatensor_torch)


if (BUILD_OMP AND APPLE)
    message(FATAL_ERROR
        "Can not enable BUILD_OMP on Apple systems, since this results in "
        "two different versions of libiomp5.dylib (one from theb system and "
        "one from Torch) being linked to the final executable, which then "
        "segfaults"
    )
endif()
