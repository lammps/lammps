// This file is needed in order to get the linker language
// for the header only submodule.
// While we set the language properties in our normal cmake
// path it does not get set in the Trilinos environment.
// Furthermore, setting LINKER_LANGUAGE is only supported
// in CMAKE 3.19 and up.
void KOKKOS_SIMD_SRC_DUMMY_PREVENT_LINK_ERROR() {}
