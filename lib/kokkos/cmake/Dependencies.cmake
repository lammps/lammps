TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
  SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
    #SubPackageName       Directory         Class    Req/Opt
    #
    # New Kokkos subpackages:
    Core                  core              PS       REQUIRED
    Containers            containers        PS       OPTIONAL
    Algorithms            algorithms        PS       OPTIONAL
    Example               example           EX       OPTIONAL
  )
