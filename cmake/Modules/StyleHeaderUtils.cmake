function(FindStyleHeaders path style_class file_pattern headers)
    file(GLOB files ${CONFIGURE_DEPENDS} "${path}/${file_pattern}*.h")
    get_property(hlist GLOBAL PROPERTY ${headers})

    foreach(file_name ${files})
        file(STRINGS ${file_name} is_style LIMIT_COUNT 1 REGEX ${style_class})
        if(is_style)
            list(APPEND hlist ${file_name})
        endif()
    endforeach()
    set_property(GLOBAL PROPERTY ${headers} "${hlist}")
endfunction(FindStyleHeaders)

function(AddStyleHeader path headers)
    get_property(hlist GLOBAL PROPERTY ${headers})
    list(APPEND hlist ${path})
    set_property(GLOBAL PROPERTY ${headers} "${hlist}")
endfunction(AddStyleHeader)

function(FindStyleHeadersExt path style_class extension headers sources)
    get_property(hlist GLOBAL PROPERTY ${headers})
    get_property(slist GLOBAL PROPERTY ${sources})
    set(ext_list)
    get_filename_component(abs_path "${path}" ABSOLUTE)

    foreach(file_name ${hlist})
        get_filename_component(basename ${file_name} NAME_WE)
        set(ext_file_name "${abs_path}/${basename}_${extension}.h")
        if(EXISTS "${ext_file_name}")
            file(STRINGS ${ext_file_name} is_style LIMIT_COUNT 1 REGEX ${style_class})
            if(is_style)
                list(APPEND ext_list ${ext_file_name})

                set(source_file_name "${abs_path}/${basename}_${extension}.cpp")
                if(EXISTS "${source_file_name}")
                    list(APPEND slist ${source_file_name})
                endif()
            endif()
        endif()
    endforeach()

    list(APPEND hlist ${ext_list})
    set_property(GLOBAL PROPERTY ${headers} "${hlist}")
    set_property(GLOBAL PROPERTY ${sources} "${slist}")
endfunction(FindStyleHeadersExt)

function(CreateStyleHeader path filename)
    set(temp "")
    if(ARGC GREATER 2)
        list(REMOVE_AT ARGV 0 1)
        set(header_list)
        foreach(FNAME ${ARGV})
            set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${FNAME}")
            get_filename_component(FNAME ${FNAME} NAME)
            list(APPEND header_list ${FNAME})
        endforeach()
        list(SORT header_list)
        foreach(FNAME ${header_list})
            set(temp "${temp}#include \"${FNAME}\"\n")
        endforeach()
    endif()
    file(WRITE "${path}/${filename}.tmp" "${temp}" )
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different "${path}/${filename}.tmp" "${path}/${filename}")
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${path}/${filename}")
endfunction(CreateStyleHeader)

function(GenerateStyleHeader path property style)
    get_property(files GLOBAL PROPERTY ${property})
    #message("${property} = ${files}")
    CreateStyleHeader("${path}" "style_${style}.h" ${files})
endfunction(GenerateStyleHeader)

function(RegisterNBinStyles search_path)
    FindStyleHeaders(${search_path} NBIN_CLASS      nbin_      NBIN      ) # nbin      ) # neighbor
endfunction(RegisterNBinStyles)

function(RegisterNPairStyles search_path)
    FindStyleHeaders(${search_path} NPAIR_CLASS     npair_     NPAIR     ) # npair     ) # neighbor
endfunction(RegisterNPairStyles)

function(RegisterNBinStyle path)
    AddStyleHeader(${path} NBIN)
endfunction(RegisterNBinStyle)

function(RegisterNPairStyle path)
    AddStyleHeader(${path} NPAIR)
endfunction(RegisterNPairStyle)

function(RegisterFixStyle path)
    AddStyleHeader(${path} FIX)
endfunction(RegisterFixStyle)

function(RegisterIntegrateStyle path)
    AddStyleHeader(${path} INTEGRATE)
endfunction(RegisterIntegrateStyle)

function(RegisterStyles search_path)
    FindStyleHeaders(${search_path} ANGLE_CLASS        angle_        ANGLE        ) # angle        ) # force
    FindStyleHeaders(${search_path} ATOM_CLASS         atom_vec_     ATOM_VEC     ) # atom         ) # atom      atom_vec_hybrid
    FindStyleHeaders(${search_path} BODY_CLASS         body_         BODY         ) # body         ) # atom_vec_body
    FindStyleHeaders(${search_path} BOND_CLASS         bond_         BOND         ) # bond         ) # force
    FindStyleHeaders(${search_path} COMMAND_CLASS      "[^.]"        COMMAND      ) # command      ) # input
    FindStyleHeaders(${search_path} COMPUTE_CLASS      compute_      COMPUTE      ) # compute      ) # modify
    FindStyleHeaders(${search_path} DIHEDRAL_CLASS     dihedral_     DIHEDRAL     ) # dihedral     ) # force
    FindStyleHeaders(${search_path} DUMP_CLASS         dump_         DUMP         ) # dump         ) # output    write_dump
    FindStyleHeaders(${search_path} FIX_CLASS          fix_          FIX          ) # fix          ) # modify
    FindStyleHeaders(${search_path} GRAN_SUB_MOD_CLASS gran_sub_mod_ GRAN_SUB_MOD ) # gran_sub_mod ) # granular_model
    FindStyleHeaders(${search_path} IMPROPER_CLASS     improper_     IMPROPER     ) # improper     ) # force
    FindStyleHeaders(${search_path} INTEGRATE_CLASS    "[^.]"        INTEGRATE    ) # integrate    ) # update
    FindStyleHeaders(${search_path} KSPACE_CLASS       "[^.]"        KSPACE       ) # kspace       ) # force
    FindStyleHeaders(${search_path} MINIMIZE_CLASS     min_          MINIMIZE     ) # minimize     ) # update
    FindStyleHeaders(${search_path} NBIN_CLASS         nbin_         NBIN         ) # nbin         ) # neighbor
    FindStyleHeaders(${search_path} NPAIR_CLASS        npair_        NPAIR        ) # npair        ) # neighbor
    FindStyleHeaders(${search_path} NSTENCIL_CLASS     nstencil_     NSTENCIL     ) # nstencil     ) # neighbor
    FindStyleHeaders(${search_path} NTOPO_CLASS        ntopo_        NTOPO        ) # ntopo        ) # neighbor
    FindStyleHeaders(${search_path} PAIR_CLASS         pair_         PAIR         ) # pair         ) # force
    FindStyleHeaders(${search_path} READER_CLASS       reader_       READER       ) # reader       ) # read_dump
    FindStyleHeaders(${search_path} REGION_CLASS       region_       REGION       ) # region       ) # domain
endfunction(RegisterStyles)

function(RegisterStylesExt search_path extension sources)
    FindStyleHeadersExt(${search_path} ANGLE_CLASS        ${extension}  ANGLE        ${sources})
    FindStyleHeadersExt(${search_path} ATOM_CLASS         ${extension}  ATOM_VEC     ${sources})
    FindStyleHeadersExt(${search_path} BODY_CLASS         ${extension}  BODY         ${sources})
    FindStyleHeadersExt(${search_path} BOND_CLASS         ${extension}  BOND         ${sources})
    FindStyleHeadersExt(${search_path} COMMAND_CLASS      ${extension}  COMMAND      ${sources})
    FindStyleHeadersExt(${search_path} COMPUTE_CLASS      ${extension}  COMPUTE      ${sources})
    FindStyleHeadersExt(${search_path} DIHEDRAL_CLASS     ${extension}  DIHEDRAL     ${sources})
    FindStyleHeadersExt(${search_path} DUMP_CLASS         ${extension}  DUMP         ${sources})
    FindStyleHeadersExt(${search_path} FIX_CLASS          ${extension}  FIX          ${sources})
    FindStyleHeadersExt(${search_path} GRAN_SUB_MOD_CLASS ${extension}  GRAN_SUB_MOD ${sources})
    FindStyleHeadersExt(${search_path} IMPROPER_CLASS     ${extension}  IMPROPER     ${sources})
    FindStyleHeadersExt(${search_path} INTEGRATE_CLASS    ${extension}  INTEGRATE    ${sources})
    FindStyleHeadersExt(${search_path} KSPACE_CLASS       ${extension}  KSPACE       ${sources})
    FindStyleHeadersExt(${search_path} MINIMIZE_CLASS     ${extension}  MINIMIZE     ${sources})
    FindStyleHeadersExt(${search_path} NBIN_CLASS         ${extension}  NBIN         ${sources})
    FindStyleHeadersExt(${search_path} NPAIR_CLASS        ${extension}  NPAIR        ${sources})
    FindStyleHeadersExt(${search_path} NSTENCIL_CLASS     ${extension}  NSTENCIL     ${sources})
    FindStyleHeadersExt(${search_path} NTOPO_CLASS        ${extension}  NTOPO        ${sources})
    FindStyleHeadersExt(${search_path} PAIR_CLASS         ${extension}  PAIR         ${sources})
    FindStyleHeadersExt(${search_path} READER_CLASS       ${extension}  READER       ${sources})
    FindStyleHeadersExt(${search_path} REGION_CLASS       ${extension}  REGION       ${sources})
endfunction(RegisterStylesExt)

function(GenerateStyleHeaders output_path)
    message(STATUS "Generating style headers...")
    GenerateStyleHeader(${output_path} ANGLE        angle        ) # force
    GenerateStyleHeader(${output_path} ATOM_VEC     atom         ) # atom      atom_vec_hybrid
    GenerateStyleHeader(${output_path} BODY         body         ) # atom_vec_body
    GenerateStyleHeader(${output_path} BOND         bond         ) # force
    GenerateStyleHeader(${output_path} COMMAND      command      ) # input
    GenerateStyleHeader(${output_path} COMPUTE      compute      ) # modify
    GenerateStyleHeader(${output_path} DIHEDRAL     dihedral     ) # force
    GenerateStyleHeader(${output_path} DUMP         dump         ) # output    write_dump
    GenerateStyleHeader(${output_path} FIX          fix          ) # modify
    GenerateStyleHeader(${output_path} GRAN_SUB_MOD gran_sub_mod ) # granular_model
    GenerateStyleHeader(${output_path} IMPROPER     improper     ) # force
    GenerateStyleHeader(${output_path} INTEGRATE    integrate    ) # update
    GenerateStyleHeader(${output_path} KSPACE       kspace       ) # force
    GenerateStyleHeader(${output_path} MINIMIZE     minimize     ) # update
    GenerateStyleHeader(${output_path} NBIN         nbin         ) # neighbor
    GenerateStyleHeader(${output_path} NPAIR        npair        ) # neighbor
    GenerateStyleHeader(${output_path} NSTENCIL     nstencil     ) # neighbor
    GenerateStyleHeader(${output_path} NTOPO        ntopo        ) # neighbor
    GenerateStyleHeader(${output_path} PAIR         pair         ) # force
    GenerateStyleHeader(${output_path} READER       reader       ) # read_dump
    GenerateStyleHeader(${output_path} REGION       region       ) # domain
endfunction(GenerateStyleHeaders)

function(DetectBuildSystemConflict lammps_src_dir)
  if(ARGC GREATER 1)
    list(REMOVE_AT ARGV 0)
    foreach(SRC_FILE ${ARGV})
        get_filename_component(FILENAME ${SRC_FILE} NAME)
        if(EXISTS ${lammps_src_dir}/${FILENAME})
            message(FATAL_ERROR "\n########################################################################\n"
                                  "Found package(s) installed by the make-based build system\n"
                                  "\n"
                                  "Please run\n"
                                  "make -C ${lammps_src_dir} no-all purge\n"
                                  "to uninstall\n"
                                  "########################################################################")
        endif()
    endforeach()
  endif()
endfunction(DetectBuildSystemConflict)


function(FindPackagesHeaders path style_class file_pattern headers)
    file(GLOB files ${CONFIGURE_DEPENDS} "${path}/${file_pattern}*.h")
    get_property(plist GLOBAL PROPERTY ${headers})

    foreach(file_name ${files})
        file(STRINGS ${file_name} is_style LIMIT_COUNT 1 REGEX ${style_class})
        if(is_style)
            list(APPEND plist ${file_name})
        endif()
    endforeach()
    set_property(GLOBAL PROPERTY ${headers} "${plist}")
endfunction(FindPackagesHeaders)

function(RegisterPackages search_path)
    FindPackagesHeaders(${search_path} ANGLE_CLASS     angle_     PKGANGLE     ) # angle     ) # force
    FindPackagesHeaders(${search_path} ATOM_CLASS      atom_vec_  PKGATOM_VEC  ) # atom      ) # atom      atom_vec_hybrid
    FindPackagesHeaders(${search_path} BODY_CLASS      body_      PKGBODY      ) # body      ) # atom_vec_body
    FindPackagesHeaders(${search_path} BOND_CLASS      bond_      PKGBOND      ) # bond      ) # force
    FindPackagesHeaders(${search_path} COMMAND_CLASS   "[^.]"     PKGCOMMAND   ) # command   ) # input
    FindPackagesHeaders(${search_path} COMPUTE_CLASS   compute_   PKGCOMPUTE   ) # compute   ) # modify
    FindPackagesHeaders(${search_path} DIHEDRAL_CLASS  dihedral_  PKGDIHEDRAL  ) # dihedral  ) # force
    FindPackagesHeaders(${search_path} DUMP_CLASS      dump_      PKGDUMP      ) # dump      ) # output    write_dump
    FindPackagesHeaders(${search_path} FIX_CLASS       fix_       PKGFIX       ) # fix       ) # modify
    FindPackagesHeaders(${search_path} IMPROPER_CLASS  improper_  PKGIMPROPER  ) # improper  ) # force
    FindPackagesHeaders(${search_path} INTEGRATE_CLASS "[^.]"     PKGINTEGRATE ) # integrate ) # update
    FindPackagesHeaders(${search_path} KSPACE_CLASS    "[^.]"     PKGKSPACE    ) # kspace    ) # force
    FindPackagesHeaders(${search_path} MINIMIZE_CLASS  min_       PKGMINIMIZE  ) # minimize  ) # update
    FindPackagesHeaders(${search_path} NBIN_CLASS      nbin_      PKGNBIN      ) # nbin      ) # neighbor
    FindPackagesHeaders(${search_path} NPAIR_CLASS     npair_     PKGNPAIR     ) # npair     ) # neighbor
    FindPackagesHeaders(${search_path} NSTENCIL_CLASS  nstencil_  PKGNSTENCIL  ) # nstencil  ) # neighbor
    FindPackagesHeaders(${search_path} NTOPO_CLASS     ntopo_     PKGNTOPO     ) # ntopo     ) # neighbor
    FindPackagesHeaders(${search_path} PAIR_CLASS      pair_      PKGPAIR      ) # pair      ) # force
    FindPackagesHeaders(${search_path} READER_CLASS    reader_    PKGREADER    ) # reader    ) # read_dump
    FindPackagesHeaders(${search_path} REGION_CLASS    region_    PKGREGION    ) # region    ) # domain
endfunction(RegisterPackages)

function(CreatePackagesHeader path filename)
  set(temp "")
  if(ARGC GREATER 2)
    list(REMOVE_AT ARGV 0 1)
    foreach(FNAME ${ARGV})
      set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${FNAME}")
      get_filename_component(DNAME ${FNAME} DIRECTORY)
      get_filename_component(DNAME ${DNAME} NAME)
      get_filename_component(FNAME ${FNAME} NAME)
      set(temp "${temp}#undef PACKAGE\n#define PACKAGE \"${DNAME}\"\n")
      set(temp "${temp}#include \"${DNAME}/${FNAME}\"\n")
    endforeach()
  endif()
  file(WRITE "${path}/${filename}.tmp" "${temp}" )
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different "${path}/${filename}.tmp" "${path}/${filename}")
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${path}/${filename}")
endfunction(CreatePackagesHeader)

function(GeneratePackagesHeader path property style)
  get_property(files GLOBAL PROPERTY ${property})
  CreatePackagesHeader("${path}" "packages_${style}.h" ${files})
endfunction(GeneratePackagesHeader)

function(GeneratePackagesHeaders output_path)
    message(STATUS "Generating package headers...")
    GeneratePackagesHeader(${output_path} PKGANGLE      angle     ) # force
    GeneratePackagesHeader(${output_path} PKGATOM_VEC   atom      ) # atom      atom_vec_hybrid
    GeneratePackagesHeader(${output_path} PKGBODY       body      ) # atom_vec_body
    GeneratePackagesHeader(${output_path} PKGBOND       bond      ) # force
    GeneratePackagesHeader(${output_path} PKGCOMMAND    command   ) # input
    GeneratePackagesHeader(${output_path} PKGCOMPUTE    compute   ) # modify
    GeneratePackagesHeader(${output_path} PKGDIHEDRAL   dihedral  ) # force
    GeneratePackagesHeader(${output_path} PKGDUMP       dump      ) # output    write_dump
    GeneratePackagesHeader(${output_path} PKGFIX        fix       ) # modify
    GeneratePackagesHeader(${output_path} PKGIMPROPER   improper  ) # force
    GeneratePackagesHeader(${output_path} PKGINTEGRATE  integrate ) # update
    GeneratePackagesHeader(${output_path} PKGKSPACE     kspace    ) # force
    GeneratePackagesHeader(${output_path} PKGMINIMIZE   minimize  ) # update
    GeneratePackagesHeader(${output_path} PKGNBIN       nbin      ) # neighbor
    GeneratePackagesHeader(${output_path} PKGNPAIR      npair     ) # neighbor
    GeneratePackagesHeader(${output_path} PKGNSTENCIL   nstencil  ) # neighbor
    GeneratePackagesHeader(${output_path} PKGNTOPO      ntopo     ) # neighbor
    GeneratePackagesHeader(${output_path} PKGPAIR       pair      ) # force
    GeneratePackagesHeader(${output_path} PKGREADER     reader    ) # read_dump
    GeneratePackagesHeader(${output_path} PKGREGION     region    ) # domain
endfunction(GeneratePackagesHeaders)

