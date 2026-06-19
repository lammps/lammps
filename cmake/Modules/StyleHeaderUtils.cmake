function(FindStyleHeaders path style_class file_pattern headers)
    file(GLOB files CONFIGURE_DEPENDS "${path}/${file_pattern}*.h")
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
    # only the array/mask-based categories that still expand the XXX_CLASS macro
    # at compile time need a style_*.h include list.  The name-keyed categories
    # are registered through the global registry from generated style_*.cpp and
    # no longer use style_*.h (see GenerateStyleSources).
    GenerateStyleHeader(${output_path} NBIN         nbin         ) # neighbor
    GenerateStyleHeader(${output_path} NPAIR        npair        ) # neighbor
    GenerateStyleHeader(${output_path} NSTENCIL     nstencil     ) # neighbor
    GenerateStyleHeader(${output_path} NTOPO        ntopo        ) # neighbor
endfunction(GenerateStyleHeaders)

# Generate a style_<style>.cpp registration translation unit by parsing the
# XxxStyle(key,Class) marker lines out of the collected style headers and
# emitting an explicit registration function that fills the global registry.
# one_arg = TRUE selects the creator signature (LAMMPS*); FALSE selects
# (LAMMPS*,int,char**).  Generated files are accumulated in the global property
# LAMMPS_STYLE_SOURCES for adding to the lammps target.
function(CreateStyleSource path property style macro base accessor regfunc one_arg includes)
    get_property(files GLOBAL PROPERTY ${property})

    # sort header files by basename for deterministic output
    set(pairs)
    foreach(fname ${files})
        get_filename_component(bname ${fname} NAME)
        list(APPEND pairs "${bname}|${fname}")
    endforeach()
    list(SORT pairs)

    set(include_block "")
    set(entry_block "")
    foreach(entry ${pairs})
        string(REGEX REPLACE "^[^|]*\\|" "" fname "${entry}")
        get_filename_component(bname ${fname} NAME)
        set(include_block "${include_block}#include \"${bname}\"\n")
        set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${fname}")

        # walk the marker block (between "#ifdef XXX_CLASS" and the matching
        # "#else"/"#endif").  Convert each XxxStyle(key,Class) marker into an
        # explicit add_builtin() and copy any nested preprocessor directives
        # verbatim so build-config-dependent markers stay guarded.
        file(STRINGS ${fname} all_lines)
        set(state 0)
        set(depth 0)
        foreach(line IN LISTS all_lines)
            if(state EQUAL 0)
                if("${line}" MATCHES "^#ifdef[ \t]+[A-Z_]+_CLASS")
                    set(state 1)
                    set(depth 1)
                endif()
            elseif(state EQUAL 1)
                if("${line}" MATCHES "^[ \t]*#[ \t]*if")
                    math(EXPR depth "${depth} + 1")
                    string(STRIP "${line}" dline)
                    set(entry_block "${entry_block}${dline}\n")
                elseif("${line}" MATCHES "^[ \t]*#[ \t]*elif")
                    if(depth GREATER 1)
                        string(STRIP "${line}" dline)
                        set(entry_block "${entry_block}${dline}\n")
                    endif()
                elseif("${line}" MATCHES "^[ \t]*#[ \t]*else")
                    if(depth EQUAL 1)
                        set(state 2)
                    else()
                        string(STRIP "${line}" dline)
                        set(entry_block "${entry_block}${dline}\n")
                    endif()
                elseif("${line}" MATCHES "^[ \t]*#[ \t]*endif")
                    math(EXPR depth "${depth} - 1")
                    if(depth EQUAL 0)
                        set(state 2)
                    else()
                        string(STRIP "${line}" dline)
                        set(entry_block "${entry_block}${dline}\n")
                    endif()
                elseif("${line}" MATCHES "${macro}\\(([^)]*)\\)")
                    # the style keyword never contains a comma, so split on the
                    # first comma; the remainder is the class name (may contain
                    # template <,>)
                    set(inner "${CMAKE_MATCH_1}")
                    string(FIND "${inner}" "," comma)
                    if(NOT comma LESS 0)
                        string(SUBSTRING "${inner}" 0 ${comma} key)
                        math(EXPR after "${comma} + 1")
                        string(SUBSTRING "${inner}" ${after} -1 cls)
                        string(STRIP "${key}" key)
                        string(STRIP "${cls}" cls)
                        set(entry_block "${entry_block}  ${accessor}().add_builtin(\"${key}\", &creator<${cls}>);\n")
                    endif()
                endif()
            endif()
        endforeach()
    endforeach()

    set(preamble "#include \"lammps.h\"\n#include \"creator_registry.h\"\n")
    foreach(inc ${includes})
        set(preamble "${preamble}#include \"${inc}\"\n")
    endforeach()

    if(one_arg)
        set(sig "LAMMPS *lmp")
        set(call "lmp")
    else()
        set(sig "LAMMPS *lmp, int narg, char **arg")
        set(call "lmp, narg, arg")
    endif()

    set(content "// auto-generated by the LAMMPS build system from style markers. DO NOT EDIT.\n\n")
    set(content "${content}${preamble}\n${include_block}\nnamespace LAMMPS_NS {\n\n")
    set(content "${content}template <typename T> static ${base} *creator(${sig})\n{\n  return new T(${call});\n}\n\n")
    set(content "${content}void ${regfunc}()\n{\n${entry_block}}\n\n")
    set(content "${content}}    // namespace LAMMPS_NS\n")

    file(WRITE "${path}/style_${style}.cpp.tmp" "${content}")
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different "${path}/style_${style}.cpp.tmp" "${path}/style_${style}.cpp")
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${path}/style_${style}.cpp")
    set_property(GLOBAL APPEND PROPERTY LAMMPS_STYLE_SOURCES "${path}/style_${style}.cpp")
endfunction(CreateStyleSource)

function(GenerateStyleSources output_path)
    message(STATUS "Generating style source files...")
    CreateStyleSource(${output_path} PAIR     pair     PairStyle     Pair     "Force::pair_styles"     register_pair_styles     TRUE "force.h;pair.h")
    CreateStyleSource(${output_path} BOND     bond     BondStyle     Bond     "Force::bond_styles"     register_bond_styles     TRUE "force.h;bond.h")
    CreateStyleSource(${output_path} ANGLE    angle    AngleStyle    Angle    "Force::angle_styles"    register_angle_styles    TRUE "force.h;angle.h")
    CreateStyleSource(${output_path} DIHEDRAL dihedral DihedralStyle Dihedral "Force::dihedral_styles" register_dihedral_styles TRUE "force.h;dihedral.h")
    CreateStyleSource(${output_path} IMPROPER improper ImproperStyle Improper "Force::improper_styles" register_improper_styles TRUE "force.h;improper.h")
    CreateStyleSource(${output_path} KSPACE   kspace   KSpaceStyle   KSpace   "Force::kspace_styles"   register_kspace_styles   TRUE "force.h;kspace.h")
    CreateStyleSource(${output_path} FIX       fix       FixStyle       Fix       "Modify::fix_styles"       register_fix_styles       FALSE "modify.h;fix.h")
    CreateStyleSource(${output_path} COMPUTE   compute   ComputeStyle   Compute   "Modify::compute_styles"   register_compute_styles   FALSE "modify.h;compute.h")
    CreateStyleSource(${output_path} INTEGRATE integrate IntegrateStyle Integrate "Update::integrate_styles" register_integrate_styles FALSE "update.h;integrate.h")
    CreateStyleSource(${output_path} MINIMIZE  minimize  MinimizeStyle  Min       "Update::minimize_styles"  register_minimize_styles  TRUE  "update.h;min.h")
    CreateStyleSource(${output_path} REGION    region    RegionStyle    Region    "Domain::region_styles"    register_region_styles    FALSE "domain.h;region.h")
    CreateStyleSource(${output_path} DUMP      dump      DumpStyle      Dump      "Output::dump_styles"      register_dump_styles      FALSE "output.h;dump.h")
    CreateStyleSource(${output_path} COMMAND   command   CommandStyle   Command   "Input::command_styles"    register_command_styles   TRUE  "input.h;command.h")
    CreateStyleSource(${output_path} ATOM_VEC  atom      AtomStyle      AtomVec   "Atom::avec_styles"        register_atom_styles      TRUE  "atom.h;atom_vec.h")
    CreateStyleSource(${output_path} BODY      body      BodyStyle      Body      "AtomVecBody::body_styles" register_body_styles      FALSE "atom_vec_body.h;body.h")
    CreateStyleSource(${output_path} READER    reader    ReaderStyle    Reader    "ReadDump::reader_styles"  register_reader_styles    TRUE  "read_dump.h;reader.h")
endfunction(GenerateStyleSources)

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
    file(GLOB files CONFIGURE_DEPENDS "${path}/${file_pattern}*.h")
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

# append "package_styles().<member>[\"key\"] = \"<package>\";" entries for one
# style category to out_var, parsed from all (enabled or not) package headers
# in the given PKG* global property.  Only the keyword and the package (the
# header's parent directory) matter here, so preprocessor directives and the
# class name are ignored.
function(AppendPackageEntries property macro member out_var)
    get_property(files GLOBAL PROPERTY ${property})
    set(pairs)
    foreach(fname ${files})
        get_filename_component(bname ${fname} NAME)
        list(APPEND pairs "${bname}|${fname}")
    endforeach()
    list(SORT pairs)

    set(block "${${out_var}}")
    foreach(entry ${pairs})
        string(REGEX REPLACE "^[^|]*\\|" "" fname "${entry}")
        get_filename_component(pkgdir ${fname} DIRECTORY)
        get_filename_component(pkgdir ${pkgdir} NAME)
        set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${fname}")
        file(STRINGS ${fname} marker_lines REGEX "${macro}\\(")
        foreach(line ${marker_lines})
            string(REGEX MATCH "${macro}\\(([^)]*)\\)" _m "${line}")
            if(NOT _m)
                continue()
            endif()
            set(inner "${CMAKE_MATCH_1}")
            string(FIND "${inner}" "," comma)
            if(comma LESS 0)
                continue()
            endif()
            string(SUBSTRING "${inner}" 0 ${comma} key)
            string(STRIP "${key}" key)
            set(block "${block}  package_styles().${member}[\"${key}\"] = \"${pkgdir}\";\n")
        endforeach()
    endforeach()
    set(${out_var} "${block}" PARENT_SCOPE)
endfunction(AppendPackageEntries)

# generate the single global package_registry.cpp (style keyword -> package name)
function(GeneratePackageRegistry output_path)
    message(STATUS "Generating package registry...")
    set(entries "")
    AppendPackageEntries(PKGANGLE     AngleStyle     angle_styles     entries)
    AppendPackageEntries(PKGATOM_VEC  AtomStyle      atom_styles      entries)
    AppendPackageEntries(PKGBODY      BodyStyle      body_styles      entries)
    AppendPackageEntries(PKGBOND      BondStyle      bond_styles      entries)
    AppendPackageEntries(PKGCOMMAND   CommandStyle   command_styles   entries)
    AppendPackageEntries(PKGCOMPUTE   ComputeStyle   compute_styles   entries)
    AppendPackageEntries(PKGDIHEDRAL  DihedralStyle  dihedral_styles  entries)
    AppendPackageEntries(PKGDUMP      DumpStyle      dump_styles      entries)
    AppendPackageEntries(PKGFIX       FixStyle       fix_styles       entries)
    AppendPackageEntries(PKGIMPROPER  ImproperStyle  improper_styles  entries)
    AppendPackageEntries(PKGINTEGRATE IntegrateStyle integrate_styles entries)
    AppendPackageEntries(PKGKSPACE    KSpaceStyle    kspace_styles    entries)
    AppendPackageEntries(PKGMINIMIZE  MinimizeStyle  minimize_styles  entries)
    AppendPackageEntries(PKGPAIR      PairStyle      pair_styles      entries)
    AppendPackageEntries(PKGREADER    ReaderStyle    reader_styles    entries)
    AppendPackageEntries(PKGREGION    RegionStyle    region_styles    entries)

    set(content "// auto-generated by the LAMMPS build system from style markers. DO NOT EDIT.\n\n")
    set(content "${content}#include \"package_registry.h\"\n\n")
    set(content "${content}#include \"lmptype.h\"\n\n")
    set(content "${content}namespace LAMMPS_NS {\n\n")
    set(content "${content}void _noopt register_package_styles()\n{\n${entries}}\n\n")
    set(content "${content}}    // namespace LAMMPS_NS\n")

    file(WRITE "${output_path}/package_registry.cpp.tmp" "${content}")
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different "${output_path}/package_registry.cpp.tmp" "${output_path}/package_registry.cpp")
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${output_path}/package_registry.cpp")
    set_property(GLOBAL APPEND PROPERTY LAMMPS_STYLE_SOURCES "${output_path}/package_registry.cpp")
endfunction(GeneratePackageRegistry)

