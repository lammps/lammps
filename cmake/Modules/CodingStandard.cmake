# use default (or custom) Python executable, if version is sufficient
if(Python_VERSION VERSION_GREATER_EQUAL 3.6)
  set(Python3_EXECUTABLE ${Python_EXECUTABLE})
endif()
find_package(Python3 COMPONENTS Interpreter)

if(Python3_EXECUTABLE)
    if(Python3_VERSION VERSION_GREATER_EQUAL 3.6)
        add_custom_target(
          check-whitespace
          ${Python3_EXECUTABLE} ${LAMMPS_TOOLS_DIR}/coding_standard/whitespace.py .
          WORKING_DIRECTORY  ${LAMMPS_DIR}
          COMMENT "Check for whitespace errors")
        add_custom_target(
          check-homepage
          ${Python3_EXECUTABLE} ${LAMMPS_TOOLS_DIR}/coding_standard/homepage.py .
          WORKING_DIRECTORY  ${LAMMPS_DIR}
          COMMENT "Check for homepage URL errors")
        add_custom_target(
          check-permissions
          ${Python3_EXECUTABLE} ${LAMMPS_TOOLS_DIR}/coding_standard/permissions.py .
          WORKING_DIRECTORY  ${LAMMPS_DIR}
          COMMENT "Check for permission errors")
        add_custom_target(
          fix-whitespace
          ${Python3_EXECUTABLE} ${LAMMPS_TOOLS_DIR}/coding_standard/whitespace.py -f .
          WORKING_DIRECTORY  ${LAMMPS_DIR}
          COMMENT "Fix whitespace errors")
        add_custom_target(
          fix-homepage
          ${Python3_EXECUTABLE} ${LAMMPS_TOOLS_DIR}/coding_standard/homepage.py -f .
          WORKING_DIRECTORY  ${LAMMPS_DIR}
          COMMENT "Fix homepage URL errors")
        add_custom_target(
          fix-permissions
          ${Python3_EXECUTABLE} ${LAMMPS_TOOLS_DIR}/coding_standard/permissions.py -f .
          WORKING_DIRECTORY  ${LAMMPS_DIR}
          COMMENT "Fix permission errors")
    endif()
endif()
