# --------------------------------------------------------------------------- #
#    Custom CMake find module for Eigen3                                      #
#                                                                             #
#    This find module is provided because the user may not have Eigen3        #
#    configuration file installed in the system default directories, which    #
#    happens to be the case for some of our main developers and testers.      #
#                                                                             #
#    Accepts the following PATH:                                              #
#                                                                             #
#    - EIGEN3_ROOT - Custom Eigen3 root directory                             #
#                                                                             #
#    Provides (at least) the following variables:                             #
#                                                                             #
#    - Eigen3_FOUND - Whether Eigen3 was found or not                         #
#    - EIGEN3_INCLUDE_DIRS - Include directories                              #
#    - Eigen3::Eigen - A target to use with target_link_libraries()           #
#                                                                             #
#                              Niccolo' Iardella                              #
#                          Operations Research Group                          #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
include(FindPackageHandleStandardArgs)

# ----- Find the headers ---------------------------------------------------- #
find_path(EIGEN3_INCLUDE_DIR Eigen/Dense
          PATHS ${EIGEN3_ROOT}
          DOC "Eigen3 include directory")

# ----- Parse the version --------------------------------------------------- #
if (EIGEN3_INCLUDE_DIR)
    file(STRINGS
         "${EIGEN3_INCLUDE_DIR}/Eigen/src/Core/util/Macros.h"
         _eigen3_version_lines REGEX "#define EIGEN_(WORLD|MAJOR|MINOR)_VERSION")

    string(REGEX REPLACE ".*EIGEN_WORLD_VERSION *\([0-9]*\).*" "\\1" _eigen3_version_major "${_eigen3_version_lines}")
    string(REGEX REPLACE ".*EIGEN_MAJOR_VERSION *\([0-9]*\).*" "\\1" _eigen3_version_minor "${_eigen3_version_lines}")
    string(REGEX REPLACE ".*EIGEN_MINOR_VERSION *\([0-9]*\).*" "\\1" _eigen3_version_patch "${_eigen3_version_lines}")

    set(EIGEN3_VERSION "${_eigen3_version_major}.${_eigen3_version_minor}.${_eigen3_version_patch}")
    unset(_eigen3_version_lines)
    unset(_eigen3_version_major)
    unset(_eigen3_version_minor)
    unset(_eigen3_version_patch)
endif ()

# ----- Handle the standard arguments --------------------------------------- #
find_package_handle_standard_args(Eigen3
                                  REQUIRED_VARS EIGEN3_INCLUDE_DIR
                                  VERSION_VAR EIGEN3_VERSION)

# ----- Export the target --------------------------------------------------- #
if (Eigen3_FOUND)
    set(EIGEN3_INCLUDE_DIRS "${EIGEN3_INCLUDE_DIR}")

    if (NOT TARGET Eigen3::Eigen)
        add_library(Eigen3::Eigen INTERFACE IMPORTED)
        set_target_properties(
                Eigen3::Eigen PROPERTIES
                INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}")
    endif ()
endif ()

mark_as_advanced(EIGEN3_INCLUDE_DIR EIGEN3_VERSION)
