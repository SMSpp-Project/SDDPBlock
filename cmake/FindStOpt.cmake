# --------------------------------------------------------------------------- #
#    CMake find module for StOpt                                              #
#                                                                             #
#    This module finds StOpt include directories and libraries.               #
#    Use it by invoking find_package() with the form:                         #
#                                                                             #
#        find_package(StOpt [version] [EXACT] [REQUIRED])                     #
#                                                                             #
#    The results are stored in the following variables:                       #
#                                                                             #
#        <LibraryName>_FOUND         - True if headers are found              #
#        <LibraryName>_INCLUDE_DIRS  - Include directories                    #
#        <LibraryName>_LIBRARIES     - Libraries to be linked                 #
#        <LibraryName>_VERSION       - Version number                         #
#                                                                             #
#    The search results are saved in these persistent cache entries:          #
#                                                                             #
#        <LibraryName>_INCLUDE_DIR   - Directory containing headers           #
#        <LibraryName>_LIBRARY       - The found library                      #
#                                                                             #
#    for each <LibraryName> in: StOpt, StOpt_geners.                          #
#                                                                             #
#    This module can read a search path from the variable:                    #
#                                                                             #
#        StOpt_ROOT          - Preferred StOpt location                       #
#                                                                             #
#    The following IMPORTED targets are also defined:                         #
#                                                                             #
#        StOpt::StOpt                                                         #
#        StOpt::geners                                                        #
#                                                                             #
#    This find module is provided because StOpt does not provide              #
#    a CMake configuration file on its own.                                   #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
include(FindPackageHandleStandardArgs)

# ----- Requirements -------------------------------------------------------- #
find_package(BZip2 REQUIRED QUIET)
find_package(ZLIB REQUIRED QUIET)
find_package(Boost REQUIRED COMPONENTS system timer QUIET)

# This will try first with Eigen3 own configuration file,
# then with the find module we provide.
find_package(Eigen3 QUIET CONFIG)
if (NOT Eigen3_FOUND)
    find_package(Eigen3 REQUIRED)
endif ()

# Check if already in cache
if (StOpt_geners_INCLUDE_DIR AND StOpt_geners_LIBRARY AND StOpt_geners_LIBRARY_DEBUG AND
    StOpt_INCLUDE_DIR AND StOpt_LIBRARY AND StOpt_LIBRARY_DEBUG)
    set(StOpt_FOUND TRUE)
else ()

    # ----- Find the geners include directory ------------------------------- #
    find_path(StOpt_geners_INCLUDE_DIR
              NAMES geners
              HINTS ${StOpt_ROOT}/geners
              DOC "geners include directory.")

    # ----- Find the geners library ----------------------------------------- #
    find_library(StOpt_geners_LIBRARY
            NAMES geners
            HINTS ${StOpt_ROOT}/lib
            DOC "geners library.")

    if (UNIX)
        set(StOpt_geners_LIBRARY_DEBUG ${StOpt_geners_LIBRARY})
    else ()
        find_library(StOpt_geners_LIBRARY_DEBUG
                     NAMES geners
                     HINTS ${StOpt_ROOT}/debug/lib
                     NO_DEFAULT_PATH
                     DOC "geners debug library.")
    endif ()

    # ----- Find the StOpt include directory -------------------------------- #
    find_path(StOpt_INCLUDE_DIR
              NAMES StOpt/sddp
              HINTS ${StOpt_ROOT}
              PATH_SUFFIXES StOpt
              DOC "StOpt include directory.")

    # ----- Find the StOpt library ------------------------------------------ #
    find_library(StOpt_LIBRARY
                 NAMES StOpt
                 HINTS ${StOpt_ROOT}/lib
                 DOC "StOpt library.")

    if (UNIX)
        set(StOpt_LIBRARY_DEBUG ${StOpt_LIBRARY})
    elseif (WIN32)
        find_library(StOpt_LIBRARY_DEBUG
                     NAMES StOpt
                     HINTS ${StOpt_ROOT}/debug/lib
                     NO_DEFAULT_PATH
                     DOC "StOpt debug library.")
    endif ()

    # ----- Parse the version ----------------------------------------------- #
    if (StOpt_INCLUDE_DIR)
        file(STRINGS
                "${StOpt_INCLUDE_DIR}/StOpt/core/utils/version.h"
                _stopt_version_lines REGEX "#define STOPT_VERSION")

        string(REGEX REPLACE ".*STOPT_VERSION *\"([0-9].[0-9]*\).*" "\\1"
                _stopt_version "${_stopt_version_lines}")

        set(StOpt_VERSION "${_stopt_version}")
        unset(_stopt_version_lines)
        unset(_stopt_version)
    endif ()

    # ----- Handle the standard arguments ----------------------------------- #
    # The following macro manages the QUIET, REQUIRED and version-related
    # options passed to find_package(). It also sets <PackageName>_FOUND if
    # REQUIRED_VARS are set.
    # REQUIRED_VARS should be cache entries and not output variables. See:
    # https://cmake.org/cmake/help/latest/module/FindPackageHandleStandardArgs.html
    find_package_handle_standard_args(
            StOpt
            REQUIRED_VARS StOpt_geners_LIBRARY StOpt_geners_INCLUDE_DIR
                          StOpt_LIBRARY        StOpt_INCLUDE_DIR)
endif ()

# ----- Export the targets -------------------------------------------------- #
if (StOpt_FOUND)
    set(StOpt_geners_INCLUDE_DIRS "${StOpt_geners_INCLUDE_DIR}")
    set(StOpt_geners_LIBRARIES "${StOpt_geners_LIBRARY}" "${StOpt_geners_LIBRARY_DEBUG}")

    if (NOT TARGET StOpt::geners)
        add_library(StOpt::geners UNKNOWN IMPORTED)
        set_target_properties(
                StOpt::geners PROPERTIES
                IMPORTED_LOCATION "${StOpt_geners_LIBRARY}"
                IMPORTED_LOCATION_DEBUG "${StOpt_geners_LIBRARY_DEBUG}"
                INTERFACE_INCLUDE_DIRECTORIES "${StOpt_geners_INCLUDE_DIRS}")
    endif ()

    set(StOpt_INCLUDE_DIRS "${StOpt_INCLUDE_DIR}")
    set(StOpt_LIBRARIES "${StOpt_LIBRARY}" "${StOpt_LIBRARY_DEBUG}")

    if (NOT TARGET StOpt::StOpt)
        add_library(StOpt::StOpt UNKNOWN IMPORTED)
        set_target_properties(
                StOpt::StOpt PROPERTIES
                IMPORTED_LOCATION "${StOpt_LIBRARY}"
                IMPORTED_LOCATION_DEBUG "${StOpt_LIBRARY_DEBUG}"
                INTERFACE_INCLUDE_DIRECTORIES "${StOpt_INCLUDE_DIRS}"
                INTERFACE_LINK_LIBRARIES "StOpt::geners;Eigen3::Eigen;BZip2::BZip2;ZLIB::ZLIB;Boost::system;Boost::timer")
    endif ()
endif ()

# Variables marked as advanced are not displayed in CMake GUIs, see:
# https://cmake.org/cmake/help/latest/command/mark_as_advanced.html
mark_as_advanced(StOpt_INCLUDE_DIR      StOpt_geners_INCLUDE_DIR
                 StOpt_LIBRARY          StOpt_geners_LIBRARY
                 StOpt_LIBRARY_DEBUG    StOpt_geners_LIBRARY_DEBUG
                 StOpt_VERSION)

# --------------------------------------------------------------------------- #
