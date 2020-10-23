# --------------------------------------------------------------------------- #
#    CMake find module for StOpt                                              #
#                                                                             #
#    Accepts the following HINTS:                                             #
#                                                                             #
#    - STOPT_INC - Custom path to StOpt headers                               #
#    - STOPT_LIB - Custom path to StOpt libraries                             #
#                                                                             #
#    Provides the following imported targets:                                 #
#                                                                             #
#    - StOpt::StOpt - the StOpt library                                       #
#    - StOpt::geners - the geners library                                     #
#                                                                             #
#                              Niccolo' Iardella                              #
#                          Operations Research Group                          #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
include(FindPackageHandleStandardArgs)

# ----- Requirements -------------------------------------------------------- #
find_package(Eigen3 REQUIRED QUIET)
find_package(BZip2 REQUIRED QUIET)
find_package(ZLIB REQUIRED QUIET)
find_package(Boost REQUIRED COMPONENTS system timer QUIET)

# ----- Find the geners library --------------------------------------------- #
find_path(geners_INCLUDE_DIR
          NAMES geners/uriUtils.hh
          HINTS ${STOPT_INC}
          DOC "geners include directories")

find_library(geners_LIBRARY
             NAMES geners
             HINTS ${STOPT_LIB}
             DOC "geners library")

mark_as_advanced(geners_INCLUDE_DIR geners_LIBRARY)

# ----- Find the StOpt library ---------------------------------------------- #
find_path(StOpt_INCLUDE_DIR
          NAMES StOpt/sddp/OptimizerSDDPBase.h
          HINTS ${STOPT_INC}
          DOC "StOpt include directories")

find_library(StOpt_LIBRARY
             NAMES StOpt
             HINTS ${STOPT_LIB}
             DOC "StOpt library")

mark_as_advanced( StOpt_INCLUDE_DIR StOpt_LIBRARY)

# ----- Handle the standard arguments --------------------------------------- #

find_package_handle_standard_args(
        StOpt
        REQUIRED_VARS
        geners_LIBRARY geners_INCLUDE_DIR StOpt_LIBRARY StOpt_INCLUDE_DIR)

# ----- Export the target(s) ------------------------------------------------ #
if (StOpt_FOUND)
    if (NOT TARGET StOpt::geners)
        set(geners_INCLUDE_DIRS "${geners_INCLUDE_DIR}")
        set(geners_LIBRARIES "${geners_LIBRARY}")

        add_library(StOpt::geners UNKNOWN IMPORTED)
        set_target_properties(
                StOpt::geners PROPERTIES
                IMPORTED_LOCATION "${geners_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${geners_INCLUDE_DIR}")
    endif ()

    if (NOT TARGET StOpt::StOpt)
        set(StOpt_INCLUDE_DIRS "${StOpt_INCLUDE_DIR}")
        set(StOpt_LIBRARIES "${StOpt_LIBRARY}")

        add_library(StOpt::StOpt UNKNOWN IMPORTED)
        set_target_properties(
                StOpt::StOpt PROPERTIES
                IMPORTED_LOCATION "${StOpt_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${StOpt_INCLUDE_DIR}"
                INTERFACE_LINK_LIBRARIES "StOpt::geners;Eigen3::Eigen;BZip2::BZip2;ZLIB::ZLIB;Boost::system;Boost::timer")
    endif ()
endif ()
