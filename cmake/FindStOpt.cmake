# --------------------------------------------------------------------------- #
#    CMake find module for StOpt                                              #
#                                                                             #
#    Accepts the following PATHS:                                             #
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
find_package(BZip2 REQUIRED QUIET)
find_package(ZLIB REQUIRED QUIET)
find_package(Boost REQUIRED COMPONENTS system timer QUIET)

# This will try first with Eigen3 own configuration file,
# then with the find module we provide.
find_package(Eigen3 QUIET NO_MODULE)
if (NOT Eigen3_FOUND)
    get_filename_component(FIND_MODULE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
    list(APPEND CMAKE_MODULE_PATH ${FIND_MODULE_DIR})
    find_package(Eigen3 REQUIRED)
    list(REMOVE_AT CMAKE_MODULE_PATH -1)
endif ()

# ----- Find the geners library --------------------------------------------- #
find_path(geners_INCLUDE_DIR
          NAMES geners/uriUtils.hh
          PATHS ${STOPT_INC}
          DOC "geners include directory")

find_library(geners_LIBRARY
             NAMES geners
             PATHS ${STOPT_LIB}
             DOC "geners library")

mark_as_advanced(geners_INCLUDE_DIR geners_LIBRARY)

# ----- Find the StOpt library ---------------------------------------------- #
find_path(StOpt_INCLUDE_DIR
          NAMES StOpt/sddp/OptimizerSDDPBase.h
          PATHS ${STOPT_INC}
          DOC "StOpt include directory")

find_library(StOpt_LIBRARY
             NAMES StOpt
             PATHS ${STOPT_LIB}
             DOC "StOpt library")

mark_as_advanced(StOpt_INCLUDE_DIR StOpt_LIBRARY)

# ----- Handle the standard arguments --------------------------------------- #
find_package_handle_standard_args(
        StOpt
        REQUIRED_VARS
        geners_LIBRARY geners_INCLUDE_DIR
        StOpt_LIBRARY StOpt_INCLUDE_DIR)

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
