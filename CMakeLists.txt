# --------------------------------------------------------------------------- #
#    Main CMake file for SMS++ SDDPBlock                                      #
#                                                                             #
#                              Niccolo' Iardella                              #
#                          Operations Research Group                          #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
cmake_minimum_required(VERSION 3.11)
project(SDDPBlock VERSION 0.3.0 LANGUAGES CXX)

# ----- Requirements -------------------------------------------------------- #
# This is to use the find modules we provide
find_package(SMSpp REQUIRED)
find_package(StochasticBlock REQUIRED)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package(StOpt REQUIRED)
list(REMOVE_AT CMAKE_MODULE_PATH -1)

# ----- Library ------------------------------------------------------------- #
add_library(SDDPBlock)
target_compile_features(SDDPBlock PUBLIC cxx_std_17)

target_sources(SDDPBlock PRIVATE
               src/SDDPBlock.cpp
               src/SDDPGreedySolver.cpp)

# Different INSTALL_INTERFACE and BUILD_INTERFACE paths are used when
# generating the target import file (***Targets.cmake), that must be different
# for local build and installation.
target_include_directories(
        SDDPBlock PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include/SMS++>)

target_link_libraries(SDDPBlock PUBLIC
                      SMS++::SMSpp
                      SMS++::StochasticBlock
                      StOpt::StOpt)

# Alias for testing
add_library(SMS++::SDDPBlock ALIAS SDDPBlock)

# ----- Install instructions ------------------------------------------------ #
# For installing the library and all its CMake modules on the system.
# Not required for local builds (see below).

include(GNUInstallDirs)

# Install library
install(TARGETS SDDPBlock
        EXPORT SDDPBlock-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Install headers
install(DIRECTORY include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SMS++)

# Install target import file
install(EXPORT SDDPBlock-targets
        FILE SDDPBlockTargets.cmake
        NAMESPACE SMS++::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SDDPBlock)

# Generate package version file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/SDDPBlockConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion)

# Generate package configuration file
configure_package_config_file(
        ${CMAKE_CURRENT_LIST_DIR}/cmake/SDDPBlockConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/SDDPBlockConfig.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SDDPBlock)

# Install package version and configuration files
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/SDDPBlockConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/SDDPBlockConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SDDPBlock)

# Shipping custom find modules should not be encouraged because third
# party developers should provide them.
# We ship FindStOpt because StOpt does not provide its own,
# and FindEigen3 because the user may want not to install Eigen3.
install(FILES
        ${CMAKE_CURRENT_LIST_DIR}/cmake/FindStOpt.cmake
        ${CMAKE_CURRENT_LIST_DIR}/cmake/FindEigen3.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SDDPBlock)

# ----- Export from the build tree ------------------------------------------ #
# This is needed so other projects can find the library in a local build.
# It does not affect installation.

if (CMAKE_EXPORT_PACKAGE_REGISTRY)

    export(EXPORT SDDPBlock-targets
           FILE ${CMAKE_CURRENT_BINARY_DIR}/SDDPBlockTargets.cmake
           NAMESPACE SMS++::)

    # Export also the find modules for the local config file to use it
    configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/FindStOpt.cmake
                   ${CMAKE_CURRENT_BINARY_DIR}/FindStOpt.cmake
                   COPYONLY)

    configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/FindEigen3.cmake
                   ${CMAKE_CURRENT_BINARY_DIR}/FindEigen3.cmake
                   COPYONLY)

    # Register package in user's package registry (~/.cmake/registry)
    export(PACKAGE SDDPBlock)
endif ()

if (BUILD_TESTING)
    # Just consuming BUILD_TESTING to avoid warnings
endif ()
