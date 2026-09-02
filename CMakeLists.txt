#-------------------------------------------------------------------------------
# morpho-pairwise/CMakeLists.txt
#-------------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.13)

# Single-config generators (Unix Makefiles, Ninja) default to an empty
# CMAKE_BUILD_TYPE, which means no -O3. Multi-config generators (Xcode, VS)
# ignore CMAKE_BUILD_TYPE; --config selects the variant at build time.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "CMAKE_BUILD_TYPE unset; defaulting to Release")
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Build type (Debug, Release, RelWithDebInfo, MinSizeRel)" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
endif()

project(morpho-pairwise)

# Build the library as a plugin
add_library(pairwise MODULE "")

# Suppress 'lib' prefix
set_target_properties(pairwise PROPERTIES PREFIX "")

#-------------------------------------------------------------------------------
# Release LTO
#-------------------------------------------------------------------------------

include(CheckIPOSupported)
check_ipo_supported(RESULT PAIRWISE_IPO_SUPPORTED OUTPUT PAIRWISE_IPO_OUTPUT LANGUAGES C)
if(PAIRWISE_IPO_SUPPORTED)
    set_property(TARGET pairwise PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        set_property(TARGET pairwise PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif()
    message(STATUS "IPO/LTO enabled for Release")
else()
    message(STATUS "IPO/LTO not enabled (${PAIRWISE_IPO_OUTPUT})")
endif()

# Add sources
add_subdirectory(src)

find_path(MORPHO_INCLUDE_DIR
    NAMES morpho.h morpho/morpho.h
    HINTS
        $ENV{MORPHO_INCLUDE_DIR}
        /usr/local/opt/morpho/include
        /opt/homebrew/opt/morpho/include
        /home/linuxbrew/.linuxbrew/include/morpho
        /home/linuxbrew/.linuxbrew/include
        /usr/local/include/morpho
        /usr/local/include
        /opt/homebrew/include
)

if(NOT MORPHO_INCLUDE_DIR)
    message(FATAL_ERROR
        "Could not find Morpho headers (morpho.h). "
        "Install Morpho or set MORPHO_INCLUDE_DIR to the include directory "
        "that contains morpho.h or morpho/morpho.h.")
endif()

get_filename_component(MORPHO_INCLUDE_PARENT ${MORPHO_INCLUDE_DIR} DIRECTORY)

# Prefer a Morpho library that matches the headers we found
# (.../include/morpho -> .../lib, or .../include -> .../lib).
get_filename_component(MORPHO_PREFIX ${MORPHO_INCLUDE_PARENT} DIRECTORY)
set(MORPHO_LIB_HINTS
    "${MORPHO_INCLUDE_PARENT}/../lib"
    "${MORPHO_PREFIX}/lib"
)

find_library(MORPHO_LIBRARY
    NAMES morpho libmorpho
    HINTS
        $ENV{MORPHO_LIBRARY}
        ${MORPHO_LIB_HINTS}
        /usr/local/opt/morpho/lib
        /opt/homebrew/opt/morpho/lib
        /home/linuxbrew/.linuxbrew/lib
        /usr/local/lib
        /opt/homebrew/lib
)

find_library(MATH_LIBRARY m)

if(NOT MORPHO_LIBRARY)
    message(FATAL_ERROR
        "Could not find libmorpho. "
        "Install Morpho or set MORPHO_LIBRARY to the full path of the Morpho library.")
endif()

# Include path should allow #include <morpho/morpho.h> and related headers.
target_include_directories(pairwise PUBLIC
    ${MORPHO_INCLUDE_DIR}
    ${MORPHO_INCLUDE_PARENT}
)

file(GLOB morpho_subdirectories LIST_DIRECTORIES true ${MORPHO_INCLUDE_DIR}/*)
foreach(dir ${morpho_subdirectories})
    if(IS_DIRECTORY ${dir})
        target_include_directories(pairwise PUBLIC ${dir})
    endif()
endforeach()

target_link_libraries(pairwise ${MORPHO_LIBRARY} ${MATH_LIBRARY})

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}" CACHE PATH "Install prefix" FORCE)
endif()

install(TARGETS pairwise LIBRARY DESTINATION lib/)
