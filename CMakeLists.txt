
cmake_minimum_required(VERSION 3.26)

project(morpho-pairwise)

include(CMakePrintHelpers)

add_library(pairwise MODULE src/pairwise.c)

# Find morpho header files 
target_include_directories(pairwise PUBLIC /usr/local/include)

# Ensure we build a shared library
set(BUILD_SHARED_LIBS ON)
set_property(TARGET pairwise PROPERTY POSITION_INDEPENDENT_CODE ON)

# Ensure the outfile filename doens't begin with lib
set_target_properties(pairwise PROPERTIES PREFIX "")

# Enable chained fixups - not platform indep
target_link_libraries(pairwise "-undefined dynamic_lookup")
