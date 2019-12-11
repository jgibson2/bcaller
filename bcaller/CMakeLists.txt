cmake_minimum_required(VERSION 3.14)
include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)
project(bcaller)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

find_package(HighFive 2.1 REQUIRED)

include_directories(${CMAKE_SOURCE_DIR}/deps/vcflib/include)
include_directories(${CMAKE_SOURCE_DIR}/deps/vcflib/tabixpp/htslib)

add_library(vcf_htslib STATIC IMPORTED)
set_property(TARGET vcf_htslib PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/deps/vcflib/tabixpp/htslib/libhts.a)

add_library(vcflib STATIC IMPORTED)
set_property(TARGET vcflib PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/deps/vcflib/lib/libvcflib.a)

add_executable(vcf_to_hdf5 vcf_to_hdf5.cpp)
target_link_libraries(vcf_to_hdf5 pthread HighFive hdf5_cpp vcflib vcf_htslib z m lzma bz2)

add_executable(hdf5_to_window_frequencies hdf5_to_window_frequencies.cpp)
target_link_libraries(hdf5_to_window_frequencies pthread HighFive hdf5_cpp vcflib vcf_htslib z m lzma bz2)

