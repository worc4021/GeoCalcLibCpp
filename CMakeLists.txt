set(CMAKE_TOOLCHAIN_FILE ${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake CACHE STRING "Vcpkg toolchain file")
cmake_minimum_required(VERSION 3.21)
project(GeoCalcLib)

include(CTest)

find_package(Matlab REQUIRED COMPONENTS MAIN_PROGRAM)
find_package(Eigen3 CONFIG REQUIRED)
find_package(Boost REQUIRED)

add_subdirectory(lrslib)
add_subdirectory(MexUtilities)
add_subdirectory(core)



matlab_add_mex(
    NAME facetEnumeration
    SRC src/facetEnumeration.cpp
    LINK_TO GeoCalcLib::Core MexUtilities::MexUtilities
    R2018a
)

matlab_add_mex(
    NAME vertexEnumeration
    SRC src/vertexEnumeration.cpp
    LINK_TO GeoCalcLib::Core MexUtilities::MexUtilities
    R2018a
)

matlab_add_unit_test(
    NAME mexTest
    UNITTEST_FILE ${CMAKE_CURRENT_SOURCE_DIR}/GeoCalcLibTest.m
)