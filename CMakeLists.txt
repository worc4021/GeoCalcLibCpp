set(CMAKE_TOOLCHAIN_FILE ${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake CACHE STRING "Vcpkg toolchain file")
cmake_minimum_required(VERSION 3.21)
project(GeoCalcLib)

set(CMAKE_CXX_STANDARD 20)

find_package(Matlab REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)
find_package(Boost REQUIRED)

add_subdirectory(lrslib)

option(BUILD_TESTS "Build tests" ON)
add_subdirectory(MexUtilities)
add_subdirectory(core)



matlab_add_mex(
    NAME facetEnumeration
    SRC src/facetEnumeration.cpp
    LINK_TO GeoCalcLib::Core MexUtilities::MexUtilities
    R2018a
)

set_target_properties(facetEnumeration 
                        PROPERTIES 
                        MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")

matlab_add_mex(
    NAME vertexEnumeration
    SRC src/vertexEnumeration.cpp
    LINK_TO GeoCalcLib::Core MexUtilities::MexUtilities
    R2018a
)

set_target_properties(vertexEnumeration 
                        PROPERTIES 
                        MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")

if (BUILD_TESTS)
    include(CTest)
    enable_testing()
    
    add_executable(GeoCalcLibTest src/gmpcpptest.cpp)
    target_link_libraries(GeoCalcLibTest PRIVATE GeoCalcLib::Core Eigen3::Eigen)
    set_target_properties(GeoCalcLibTest
                        PROPERTIES 
                        MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")

    matlab_add_unit_test(
        NAME mexTest
        UNITTEST_FILE ${CMAKE_CURRENT_SOURCE_DIR}/GeoCalcLibTest.m
        ADDITIONAL_PATH $<TARGET_FILE_DIR:vertexEnumeration>
    )
endif()