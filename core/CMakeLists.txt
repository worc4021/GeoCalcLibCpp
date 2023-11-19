cmake_minimum_required(VERSION 3.21)
project(GeoCalcLibCore)

set(libsource)

option(USE_GMPXX "Use GMPXX" OFF)

list(APPEND libsource
            src/lrsinterface.cpp
            src/mpq.cpp
            src/mpz.cpp
        )

if (USE_GMPXX)
    list(APPEND libsource
        src/gmpMatrix.cpp
        )
else()
    list(APPEND libsource
        src/mpqMatrix.cpp
        )
endif()

add_library(GeoCalcLib ${libsource})

target_include_directories(GeoCalcLib PRIVATE ${lrslib_SOURCE_DIR})
set_target_properties(GeoCalcLib 
        PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(GeoCalcLib PUBLIC Eigen3::Eigen Boost::boost lrslib Matlab::mex)

if (USE_GMPXX)
    target_link_libraries(GeoCalcLib PUBLIC GMP::GMPXX GMP::GMP)
    target_compile_definitions(GeoCalcLib PUBLIC GMPMATRIX)
endif()

target_include_directories(GeoCalcLib 
        PUBLIC 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)

target_compile_features(GeoCalcLib INTERFACE cxx_std_17)

install(TARGETS GeoCalcLib EXPORT GeoCalcLibConfig
    ARCHIVE  DESTINATION lib
    LIBRARY  DESTINATION lib
    RUNTIME  DESTINATION bin
)

install(EXPORT GeoCalcLibConfig DESTINATION share/GeoCalcLib/cmake)
        