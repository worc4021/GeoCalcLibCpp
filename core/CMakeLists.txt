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

add_library(${PROJECT_NAME} ${libsource})
add_library(GeoCalcLib::Core ALIAS ${PROJECT_NAME})

target_include_directories(${PROJECT_NAME} PRIVATE ${lrslib_SOURCE_DIR})
set_target_properties(${PROJECT_NAME}
                    PROPERTIES 
                    POSITION_INDEPENDENT_CODE ON
                    MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")

target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen Boost::boost lrslib Matlab::mex)

if (USE_GMPXX)
    target_link_libraries(${PROJECT_NAME} PUBLIC PkgConfig::gmpxx)
    target_compile_definitions(${PROJECT_NAME} PUBLIC GMPMATRIX)
endif()

target_include_directories(${PROJECT_NAME} 
        PUBLIC 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)

install(TARGETS ${PROJECT_NAME} EXPORT GeoCalcLibConfig
    ARCHIVE  DESTINATION lib
    LIBRARY  DESTINATION lib
    RUNTIME  DESTINATION bin
)

install(EXPORT GeoCalcLibConfig DESTINATION share/GeoCalcLib/cmake)
        