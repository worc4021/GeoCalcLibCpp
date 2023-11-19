cmake_minimum_required(VERSION 3.21)
project(lrslib)

set(LRS_VERSION "072" CACHE STRING "LRS version")

include(FetchContent)
FetchContent_Declare(
    lrslibsrc
    URL "https://cgm.cs.mcgill.ca/~avis/C/lrslib/archive/lrslib-${LRS_VERSION}.tar.gz"
)

FetchContent_MakeAvailable(lrslibsrc)

find_package(PkgConfig REQUIRED)
pkg_check_modules(gmp REQUIRED IMPORTED_TARGET gmp)

add_library(lrslib 
        ${lrslibsrc_SOURCE_DIR}/lrs.c 
        ${lrslibsrc_SOURCE_DIR}/lrslib.c 
        ${lrslibsrc_SOURCE_DIR}/lrsgmp.c 
        ${lrslibsrc_SOURCE_DIR}/lrsdriver.c)

set_target_properties(lrslib 
            PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_link_libraries(lrslib PUBLIC PkgConfig::gmp)
target_compile_definitions(lrslib PUBLIC GMP TIMES SIGNALS B64 LRS_QUIET NOINFO DIRECT)

target_include_directories(lrslib PUBLIC
    $<BUILD_INTERFACE:${lrslibsrc_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>)

install(DIRECTORY ${lrslibsrc_SOURCE_DIR} DESTINATION include FILES_MATCHING PATTERN "*.h")

install(TARGETS lrslib EXPORT lrslibConfig
    ARCHIVE  DESTINATION lib
    LIBRARY  DESTINATION lib
    RUNTIME  DESTINATION bin
)

install(EXPORT lrslibConfig DESTINATION share/lrslib/cmake)
