cmake_minimum_required(VERSION 3.24)
project(lrslib)

set(LRS_VERSION "072" CACHE STRING "LRS version")

include(FetchContent)
FetchContent_Declare(
    lrslibsrc
    URL "https://cgm.cs.mcgill.ca/~avis/C/lrslib/archive/lrslib-${LRS_VERSION}.tar.gz"
    DOWNLOAD_EXTRACT_TIMESTAMP true
    )
FetchContent_MakeAvailable(lrslibsrc)

# Protect against libgen.h being included on Windows. We provide a dummy to obtain basename and dirname. Imagine if windows had a patch command..!
file(READ ${lrslibsrc_SOURCE_DIR}/lrslib.c CONTENT)
string(REPLACE "#include <libgen.h>" "#ifndef WIN32\n#include <libgen.h>\n#endif" MODIFIED_CONTENT "${CONTENT}")
file(WRITE ${lrslibsrc_SOURCE_DIR}/lrslib.c "${MODIFIED_CONTENT}")

file(READ ${lrslibsrc_SOURCE_DIR}/lrsgmp.h CONTENT)
string(REPLACE "#include <unistd.h>" "#ifdef WIN32\n#include <io.h>\n#else\n#include <libgen.h>\n#endif" MODIFIED_CONTENT "${CONTENT}")
file(WRITE ${lrslibsrc_SOURCE_DIR}/lrsgmp.h "${MODIFIED_CONTENT}")

find_package(PkgConfig REQUIRED)
pkg_check_modules(gmp REQUIRED IMPORTED_TARGET gmp)

set(lrssrc)

list(APPEND lrssrc  ${lrslibsrc_SOURCE_DIR}/lrslib.c 
                    ${lrslibsrc_SOURCE_DIR}/lrsgmp.c 
                    ${lrslibsrc_SOURCE_DIR}/lrsdriver.c)

if (WIN32)
list(APPEND lrssrc src/libgen-mockup.cpp)
endif()

add_library(lrslib 
        ${lrssrc})

set_target_properties(lrslib 
            PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_link_libraries(lrslib PUBLIC PkgConfig::gmp)
target_compile_definitions(lrslib PUBLIC GMP TIMES SIGNALS B64 LRS_QUIET NOINFO DIRECT)

if (WIN32)
target_compile_definitions(lrslib PRIVATE PATH_MAX=1024 )
target_compile_definitions(lrslib PUBLIC __STDC__)

endif()

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
