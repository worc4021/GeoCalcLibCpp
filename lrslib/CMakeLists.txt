cmake_minimum_required(VERSION 3.24)
project(lrslib)

set(LRS_VERSION "072" CACHE STRING "LRS version")
option(LRS_BUILD_EXECUTABLE "Build the lrs executable" ON)

include(FetchContent)
FetchContent_Declare(
    lrslibsrc
    URL "https://cgm.cs.mcgill.ca/~avis/C/lrslib/archive/lrslib-${LRS_VERSION}.tar.gz"
    DOWNLOAD_EXTRACT_TIMESTAMP true
    )
FetchContent_MakeAvailable(lrslibsrc)

# Protect against libgen.h being included on Windows. We provide a dummy to obtain basename and dirname. Imagine if windows had a patch command..!
file(READ ${lrslibsrc_SOURCE_DIR}/lrslib.c CONTENT)
string(REPLACE "#include <libgen.h>" "#ifdef WIN32\n#include <io.h>\n#include \"mockup-header.h\"\n#else\n#include <libgen.h>\n#endif" MODIFIED_CONTENT "${CONTENT}")
file(WRITE ${lrslibsrc_SOURCE_DIR}/lrslib.c "${MODIFIED_CONTENT}")

file(READ ${lrslibsrc_SOURCE_DIR}/lrsgmp.h CONTENT)
string(REPLACE "#include <unistd.h>" "#ifdef WIN32\n#include <io.h>\n#include \"mockup-header.h\"\n#else\n#include <libgen.h>\n#endif" MODIFIED_CONTENT "${CONTENT}")
file(WRITE ${lrslibsrc_SOURCE_DIR}/lrsgmp.h "${MODIFIED_CONTENT}")

find_package(PkgConfig REQUIRED)
pkg_check_modules(gmp REQUIRED IMPORTED_TARGET gmp)

add_library(lrslong1 ${lrslibsrc_SOURCE_DIR}/lrslong.c $<$<BOOL:WIN32>:src/libgen-mockup.cpp>)
target_compile_definitions(lrslong1 PRIVATE MA SAFE LRSLONG OMIT_SIGNALS TIMES)
add_library(lrslong2 ${lrslibsrc_SOURCE_DIR}/lrslong.c $<$<BOOL:WIN32>:src/libgen-mockup.cpp>)
target_compile_definitions(lrslong2 PRIVATE B128 MA SAFE LRSLONG OMIT_SIGNALS TIMES)
add_library(lrslib1 ${lrslibsrc_SOURCE_DIR}/lrslib.c $<$<BOOL:WIN32>:src/libgen-mockup.cpp>)
target_compile_definitions(lrslib1 PRIVATE MA SAFE LRSLONG OMIT_SIGNALS SIGNALS TIMES)
add_library(lrslib2 ${lrslibsrc_SOURCE_DIR}/lrslib.c $<$<BOOL:WIN32>:src/libgen-mockup.cpp>)
target_compile_definitions(lrslib2 PRIVATE B128 MA SAFE LRSLONG OMIT_SIGNALS SIGNALS TIMES)
add_library(lrslibgmp ${lrslibsrc_SOURCE_DIR}/lrslib.c $<$<BOOL:WIN32>:src/libgen-mockup.cpp>)
target_compile_definitions(lrslibgmp PRIVATE MA GMP OMIT_SIGNALS SIGNALS TIMES)
target_link_libraries(lrslibgmp PUBLIC PkgConfig::gmp)
add_library(lrslibmp ${lrslibsrc_SOURCE_DIR}/lrslib.c $<$<BOOL:WIN32>:src/libgen-mockup.cpp>)
target_compile_definitions(lrslibmp PRIVATE MA MP OMIT_SIGNALS TIMES)
add_library(lrs1 ${lrslibsrc_SOURCE_DIR}/lrs.c $<$<BOOL:WIN32>:src/libgen-mockup.cpp>)
target_compile_definitions(lrs1 PRIVATE B128 MA OMIT_SIGNALS SIGNALS TIMES)


set(lrssrc)
list(APPEND lrssrc ${lrslibsrc_SOURCE_DIR}/lrsdriver.c)

if (WIN32)
list(APPEND lrssrc src/libgen-mockup.cpp)
endif()

add_library(lrslib 
        ${lrssrc})


add_executable(lrs ${lrslibsrc_SOURCE_DIR}/lrsdriver.c $<$<BOOL:WIN32>:src/libgen-mockup.cpp>)
target_link_libraries(lrs PRIVATE lrslong1 lrslong2 lrslib1 lrslib2 lrslibgmp lrs1)
target_include_directories(lrs PRIVATE ${lrslibsrc_SOURCE_DIR})
target_compile_definitions(lrs PUBLIC B128 MA OMIT_SIGNALS SIGNALS TIMES)

target_link_libraries(lrslib PUBLIC PkgConfig::gmp)
target_compile_definitions(lrslib PUBLIC GMP TIMES SIGNALS B64 LRS_QUIET NOINFO DIRECT)



if (LRS_BUILD_EXECUTABLE)
    add_executable(lrsgmp  
        ${lrslibsrc_SOURCE_DIR}/lrs.c)
        
    target_link_libraries(lrsgmp PUBLIC lrslib)
    install(TARGETS lrsgmp DESTINATION bin)
    
endif()

get_property(current_targets DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY BUILDSYSTEM_TARGETS)
foreach(TARGET ${current_targets})
    set_target_properties(${TARGET}
            PROPERTIES 
            POSITION_INDEPENDENT_CODE ON
            MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>"
            WINDOWS_EXPORT_ALL_SYMBOLS ON)
    if (WIN32)
        target_compile_definitions(${TARGET} PUBLIC PATH_MAX=1024 __STDC__)
        target_include_directories(${TARGET} AFTER PRIVATE include)
    endif()
endforeach()


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
