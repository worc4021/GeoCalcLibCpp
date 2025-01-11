cmake_minimum_required(VERSION 3.24)
project(lrslib)

set(LRS_VERSION "073" CACHE STRING "LRS version")
option(LRS_BUILD_EXECUTABLE "Build the lrs executable" ON)

include(FetchContent)
FetchContent_Declare(
    lrslibsrc
    URL "https://cgm.cs.mcgill.ca/~avis/C/lrslib/archive/lrslib-${LRS_VERSION}.tar.gz"
    )
FetchContent_MakeAvailable(lrslibsrc)
message(INFO "LRS libsrc is ${lrslibsrc_SOURCE_DIR}")

# Protect against libgen.h being included on Windows. We provide a dummy to obtain basename and dirname. Imagine if windows had a patch command..!
# file(READ ${lrslibsrc_SOURCE_DIR}/lrslib.c CONTENT)
# string(REPLACE "#include <libgen.h>" "#ifndef WIN32\n#include <libgen.h>\n#endif" MODIFIED_CONTENT "${CONTENT}")
# file(WRITE ${lrslibsrc_SOURCE_DIR}/lrslib.c "${MODIFIED_CONTENT}")

# file(READ ${lrslibsrc_SOURCE_DIR}/lrsgmp.h CONTENT)
# string(REPLACE "#include <unistd.h>" "#ifdef WIN32\n#include <io.h>\n#else\n#include <libgen.h>\n#endif" MODIFIED_CONTENT "${CONTENT}")
# file(WRITE ${lrslibsrc_SOURCE_DIR}/lrsgmp.h "${MODIFIED_CONTENT}")

find_package(PkgConfig REQUIRED)
pkg_check_modules(gmp REQUIRED IMPORTED_TARGET gmp)

include(FindOpenMP)

set(ARITH lrsarith-011)

# set(lrssrc)

# list(APPEND lrssrc  ${lrslibsrc_SOURCE_DIR}/lrslib.c 
#                     ${lrslibsrc_SOURCE_DIR}/lrsgmp.c 
#                     ${lrslibsrc_SOURCE_DIR}/lrsdriver.c)

if (WIN32)
set(WIN_CFLAGS __STDC__)
list(APPEND WIN_CFLAGS PATH_MAX=1024)
# list(APPEND lrssrc src/libgen-mockup.cpp)
add_library(mockup src/libgen-mockup.cpp
                   src/tempfile.cpp)
set_target_properties(mockup PROPERTIES 
                   POSITION_INDEPENDENT_CODE ON
                   MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
target_include_directories(mockup PUBLIC include)
endif()

add_library(lrs.o OBJECT ${lrslibsrc_SOURCE_DIR}/lrs.c)
set_target_properties(lrs.o PROPERTIES 
                            POSITION_INDEPENDENT_CODE ON
                            MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
target_compile_definitions(lrs.o PRIVATE MA B128 SIGNALS TIMES LRS_QUIET NOINFO DIRECT ${WIN_CFLAGS})
target_include_directories(lrs.o PRIVATE ${lrslibsrc_SOURCE_DIR}/${ARITH})

add_library(lrslong1.o OBJECT ${lrslibsrc_SOURCE_DIR}/${ARITH}/lrslong.c)
set_target_properties(lrslong1.o PROPERTIES 
                            POSITION_INDEPENDENT_CODE ON
                            MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
target_compile_definitions(lrslong1.o PRIVATE MA SAFE LRSLONG SIGNALS TIMES LRS_QUIET NOINFO DIRECT ${WIN_CFLAGS})
target_include_directories(lrslong1.o PRIVATE ${lrslibsrc_SOURCE_DIR}/${ARITH})

add_library(lrslong2.o OBJECT ${lrslibsrc_SOURCE_DIR}/${ARITH}/lrslong.c)
set_target_properties(lrslong2.o PROPERTIES 
                            POSITION_INDEPENDENT_CODE ON
                            MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
target_compile_definitions(lrslong2.o PRIVATE MA SAFE LRSLONG B128 SIGNALS TIMES LRS_QUIET NOINFO DIRECT ${WIN_CFLAGS})
target_include_directories(lrslong2.o PRIVATE ${lrslibsrc_SOURCE_DIR}/${ARITH})


add_library(lrslib1.o OBJECT ${lrslibsrc_SOURCE_DIR}/lrslib.c)
set_target_properties(lrslib1.o PROPERTIES 
                            POSITION_INDEPENDENT_CODE ON
                            MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
target_compile_definitions(lrslib1.o PRIVATE MA SAFE LRSLONG PLRS SIGNALS TIMES LRS_QUIET NOINFO DIRECT ${WIN_CFLAGS})
target_include_directories(lrslib1.o PRIVATE ${lrslibsrc_SOURCE_DIR}/${ARITH})
target_link_libraries(lrslib1.o PRIVATE OpenMP::OpenMP_C)


add_library(lrslib2.o OBJECT ${lrslibsrc_SOURCE_DIR}/lrslib.c)
set_target_properties(lrslib2.o PROPERTIES 
                            POSITION_INDEPENDENT_CODE ON
                            MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
target_compile_definitions(lrslib2.o PRIVATE MA SAFE LRSLONG B128 PLRS SIGNALS TIMES LRS_QUIET NOINFO DIRECT ${WIN_CFLAGS})
target_include_directories(lrslib2.o PRIVATE ${lrslibsrc_SOURCE_DIR}/${ARITH})
target_link_libraries(lrslib2.o PRIVATE OpenMP::OpenMP_C)

add_library(lrslibgmp.o OBJECT ${lrslibsrc_SOURCE_DIR}/lrslib.c)
set_target_properties(lrslibgmp.o PROPERTIES 
                            POSITION_INDEPENDENT_CODE ON
                            MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
target_compile_definitions(lrslibgmp.o PRIVATE MA GMP PLRS SIGNALS TIMES LRS_QUIET NOINFO DIRECT ${WIN_CFLAGS})
target_include_directories(lrslibgmp.o PRIVATE ${lrslibsrc_SOURCE_DIR}/${ARITH})
target_link_libraries(lrslibgmp.o PUBLIC PkgConfig::gmp OpenMP::OpenMP_C)

add_library(lrsgmp.o OBJECT ${lrslibsrc_SOURCE_DIR}/${ARITH}/lrsgmp.c)
set_target_properties(lrsgmp.o PROPERTIES 
                            POSITION_INDEPENDENT_CODE ON
                            MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
target_compile_definitions(lrsgmp.o PRIVATE MA GMP SIGNALS TIMES LRS_QUIET NOINFO DIRECT ${WIN_CFLAGS})
target_include_directories(lrsgmp.o PRIVATE ${lrslibsrc_SOURCE_DIR} ${lrslibsrc_SOURCE_DIR}/${ARITH})
target_link_libraries(lrsgmp.o PUBLIC PkgConfig::gmp)

# lrs.o lrslong1.o lrslong2.o lrslib1.o lrslib2.o lrslibgmp.o lrsgmp.o lrsdriver.o

if (WIN32)
    target_link_libraries(lrslib1.o PRIVATE mockup)
    target_link_libraries(lrslib2.o PRIVATE mockup)
    target_link_libraries(lrslibgmp.o PRIVATE mockup)
endif()

add_executable(lrs ${lrslibsrc_SOURCE_DIR}/lrsdriver.c)
set_target_properties(lrs   PROPERTIES 
                            POSITION_INDEPENDENT_CODE ON
                            MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
target_include_directories(lrs PRIVATE ${lrslibsrc_SOURCE_DIR} ${lrslibsrc_SOURCE_DIR}/${ARITH})
target_link_libraries(lrs PRIVATE lrs.o lrslong1.o lrslong2.o lrslib1.o lrslib2.o lrslibgmp.o lrsgmp.o OpenMP::OpenMP_C PkgConfig::gmp)
target_compile_definitions(lrs PRIVATE ${WIN_CFLAGS})

# add_library(lrslib 
#         ${lrssrc})


# set_target_properties(lrslib 
#             PROPERTIES 
#             POSITION_INDEPENDENT_CODE ON
#             MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")

# target_link_libraries(lrslib PUBLIC PkgConfig::gmp)
# target_compile_definitions(lrslib PUBLIC GMP TIMES SIGNALS B64 LRS_QUIET NOINFO DIRECT)

# if (WIN32)
# target_compile_definitions(lrslib PUBLIC PATH_MAX=1024 __STDC__ SIGNALS)
# endif()

# if (LRS_BUILD_EXECUTABLE)
#     add_executable(lrsgmp  
#         ${lrslibsrc_SOURCE_DIR}/lrs.c)
        
#     target_link_libraries(lrsgmp PUBLIC lrslib)
#     install(TARGETS lrsgmp DESTINATION bin)
#     set_target_properties(lrsgmp 
#                         PROPERTIES 
#                         MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
# endif()

# target_include_directories(lrslib PUBLIC
#     $<BUILD_INTERFACE:${lrslibsrc_SOURCE_DIR}>
#     $<INSTALL_INTERFACE:include>)

# install(DIRECTORY ${lrslibsrc_SOURCE_DIR} DESTINATION include FILES_MATCHING PATTERN "*.h")

# install(TARGETS lrslib EXPORT lrslibConfig
#     ARCHIVE  DESTINATION lib
#     LIBRARY  DESTINATION lib
#     RUNTIME  DESTINATION bin
# )

# install(EXPORT lrslibConfig DESTINATION share/lrslib/cmake)
