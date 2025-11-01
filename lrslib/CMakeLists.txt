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

set(LRS_ADDITIONAL_DEFINES)

include(CheckSourceCompiles)
check_source_compiles(
    C
    "#include <signal.h>
    #include <unistd.h>
    int main() { 
        signal(0, 0); 
        return 0; 
    }"
    HAS_SIGNALS
)
if (NOT HAS_SIGNALS)
    list(APPEND LRS_ADDITIONAL_DEFINES SIGNALS)
endif()

check_source_compiles(
    C
    "int main() {
        get_time();
        double t = get_time();
        return 0;
    }"
    HAS_TIMES
)
if (NOT HAS_TIMES)
    list(APPEND LRS_ADDITIONAL_DEFINES TIMES)
endif()
check_source_compiles(
    C
    "
    #include <signal.h>
    #include <stdlib.h>		/* llabs */
    #include <unistd.h>
    int main() { 
        long long x = llabs(-5); 
        return 0; 
    }"
    HAS_SIGNALS
)



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
add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
list(APPEND LRS_ADDITIONAL_DEFINES __STDC__ PATH_MAX=1024 WIN)
add_library(mockup src/libgen-mockup.cpp
                   src/tempfile.cpp
                   src/time.cpp
                   src/int128.c)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
endif()

set(MPLRSOBJ2 lrslib2-mplrs.o lrslong2-mplrs.o)
set(LRSOBJ2 lrslib2.o lrslong2.o)
set(LRSOBJ lrs.o lrslong1.o lrslib1.o lrslibgmp.o lrsgmp.o)
set(LRSOBJMP lrs.o lrslong1.o lrslong2.o lrslib1.o lrslib2.o lrslibmp.o lrsmp.o)
set(MPLRSOBJ lrslong1-mplrs.o lrslib1-mplrs.o lrslibgmp-mplrs.o lrsgmp-mplrs.o mplrs.o)
set(MPLRSOBJ64 lrslong1-mplrs.o lrslib1-mplrs.o lrslibgmp-mplrs.o lrsgmp-mplrs.o mplrs64.o)


# lrs.o: lrs.c
#	 $(CC) ${CFLAGS} -DMA ${BITS} -c -o lrs.o lrs.c
add_library(lrs.o OBJECT ${lrslibsrc_SOURCE_DIR}/lrs.c)
target_compile_definitions(lrs.o PRIVATE MA B128 ${LRS_ADDITIONAL_DEFINES})
target_include_directories(lrs.o PRIVATE ${lrslibsrc_SOURCE_DIR}/${ARITH})

# lrs64.o: lrs.c
#     $(CC) ${CFLAGS} -DMA -c -o lrs64.o lrs.c
add_library(lrs64.o OBJECT ${lrslibsrc_SOURCE_DIR}/lrs.c)
target_compile_definitions(lrs64.o PRIVATE MA ${LRS_ADDITIONAL_DEFINES})
target_include_directories(lrs64.o PRIVATE ${lrslibsrc_SOURCE_DIR}/${ARITH})

# lrslong1.o: ${ARITH}lrslong.c ${ARITH}lrslong.h
#     $(CC) ${CFLAGS} -DMA -DSAFE -DLRSLONG -c -o lrslong1.o ${ARITH}lrslong.c
add_library(lrslong1.o OBJECT ${lrslibsrc_SOURCE_DIR}/${ARITH}/lrslong.c)
target_compile_definitions(lrslong1.o PRIVATE MA SAFE LRSLONG ${LRS_ADDITIONAL_DEFINES})
target_include_directories(lrslong1.o PRIVATE ${lrslibsrc_SOURCE_DIR}/${ARITH})
set_target_properties(lrslong1.o PROPERTIES C_STANDARD 23 C_STANDARD_REQUIRED YES)

# lrslib1.o: lrslib.c lrslib.h
#     $(CC) ${CFLAGS} ${PLRSFLAGS} -DMA -DSAFE -DLRSLONG -c -o lrslib1.o lrslib.c
add_library(lrslib1.o OBJECT ${lrslibsrc_SOURCE_DIR}/lrslib.c)
target_compile_definitions(lrslib1.o PRIVATE MA SAFE LRSLONG PLRS ${LRS_ADDITIONAL_DEFINES})
target_include_directories(lrslib1.o PRIVATE ${lrslibsrc_SOURCE_DIR}/${ARITH})
target_link_libraries(lrslib1.o PRIVATE OpenMP::OpenMP_C)

# lrslong2.o: ${ARITH}lrslong.c ${ARITH}lrslong.h
#     $(CC) ${CFLAGS} -DMA -DSAFE ${BITS} -DLRSLONG -c -o lrslong2.o ${ARITH}lrslong.c
add_library(lrslong2.o OBJECT ${lrslibsrc_SOURCE_DIR}/${ARITH}/lrslong.c)
target_compile_definitions(lrslong2.o PRIVATE MA SAFE LRSLONG B128 ${LRS_ADDITIONAL_DEFINES})
target_include_directories(lrslong2.o PRIVATE ${lrslibsrc_SOURCE_DIR}/${ARITH})
set_target_properties(lrslong2.o PROPERTIES C_STANDARD 23 C_STANDARD_REQUIRED YES)

# lrslib2.o: lrslib.c lrslib.h
#     $(CC) ${CFLAGS} ${PLRSFLAGS} -DMA -DSAFE ${BITS} -DLRSLONG -c -o lrslib2.o lrslib.c
add_library(lrslib2.o OBJECT ${lrslibsrc_SOURCE_DIR}/lrslib.c)
target_compile_definitions(lrslib2.o PRIVATE MA SAFE LRSLONG B128 PLRS ${LRS_ADDITIONAL_DEFINES})
target_include_directories(lrslib2.o PRIVATE ${lrslibsrc_SOURCE_DIR}/${ARITH})
target_link_libraries(lrslib2.o PRIVATE OpenMP::OpenMP_C)

# lrslibgmp.o: lrslib.c lrslib.h
#     $(CC) ${CFLAGS} ${PLRSFLAGS}  -DMA -I${INCLUDEDIR} -c -o lrslibgmp.o lrslib.c ${GMP}
add_library(lrslibgmp.o OBJECT ${lrslibsrc_SOURCE_DIR}/lrslib.c)
target_compile_definitions(lrslibgmp.o PRIVATE MA GMP PLRS ${LRS_ADDITIONAL_DEFINES})
target_include_directories(lrslibgmp.o PRIVATE ${lrslibsrc_SOURCE_DIR}/${ARITH})
target_link_libraries(lrslibgmp.o PUBLIC PkgConfig::gmp OpenMP::OpenMP_C)

# lrslibmp.o: lrslib.c lrslib.h
#     $(CC) ${CFLAGS} ${PLRSFLAGS} -DMA -DMP -c -o lrslibmp.o lrslib.c
add_library(lrslibmp.o OBJECT ${lrslibsrc_SOURCE_DIR}/lrslib.c)
target_compile_definitions(lrslibmp.o PRIVATE MA MP PLRS ${LRS_ADDITIONAL_DEFINES})
target_include_directories(lrslibmp.o PRIVATE ${lrslibsrc_SOURCE_DIR}/${ARITH})
target_link_libraries(lrslibmp.o PRIVATE OpenMP::OpenMP_C)

# lrsgmp.o: ${ARITH}lrsgmp.c ${ARITH}lrsgmp.h
#     $(CC) ${CFLAGS} -DMA -I${INCLUDEDIR} -c -o lrsgmp.o ${ARITH}lrsgmp.c ${GMP}
add_library(lrsgmp.o OBJECT ${lrslibsrc_SOURCE_DIR}/${ARITH}/lrsgmp.c)
target_compile_definitions(lrsgmp.o PRIVATE MA GMP ${LRS_ADDITIONAL_DEFINES})
target_include_directories(lrsgmp.o PRIVATE ${lrslibsrc_SOURCE_DIR} ${lrslibsrc_SOURCE_DIR}/${ARITH})
target_link_libraries(lrsgmp.o PUBLIC PkgConfig::gmp)

# lrsmp.o: ${ARITH}lrsmp.c ${ARITH}lrsmp.h
#     $(CC) ${CFLAGS} -DMA -DMP -c -o lrsmp.o ${ARITH}lrsmp.c
add_library(lrsmp.o OBJECT ${lrslibsrc_SOURCE_DIR}/${ARITH}/lrsmp.c)
target_compile_definitions(lrsmp.o PRIVATE MA MP ${LRS_ADDITIONAL_DEFINES})
target_include_directories(lrsmp.o PRIVATE ${lrslibsrc_SOURCE_DIR} ${lrslibsrc_SOURCE_DIR}/${ARITH})



add_executable(lrs ${lrslibsrc_SOURCE_DIR}/lrsdriver.c)
target_include_directories(lrs PRIVATE ${lrslibsrc_SOURCE_DIR} ${lrslibsrc_SOURCE_DIR}/${ARITH})
target_link_libraries(lrs PRIVATE 
    $<TARGET_OBJECTS:lrs.o> 
    $<TARGET_OBJECTS:lrslong1.o> 
    $<TARGET_OBJECTS:lrslong2.o> 
    $<TARGET_OBJECTS:lrslib1.o> 
    $<TARGET_OBJECTS:lrslib2.o> 
    $<TARGET_OBJECTS:lrslibgmp.o> 
    $<TARGET_OBJECTS:lrsgmp.o> 
    OpenMP::OpenMP_C PkgConfig::gmp $<$<PLATFORM_ID:Windows>:mockup>
    )
target_compile_definitions(lrs PRIVATE ${WIN_CFLAGS})


add_executable(lrsMP ${lrslibsrc_SOURCE_DIR}/lrsdriver.c)
target_include_directories(lrsMP PRIVATE ${lrslibsrc_SOURCE_DIR} ${lrslibsrc_SOURCE_DIR}/${ARITH})
target_link_libraries(lrsMP PRIVATE
    $<TARGET_OBJECTS:lrs.o> 
    $<TARGET_OBJECTS:lrslong1.o> 
    $<TARGET_OBJECTS:lrslong2.o> 
    $<TARGET_OBJECTS:lrslib1.o> 
    $<TARGET_OBJECTS:lrslib2.o> 
    $<TARGET_OBJECTS:lrslibmp.o> 
    $<TARGET_OBJECTS:lrsmp.o>
    OpenMP::OpenMP_C PkgConfig::gmp $<$<PLATFORM_ID:Windows>:mockup>
    )
target_compile_definitions(lrsMP PRIVATE ${WIN_CFLAGS} )

add_executable(lrsgmp  
    ${lrslibsrc_SOURCE_DIR}/lrs.c 
    ${lrslibsrc_SOURCE_DIR}/lrslib.c 
    ${lrslibsrc_SOURCE_DIR}/${ARITH}/lrsgmp.c 
    ${lrslibsrc_SOURCE_DIR}/lrsdriver.c )
target_include_directories(lrsgmp PRIVATE ${lrslibsrc_SOURCE_DIR} ${lrslibsrc_SOURCE_DIR}/${ARITH})
target_link_libraries(lrsgmp PRIVATE PkgConfig::gmp OpenMP::OpenMP_C $<$<PLATFORM_ID:Windows>:mockup>)
target_compile_definitions(lrsgmp PRIVATE GMP PLRS ${WIN_CFLAGS} ${LRS_ADDITIONAL_DEFINES})

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
