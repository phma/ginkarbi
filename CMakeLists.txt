project(ginkarbi)
cmake_minimum_required(VERSION 2.8)

add_executable(ancescomp ancescomp.cpp snpfile.cpp ethnicity.cpp)

install(TARGETS ancescomp DESTINATION bin)
include(TestBigEndian)
include(CheckIncludeFiles)

check_include_files(time.h HAVE_TIME_H)
check_include_files(sys/time.h HAVE_SYS_TIME_H)
check_include_files(sys/resource.h HAVE_SYS_RESOURCE_H)
check_include_files(windows.h HAVE_WINDOWS_H)
test_big_endian(BIGENDIAN)

set(GINKARBI_MAJOR_VERSION 0)
set(GINKARBI_MINOR_VERSION 1)
set(GINKARBI_PATCH_VERSION 1)
set(GINKARBI_VERSION ${GINKARBI_MAJOR_VERSION}.${GINKARBI_MINOR_VERSION}.${GINKARBI_PATCH_VERSION})

include_directories(${PROJECT_BINARY_DIR})
configure_file (config.h.in config.h)

set(CPACK_PACKAGE_VERSION_MAJOR ${GINKARBI_MAJOR_VERSION})
set(CPACK_PACKAGE_VERSION_MINOR ${GINKARBI_MINOR_VERSION})
set(CPACK_PACKAGE_VERSION_PATCH ${GINKARBI_PATCH_VERSION})
set(CPACK_SOURCE_IGNORE_FILES /\\\\.git;.*~)
include(CPack)
