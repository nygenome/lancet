# Install script for directory: /Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/build/src/libbamtools.a")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libbamtools.a" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libbamtools.a")
    execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/ranlib" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libbamtools.a")
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/build/src/bamtools")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/bamtools" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/bamtools")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/bamtools")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelopmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES
    "/Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/src/api/BamAlgorithms.h"
    "/Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/src/api/BamAlignment.h"
    "/Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/src/api/BamAux.h"
    "/Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/src/api/BamConstants.h"
    "/Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/src/api/BamIndex.h"
    "/Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/src/api/BamMultiReader.h"
    "/Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/src/api/BamReader.h"
    "/Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/src/api/BamWriter.h"
    "/Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/src/api/IBamIODevice.h"
    "/Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/src/api/SamConstants.h"
    "/Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/src/api/SamHeader.h"
    "/Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/src/api/SamProgram.h"
    "/Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/src/api/SamProgramChain.h"
    "/Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/src/api/SamReadGroup.h"
    "/Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/src/api/SamReadGroupDictionary.h"
    "/Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/src/api/SamSequence.h"
    "/Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/src/api/SamSequenceDictionary.h"
    "/Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/src/api/api_global.h"
    "/Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/build/src/api/bamtools_api_export.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelopmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api/algorithms" TYPE FILE FILES "/Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/src/api/algorithms/Sort.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelopmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/shared" TYPE FILE FILES "/Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/src/shared/bamtools_global.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelopmentx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig" TYPE FILE FILES "/Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/build/src/bamtools-1.pc")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/build/src/third_party/cmake_install.cmake")
  include("/Users/rmusunuri/CLionProjects/lancet/bamtools-2.5.2/build/src/api/cmake_install.cmake")

endif()

