cmake_minimum_required(VERSION 2.8)
project(euler)
option(OPTIMIZE "Enable platform-specific optimizations")
set(SANITIZER "" CACHE STRING "Enable specified sanitizer")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
include_directories(.)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-result")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math")

if(${OPTIMIZE})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

if(NOT ${SANITIZER} STREQUAL "")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=${SANITIZER} -fno-omit-frame-pointer")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "" OR CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Og")
endif()

set(SOURCES
 main.cxx
 misc/terminal.cxx
 misc/file.cxx
 misc/rng.cxx
)

add_executable(euler ${SOURCES})
