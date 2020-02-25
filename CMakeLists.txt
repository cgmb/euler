cmake_minimum_required(VERSION 2.8.11)
project(euler)
option(OVERDRIVE "Enable platform-specific optimizations" ON)
option(SHERLOCK "Enable expensive debug checks" OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

if(${OVERDRIVE})
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
endif()

if(${SHERLOCK})
  add_definitions(-DAE_FLTDEBUG)
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99 -ffast-math")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror=return-type")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror=switch-enum")

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Og")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_FORTIFY_SOURCE=2")
endif()

add_executable(euler
  main.c
  misc/terminal.c
  misc/file.c
  misc/rng.c
  misc/debug.c
)
target_include_directories(euler PRIVATE .)
target_link_libraries(euler m)
