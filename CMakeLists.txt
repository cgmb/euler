cmake_minimum_required(VERSION 2.8.11)
project(euler)
option(OVERDRIVE "Enable platform-specific optimizations" ON)
option(SHERLOCK "Enable expensive debug checks" OFF)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -ffast-math")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-exceptions -fno-rtti")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-result")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=return-type")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=reorder")

if(${OVERDRIVE})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

if(${SHERLOCK})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DAE_FLTDEBUG")
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Developer" OR CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_FORTIFY_SOURCE=2")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_ASSERTIONS")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_DEBUG")
  if(CMAKE_BUILD_TYPE STREQUAL "Developer")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Og")
  endif()
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")
endif()

add_executable(euler
  main.cxx
  world_grid.cxx
  misc/color.cxx
  misc/debug.cxx
  misc/terminal.cxx
  misc/file.cxx
  misc/rng.cxx
)
target_include_directories(euler PRIVATE .)
