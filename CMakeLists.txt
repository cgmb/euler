cmake_minimum_required(VERSION 2.8.11)
project(euler)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99 -ffast-math")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror=return-type")

if(CMAKE_BUILD_TYPE STREQUAL "" OR CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Og")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_FORTIFY_SOURCE=2")
endif()

add_executable(euler
  main.c
  misc/terminal.c
  misc/file.c
  misc/rng.c
)
target_include_directories(euler PRIVATE .)
target_link_libraries(euler m)
