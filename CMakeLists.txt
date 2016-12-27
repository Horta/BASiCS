cmake_minimum_required (VERSION 3.2.2)

project(Basics CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

file(STRINGS "VERSION" BASICS_VERSION)

add_subdirectory(basics)

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
add_library(basics SHARED $<TARGET_OBJECTS:basics_objs>)
target_link_libraries(basics ${ARMADILLO_LIBRARIES})
