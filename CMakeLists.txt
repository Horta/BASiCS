cmake_minimum_required (VERSION 3.2.2)

project(Basics CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

file(STRINGS "VERSION" BASICS_VERSION)

add_subdirectory(basics)

add_library(basics SHARED $<TARGET_OBJECTS:basics_objs>)
