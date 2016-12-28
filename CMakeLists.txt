cmake_minimum_required (VERSION 3.2.2)

project(Basics CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

file(STRINGS "VERSION" BASICS_VERSION)

include_directories(.)

add_subdirectory(basics)

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
add_library(basics SHARED $<TARGET_OBJECTS:basics_objs>)
target_link_libraries(basics ${ARMADILLO_LIBRARIES})

set_target_properties(basics
  PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    VERSION ${BASICS_VERSION}
    PUBLIC_HEADER "basics/basics.h")

add_subdirectory(test)

enable_testing()
add_test(NAME TestPosteriors COMMAND test_posteriors)
