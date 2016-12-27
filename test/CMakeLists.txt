include_directories (${TEST_SOURCE_DIR}/src)

add_executable (test_nu test_nu.cpp)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  target_link_libraries (test_nu LINK_PUBLIC liknorm)
else()
  target_link_libraries (test_nu LINK_PUBLIC liknorm m)
endif()
