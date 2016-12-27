include_directories (${TEST_SOURCE_DIR}/src)

add_executable (test_posteriors test_posteriors.cpp)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  target_link_libraries (test_posteriors LINK_PUBLIC basics)
else()
  target_link_libraries (test_posteriors LINK_PUBLIC basics m)
endif()
