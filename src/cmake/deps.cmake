if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
  # Require at least gcc 4.7
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.7)
    message(FATAL_ERROR "SPAdes requires gcc version 4.7 or later")
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.2)
    message(FATAL_ERROR "SPAdes requires clang version 3.2 or later")
  endif()
else()
  message(WARNING "Unsupported compiler is detected. SPAdes compilation was not tested on it and may fail")
endif()

find_package(OpenMP QUIET)
link_libraries(z)
link_libraries(bz2)
link_libraries(pthread)

add_definitions(-DSEQAN_HAS_ZLIB=1)
add_definitions(-DSEQAN_HAS_BZIP2=1)

set(BOOST_ROOT "${EXT_DIR}/include")
set(Boost_USE_MULTITHREADED ON)
find_package(Boost REQUIRED)
