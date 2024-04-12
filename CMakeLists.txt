cmake_minimum_required(VERSION 3.16.3)

project(symbolic_hilbert)

#-----------------------------------------------------------------------
# set(MY_LIBRARY_SEARCH_PATHS 
#     /usr/lib /usr/lib64 /usr/local/lib /opt/homebrew/lib 
#     /usr/lib/x86_64-linux-gnu /usr/local/lib/x86_64-linux-gnu )

# set(MY_HEADER_SEARCH_PATHS 
#     /usr/include /usr/local/include /opt/homebrew/include )

find_library(MEDDLY_LIBRARY NAMES meddly )
find_path(MEDDLY_HEADERS_PATH NAMES meddly/meddly.h )
if (MEDDLY_LIBRARY)
  add_compile_definitions(HAS_MEDDLY)
endif()

# message("MEDDLY_LIBRARY = ${MEDDLY_LIBRARY}")
# message("MEDDLY_HEADERS_PATH = ${MEDDLY_HEADERS_PATH}")

find_path(BOOST_CPP_HEADERS_PATH NAMES boost/graph/sloan_ordering.hpp )
if (BOOST_CPP_HEADERS_PATH)
  add_compile_definitions(HAS_BOOST_CPP)
endif()

#-----------------------------------------------------------------------

set(CMAKE_CXX_STANDARD 17)
# set(_GLIBCXX_DEBUG 1)

# Declare our build target and its source files
add_executable(sym_hilbert   
                    src/dd_operations.cpp  
                    src/dd_to_dot.cpp  
                    src/main.cpp  
                    src/matrix.cpp  
                    src/variable_order.cpp 
                    src/sym_pottier.cpp
                    src/unit_tests.cpp
                )
# Compile and link sym_hilbert to meddly
target_include_directories(sym_hilbert  PRIVATE  ${MEDDLY_HEADERS_PATH})
target_include_directories(sym_hilbert  PRIVATE  ${BOOST_CPP_HEADERS_PATH})
target_link_libraries(sym_hilbert  ${MEDDLY_LIBRARY})

# Type of build
# set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

# Build for profiling
# set(CMAKE_CXX_FLAGS "-pg")
# set(CMAKE_CXX_FLAGS_DEBUG "-pg")
# set(CMAKE_CXX_LINK_FLAGS "-pg")
# message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
# message("CMAKE_CXX_FLAGS_DEBUG = ${CMAKE_CXX_FLAGS_DEBUG}")
# message("CMAKE_CXX_LINK_FLAGS = ${CMAKE_CXX_LINK_FLAGS}")
