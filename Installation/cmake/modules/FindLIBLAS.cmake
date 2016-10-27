# - Try to find LIBLAS
# Once done this will define
#
#  LIBLAS_FOUND = LIBLAS_FOUND - TRUE
#  LIBLAS_INCLUDE_DIR - include directory for libLAS
#  LIBLAS_LIBRARIES   - the libraries (as targets)

# first look in user defined locations
find_path(LIBLAS_INCLUDE_DIR
          NAMES liblas/liblas.hpp
          PATHS /usr/local/include
         )

find_library(LIBLAS_LIBRARIES
             NAMES las las_c liblas liblas_c
             PATHS ENV LD_LIBRARY_PATH
                   ENV LIBRARY_PATH
                   /usr/local/lib
            )

if(LIBLAS_LIBRARIES)
  set(LIBLAS_FOUND TRUE)
endif()

