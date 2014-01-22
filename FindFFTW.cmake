# Find the FFTW includes and library.
# 
# This module defines
# FFTW_INCLUDE_DIR, where to locate fftw3.h file
# FFTW_LIBRARIES, the libraries to link against to use fftw3
# FFTW_FOUND.  If false, you cannot build anything that requires fftw3.
# FFTW_LIBRARY, where to find the libfftw3 library.
# FFTW_THREADS_FOUND, If true, can use multithreaded FFTW routines.
# FFTW_THREADS_LIBRARY, where to find the libfftw3_threads library.

set(FFTW_FOUND FALSE)
set(FFTW_THREADS_FOUND FALSE)
if(FFTW_LIBRARY AND FFTW_INCLUDE_DIR)
  set(FFTW_FIND_QUIETLY TRUE)
endif()

find_path(FFTW_INCLUDE_DIR fftw3.h
  $ENV{FFTW_DIR}/include
  $ENV{FFTW3} $ENV{FFTW3}/include $ENV{FFTW3}/api
  /usr/local/include
  /usr/include
  /opt/fftw3/include
  DOC "Specify the directory containing fftw3.h"
)

find_library(FFTW_LIBRARY NAMES fftw3 fftw3-3 PATHS
  $ENV{FFTW_DIR}/lib
  $ENV{FFTW3} $ENV{FFTW3}/lib $ENV{FFTW3}/.libs
  /usr/local/lib
  /usr/lib 
  /opt/local/lib
  /opt/fftw3/lib
  DOC "Specify the fftw3 library here."
)

find_library(FFTW_THREADS_LIBRARY NAMES fftw3_threads fftw3-3_threads PATHS
  $ENV{FFTW_DIR}/lib
  $ENV{FFTW3} $ENV{FFTW3}/lib $ENV{FFTW3}/.libs
  /usr/local/lib
  /usr/lib 
  /opt/local/lib
  /opt/fftw3/lib
  DOC "Specify the fftw3 threads library here."
)

# make sure pthreads is present, and add to list of threads libraries
if(FFTW_THREADS_LIBRARY)
  find_package(Threads)
  if(CMAKE_USE_PTHREADS_INIT)
    set(FFTW_THREADS_LIBRARY ${FFTW_THREADS_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})
  else()
    message(WARNING "Found FFTW threads but no pthreads! Will not use FFTW threads.")
    unset(FFTW_THREADS_LIBRARY)
  endif()
endif()


if(FFTW_INCLUDE_DIR AND FFTW_LIBRARY)
  set(FFTW_FOUND TRUE)
  if(NOT FFTW_FIND_QUIETLY)
     message(STATUS "Found fftw3 includes at ${FFTW_INCLUDE_DIR}")
     message(STATUS "Found fftw3 library at ${FFTW_LIBRARY}")
  endif()
  
  if(FFTW_THREADS_LIBRARY)
    set(FFTW_THREADS_FOUND TRUE)
    if(NOT FFTW_FIND_QUIETLY)
      message(STATUS "Found fftw3 threads and pthreads libraries at ${FFTW_THREADS_LIBRARY}")
    endif()
  endif()
endif()

mark_as_advanced(FFTW_FOUND FFTW_LIBRARY FFTW_INCLUDE_DIR FFTW_THREADS_FOUND)
if (FFTW_THREADS_LIBRARY)
  mark_as_advanced(FFTW_THREADS_LIBRARY)
  set(FFTW_LIBRARIES ${FFTW_LIBRARY} ${FFTW_THREADS_LIBRARY})
else()
  set(FFTW_LIBRARIES ${FFTW_LIBRARY})
endif()
