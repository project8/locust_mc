# Install the header files including the generated Config.h (which is in
# the build tree).
file (GLOB HEADERS *.hpp)
install (FILES ${HEADERS} DESTINATION ${INCLUDE_INSTALL_DIR}/RandomLib)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/Config.h
  DESTINATION include/RandomLib)
