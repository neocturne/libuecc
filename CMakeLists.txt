cmake_minimum_required(VERSION 2.6)
project(LIBUECC C)

set(CMAKE_MODULE_PATH ${LIBUECC_SOURCE_DIR})

add_subdirectory(src)

FILE(GLOB include_files "${CMAKE_CURRENT_SOURCE_DIR}/include/libuecc/*.h")
install(FILES ${include_files} DESTINATION include/libuecc)
