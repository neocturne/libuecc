include_directories(${LIBUECC_SOURCE_DIR}/include)

set(UECC_SRC ec25519.c ec25519_gf.c)
set(UECC_ABI 0)

add_library(uecc_shared SHARED ${UECC_SRC})
set_target_properties(uecc_shared PROPERTIES
  OUTPUT_NAME uecc
  SOVERSION ${UECC_ABI}
  VERSION ${UECC_ABI}.${PROJECT_VERSION}
)

add_library(uecc_static STATIC ${UECC_SRC})
set_target_properties(uecc_static PROPERTIES
  OUTPUT_NAME uecc
)

install(TARGETS uecc_shared uecc_static
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
)
