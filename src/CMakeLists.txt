include_directories(${LIBUECC_SOURCE_DIR}/include)

add_library(uecc STATIC ec25519.c ec25519_secret.c)

install(TARGETS uecc
  ARCHIVE DESTINATION lib
)
