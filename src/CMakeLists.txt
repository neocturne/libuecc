include_directories(${LIBUECC_SOURCE_DIR}/include)

add_library(uecc STATIC ec25519.c)

install(TARGETS uecc
  ARCHIVE DESTINATION lib
)
