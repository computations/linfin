
set(BIN_SOURCE_FILES 
  main.cpp
)

add_executable(linass ${BIN_SOURCE_FILES})

target_link_libraries(linass PUBLIC corax logger CLI11)

target_compile_options(linass PUBLIC $<$<CONFIG:Debug>:-fsanitize=address>)
target_link_options(linass PUBLIC $<$<CONFIG:Debug>:-fsanitize=address>)
