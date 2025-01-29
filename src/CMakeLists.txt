set(BIN_SOURCE_FILES 
  main.cpp
)

add_executable(linfin ${BIN_SOURCE_FILES})

target_link_libraries(linfin PUBLIC corax logger CLI11 yaml-cpp)

target_compile_options(linfin PUBLIC $<$<CONFIG:Debug>:-fsanitize=address>)
target_link_options(linfin PUBLIC $<$<CONFIG:Debug>:-fsanitize=address>)
