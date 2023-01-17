add_library(Formatting)
target_sources(Formatting PRIVATE Formatting.cpp)
target_include_directories(Formatting PUBLIC include)
target_link_libraries(Formatting PUBLIC ${ROOT_LIBRARIES})