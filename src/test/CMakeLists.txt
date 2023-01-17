add_executable(TestSensitivity TestSensitivity.cpp)
target_link_libraries(TestSensitivity PRIVATE Experiment ${ROOT_LIBRARIES})
target_include_directories(TestSensitivity PRIVATE ${PROJECT_SOURCE_DIR}/include)