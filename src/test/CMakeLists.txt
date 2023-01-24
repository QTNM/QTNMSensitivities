add_executable(TestSensitivity TestSensitivity.cpp)
target_link_libraries(TestSensitivity PRIVATE Experiment Formatting ${ROOT_LIBRARIES})
target_include_directories(TestSensitivity PRIVATE ${PROJECT_SOURCE_DIR}/include)

add_executable(Project8Sensitivities Project8Sensitivities.cpp)
target_link_libraries(Project8Sensitivities PRIVATE Experiment Formatting ${ROOT_LIBRARIES})
target_include_directories(Project8Sensitivities PRIVATE ${PROJECT_SOURCE_DIR}/include)