add_executable(TestSensitivity TestSensitivity.cpp)
target_link_libraries(TestSensitivity PRIVATE Experiment Formatting ${ROOT_LIBRARIES})
target_include_directories(TestSensitivity PRIVATE ${PROJECT_SOURCE_DIR}/include)

add_executable(Project8Sensitivities Project8Sensitivities.cpp)
target_link_libraries(Project8Sensitivities PRIVATE Experiment Formatting ${ROOT_LIBRARIES})
target_include_directories(Project8Sensitivities PRIVATE ${PROJECT_SOURCE_DIR}/include)

add_executable(DopplerBroadening DopplerBroadening.cpp)
target_link_libraries(DopplerBroadening PRIVATE Experiment Formatting ${ROOT_LIBRARIES})
target_include_directories(DopplerBroadening PRIVATE ${PROJECT_SOURCE_DIR}/include)

add_executable(StatsOnly StatsOnly.cpp)
target_link_libraries(StatsOnly PRIVATE Experiment Formatting ${ROOT_LIBRARIES})
target_include_directories(StatsOnly PRIVATE ${PROJECT_SOURCE_DIR}/include)

add_executable(WhitePaperPlots WhitePaperPlots.cpp)
target_link_libraries(WhitePaperPlots PRIVATE Experiment Formatting ${ROOT_LIBRARIES})
target_include_directories(WhitePaperPlots PRIVATE ${PROJECT_SOURCE_DIR}/include)