add_library(Experiment)

target_sources(Experiment PRIVATE Experiment.cpp SystematicEffect.cpp)

target_include_directories(Experiment PUBLIC include)