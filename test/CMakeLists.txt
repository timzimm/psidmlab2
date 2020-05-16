enable_testing()
find_package(GTest MODULE REQUIRED)
 
link_libraries(GTest::GTest GTest::Main Boost::boost blaze::blaze
               nlohmann_json::nlohmann_json GSL::gsl
)

if(NOT DEFINED ENV{MKLROOT})
    link_libraries(FFTW3::fftw3)
endif()

add_executable(entropy entropy.cc)
target_sources(entropy 
    PRIVATE ../src/cosmology.cc 
            ../src/state.cc
            ../src/domain.cc
            ../src/fftw.cc
            ../src/convolution_functions.cc
            ../src/observables/entropy.cc
            ../src/observables/phasespace_distribution.cc
)

add_executable(phasespace phasespace_distribution.cc)
target_sources(phasespace
    PRIVATE ../src/cosmology.cc 
            ../src/state.cc
            ../src/domain.cc
            ../src/fftw.cc
            ../src/convolution_functions.cc
            ../src/observables/phasespace_distribution.cc
)

add_executable(periodic_convolution periodic_convolution.cc)
target_sources(periodic_convolution
    PRIVATE ../src/interaction/periodic_convolution.cc
            ../src/state.cc
            ../src/domain.cc
            ../src/fftw.cc
)

add_executable(state state.cc)
target_sources(state
    PRIVATE ../src/state.cc
            ../src/domain.cc
            ../src/convolution_functions.cc
            ../src/fftw.cc
)

gtest_discover_tests(entropy phasespace periodic_convolution state)
