enable_testing()
find_package(GTest MODULE REQUIRED)
 
link_libraries(GTest::GTest GTest::Main Boost::boost blaze::blaze
               nlohmann_json::nlohmann_json
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

add_executable(poisson_epsilon poisson_fft_epsilon.cc)
target_sources(poisson_epsilon
    PRIVATE ../src/interaction/poisson_fft_epsilon.cc
            ../src/state.cc
            ../src/domain.cc
            ../src/fftw.cc
)

gtest_discover_tests(entropy phasespacem poisson_epsilon)
