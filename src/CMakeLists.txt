add_executable(psidm 
    main.cc 
    state.cc 
    cosmology.cc
    ic.cc
    observables.cc
    convolution_functions.cc
    interaction/poisson_fft.cc
    interaction/gp.cc
    evolution/kinetic.cc
    evolution/interaction_potential_trapezodial.cc
    evolution/concrete_splittings.cc
    evolution/uso_kdk_itp.cc
)

#Link time optimization for release build - if supported by the compiler
include(CheckIPOSupported)
check_ipo_supported(RESULT result OUTPUT output)
if(result)
    set_property(TARGET psidm PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
else()
    message(WARNING "IPO is not supported: ${output}")
endif()

target_compile_options(psidm
    PRIVATE $<$<CONFIG:Release>:
            -ffast-math
            -march=native>
    PRIVATE $<$<CONFIG:Debug>:
            -fno-omit-frame-pointer
            -fsanitize=undefined
            -fsanitize=address>
)

target_link_options(psidm
    PRIVATE $<$<CONFIG:Debug>:
            -fsanitize=undefined
            -fsanitize=address>
)

target_include_directories(psidm 
    PRIVATE ${CMAKE_SOURCE_DIR}/inc
            $<$<CONFIG:Release>:
            ${CBLAS_INCLUDE_DIRS}>
)

target_link_libraries(psidm 
    PRIVATE ${BLAS_LIBRARIES}
            Boost::boost
            blaze::blaze
            $<$<BOOL:OpenMP_CXX_FOUND>:OpenMP::OpenMP_CXX>
            FFTW3::fftw3
            hdf5::hdf5-static 
            hdf5::hdf5_hl-static
            nlohmann_json::nlohmann_json
)
