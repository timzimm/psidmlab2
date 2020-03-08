file(GLOB OBSERVABLES_CXX CONFIGURE_DEPENDS, observables/*.cc)

add_executable(psidm 
    main.cc 
    state.cc 
    fftw.cc 
    cosmology.cc
    ic.cc
    convolution_functions.cc
    interaction/poisson_fft.cc
    evolution/kinetic.cc
    evolution/interaction_potential_trapezodial.cc
    evolution/concrete_splittings.cc
    "${OBSERVABLES_CXX}"
)

target_compile_features(psidm PRIVATE cxx_std_17)

# don't use -isystem, just -I
set_target_properties(psidm PROPERTIES 
    NO_SYSTEM_FROM_IMPORTED TRUE
)

target_compile_options(psidm
    PRIVATE $<$<CONFIG:Release>:
            -ffast-math
            -march=native>
            $<$<CONFIG:Debug>:
            -fsanitize=undefined
            -fsanitize=address>
)

target_link_options(psidm
    PRIVATE $<$<CONFIG:Debug>:
            -fsanitize=undefined
            -fsanitize=address>
)

if(${PSIDMLAB_SMP})
    target_compile_definitions(psidm PRIVATE PSIDMLAB_SMP) 
endif()

target_include_directories(psidm PRIVATE inc)

target_link_libraries(psidm 
            Boost::boost
            blaze::blaze
            hdf5::hdf5-static 
            hdf5::hdf5_hl-static
            nlohmann_json::nlohmann_json
)
