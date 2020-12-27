file(GLOB OBSERVABLES_CXX CONFIGURE_DEPENDS, observables/*.cc)

add_executable(psidm main.cc)

target_sources(psidm
    PRIVATE state.cc 
            domain.cc 
            fftw.cc 
            cosmology.cc
            ic.cc
            convolution_functions.cc
            interaction/poisson_fc.cc
            interaction/screened_poisson_fc.cc
            interaction/periodic_convolution.cc
            interaction/nonlocal_gp.cc
            evolution/kinetic.cc
            evolution/bfsp.cc
            evolution/bffp.cc
            evolution/interaction_potential.cc
            evolution/interaction_potential_magnus.cc
            evolution/interaction_potential_cap_trapezodial.cc
            evolution/concrete_splittings.cc
            "${OBSERVABLES_CXX}"
)
# don't use -isystem, just -I
set_target_properties(psidm PROPERTIES 
    NO_SYSTEM_FROM_IMPORTED TRUE
)

target_compile_options(psidm
    PRIVATE $<$<CONFIG:Release>:
            -ffast-math
            -march=native>
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
            GSL::gsl
)

if(NOT DEFINED ENV{MKLROOT})
    target_link_libraries(psidm FFTW3::fftw3)
endif()
