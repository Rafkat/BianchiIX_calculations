cmake_minimum_required(VERSION 3.28)
project(odeCalculations)

set(CMAKE_CXX_STANDARD 17)

add_executable(odeCalculations main.cpp
        nr/nr3.h
        nr/odeint.h
        nr/stepperdopr5.h
        nr/stepper.h
        nr/stepperdopr853.h
        nr/stepperross.h
        nr/ludcmp.h
        nr/stepperbs.h
        constants.h
        bianchiIX/bianchiIX_isotropy.h
        bianchiIX/bianchiIX_isotropy.cpp
        bianchiIX/bianchiIX_anisotropy.h
        bianchiIX/bianchiIX_anisotropy.cpp
        initial_conditions/BIX_isotropic_initial_conditions.h
        initial_conditions/BIX_anisotropic_initial_conditions.h
        output/output_writer.h
        bianchiV/bianchiV_isotropy.h
        bianchiV/bianchiV_isotropy.cpp
        bianchiV/bianchiV_anisotropy.h
        bianchiV/bianchiV_anisotropy.cpp
        initial_conditions/BV_isotropic_initial_conditions.h
        initial_conditions/BV_anisotropic_initial_conditions.h
)
