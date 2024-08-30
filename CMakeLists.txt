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
        bianchiIX/bianchiIX_anisotropy.h
        bianchiIX/bianchiIX_anisotropy.cpp
        output/output_writer.h
        anisotropic_initial_conditions.h
)
