cmake_minimum_required(VERSION 3.16)
project(DPSpeedSolver C)

set(CMAKE_C_STANDARD 99)

add_executable(DPSolver2D main.c src/SolverAlgorithm.c inc/SolverAlgorithm.h src/BoundaryLine.c inc/BoundaryLine.h src/MathFunction.c inc/MathFunction.h src/PrintResult.c inc/PrintResult.h inc/SolverStruct.h src/SystemDynamic.c inc/SystemDynamic.h inc/AdaptiveGrid.h src/AdaptiveGrid.c inc/SpeedSolver.h src/SpeedSolver.c inc/ThermalSolver.h src/ThermalSolver.c)