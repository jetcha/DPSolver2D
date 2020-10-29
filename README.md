# DPSolver2D

## Software Information:

- Matlab version: 2019b
- C codes IDE: JetBrains CLion v2020.3 (https://www.jetbrains.com/clion/)
- Compiler: MinGW-w64 C/C++ v6.3

## Hardware Information (that used to evaluate the solver):

- CPU: i7-8550U
- RAM: 16GB
- System type: 64-bit Winddows 10 operating system

## File Explanations:

- ./MatlabFiles
  - /createSfunc.m: Integrate the C function that contains the DP algorithm (MagicBox) into a Simulink model
  - /inputSolver.m: Input parameters, including model parameters and constraints, and an example scenario.
  - /showResult.m: Figures of the generated solutions.
  - /ModelSolver.slx: Simulink model of the DP solver
- ./inc (header files)
- ./src (source files)
  - /SolverStruct.h: Public structures used in the entire program.
  - /SolverAlgorithm.h/c: The main DP algorithm - function MagicBox connects everything.
  - /MathFunction.h/c: General mathematical functions (interpolation, sorting, etc.).
  - /SystemDynamics.h/c: System dynamics calculation.
  - /BoundaryLine.h/c: Boundary line calculation for the vechile speed.
  - /AdaptiveGrid.h/c: Generate state boxes for the adaptive state grid method.
  - /PrintResult.h/c (Not used): Record the output results into a txt file.
- ./main.c: It is only used when running the codes on the C platform.

## Example of using this solver

1. Run createSfunc.m to compile the C codes into a S-function.
2. Run inputSolver.m to set the input parameters and the input scenario.
3. Run ModelSolver.slx to generate the solutions.
4. Run showResult.m to see the results.

## Important Notes:

- The S-function <b>has to be recompiled</b> when the following parameter(s) is(are) changed:
  - Grid sizes: Nv, Nf, Nt, Nq (in inputSolver.m) / NV, NF, NT, NQ (in SolverStruct.c).
  - Length of horizon: Nhrz (in inputSolver.m)/ HORIZON (in SolverStruct.c)
  - Number of iterations: NUM_IDP (in SolverStruct.c)
  - Shrinking rate: GAMMA (in SolverStruct.c)
- The parameters of the Grid sizes and the Length of horizon need to be kept the <b>same</b> in inputSolver.m and SolverStruct.c.
- The boundary line calculation highly depends on the model equations, so it may not be easy to modify this part when the models are very different. It can be turned off by using '#define NOBOUND' in SolverStruct.c.
- To the best of my knowledge, the C codes <b>does not</b> work with Visual Studio compiler due to the use of Variable Length Array (VLA).
