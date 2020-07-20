#include "stdlib.h"
#include "stdio.h"
#include "SolverStruct.h"

/*--- Public Functions ---*/
// Print out the input settings to the screen: Inital state, grids, etc.
void printInputInfo(SolverInput *InputPtr, real_T X0, real_T X0_round, uint16_t startIdx, real_T *StateVec, real_T *ControlVec, uint16_t Nx, uint16_t Nu);

// Print out the solution output to the screen: Optimal speed trajectory, control policy, cost.
void printSpeedSolution(SolverInput *InputPtr, real_T X0_round, SolverOutput *OutputPtr);

// Print out the solution output to the screen: Optimal temperature trajectory, control policy, cost.
void printThermalSolution(SolverInput *InputPtr, real_T X0_round, SolverOutput *OutputPtr);