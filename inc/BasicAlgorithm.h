#include "SolverStruct.h"
#include "MathFunction.h"
#include "SpeedSolver.h"
#include "ThermalSolver.h"
#include "BoundaryLine.h"
#include "AdaptiveGrid.h"
#include "PrintResult.h"

/*--- Public Functions ---*/
// The main algorithm running DP ---
// SolverInput *InputPtr:		Pointer to input structure
// DynParameter *ParaPtr	    Pointer to parameter structure
// SolverOutput *OutputPtr:		Pointer to output structure in which the solution is stored
// real_T V0:					The initial speed
// real_T T0:                   The initial temperature
// real_T Vfmax:				Upper bound for the final speed
// real_T Vfmin:				Lower bound for the final speed
// real_T Tfmax:				Upper bound for the final temperature
// real_T Tfmin:				Lower bound for the final temperature

void
MagicBox(SolverInput *InputPtr, DynParameter *ParaPtr, EnvFactor *EnvPtr, SolverOutput *OutputPtr, real_T V0, real_T T0,
         real_T Vfmin, real_T Vfmax, real_T Tfmin, real_T Tfmax);