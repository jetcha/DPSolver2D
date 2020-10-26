#include "SolverStruct.h"
#include "MathFunction.h"
#include "SystemDynamic.h"
#include "BoundaryLine.h"
#include "AdaptiveGrid.h"
#include "PrintResult.h"
#include "time.h"

/*--- Public Functions ---*/
// The main algorithm running DP ---
// SolverInput *InputPtr:		Pointer to input structure
// DynParameter *ParaPtr	    Pointer to parameter structure
// SolverOutput *OutputPtr:		Pointer to output structure in which the solution is stored
// real_T V0:					The initial speed
// real_T T0:                   The initial temperature

void MagicBox(SolverInput *InputPtr, DynParameter *ParaPtr, EnvFactor *EnvPtr, SolverOutput *OutputPtr, real_T V0, real_T T0);