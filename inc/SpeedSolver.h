#include "SolverStruct.h"
#include "MathFunction.h"
#include "SystemDynamic.h"
#include "BoundaryLine.h"
#include "AdaptiveGrid.h"
#include "PrintResult.h"

/*--- Public Functions ---*/
// Speed Solver part
void speedSolver(SolverInput *InputPtr, DynParameter *ParaPtr, EnvFactor *EnvPtr, SolverOutput *OutputPtr, real_T V0, real_T T0,
                 real_T Vfmin, real_T Vfmax, real_T Tfmin, real_T Tfmax);