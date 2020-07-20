#include "SolverStruct.h"
#include "MathFunction.h"
#include "SystemDynamic.h"
#include "BoundaryLine.h"
#include "AdaptiveGrid.h"
#include "PrintResult.h"

/*--- Public Functions ---*/
// Speed Solver part
void speedSolver(SolverInput *InputPtr, DynParameter *ParaPtr, EnvFactor *EnvPtr, SolverOutput *OutputPtr, real_T X0,
                 real_T Xfmin, real_T Xfmax);