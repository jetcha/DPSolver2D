#include "SolverStruct.h"
#include "MathFunction.h"
#include "SystemDynamic.h"
#include "BoundaryLine.h"
#include "AdaptiveGrid.h"
#include "PrintResult.h"

/*--- Public Functions ---*/
// Thermal Solver part
void thermalSolver(SolverInput *InputPtr, DynParameter *ParaPtr, EnvFactor *EnvPtr, SolverOutput *OutputPtr,
                   Bridge *BridgePtr, real_T X0, real_T Xfmin, real_T Xfmax);