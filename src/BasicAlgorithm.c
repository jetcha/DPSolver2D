#include "../inc/BasicAlgorithm.h"

/*--- Pubic Function Definition ---*/
void
MagicBox(SolverInput *InputPtr, DynParameter *ParaPtr, EnvFactor *EnvPtr, SolverOutput *OutputPtr, real_T V0, real_T T0,
         real_T Vfmin, real_T Vfmax, real_T Tfmin, real_T Tfmax) {

    /*--- Cascading DP solver fashion ---*/
    // Speed Solver part
    speedSolver(InputPtr, ParaPtr, EnvPtr, OutputPtr, V0, Vfmin, Vfmax);

    // Bridge Connection
    Bridge BridgeStruct;
    initBridge(&BridgeStruct);
    bridgeConnection(&BridgeStruct, OutputPtr, V0);

    // Counter Reset
    counterDynamics = 0;
    counterInterpo = 0;
    counterBound = 0;

    // Thermal Solver part
    thermalSolver(InputPtr, ParaPtr, EnvPtr, OutputPtr, &BridgeStruct, T0, Tfmin, Tfmax);

    // Destruct the Bridge
    freeBridge(&BridgeStruct);
}