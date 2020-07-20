#include "SolverStruct.h"

#ifndef BOUNDARY
#define BOUNDARY

typedef struct {
    real_T *upperBound;                         // Upper Boundary Line
    real_T *lowerBound;                         // Lower Boundary Line

#ifdef BOUNDCALIBRATION
    real_T (*boundMemo)[4];                     // Record the actual bounds over time
    // ... [lower bound], [upper bound], [lower index], [upper index]
#endif
}
        Boundary;

#endif // !BOUNDARY

#ifndef BRIDGE
#define BRIDGE

typedef struct {
    real_T *Pdc;                                    // Vector of Pdc based on the optimal Speed trajectory
    real_T *tDelta;                                 // Vector of delta t based on the optimal Speed trajectory
}
        Bridge;
#endif

/*--- Public Functions ---*/
#if defined(NORMALBOUND) || defined(CUSTOMBOUND)

// Allocate memory to boundary lines
void initSpeedBoundary(Boundary *BoundaryPtr);
void initThermalBoundary(Boundary *BoundaryPtr);

// Copy the boundary line results to the output pointer
void copySpeedBoundary(Boundary *BoundaryPtr, SolverOutput *OutputPtr);

// Copy the boundary line results to the output pointer
void copyThermalBoundary(Boundary *BoundaryPtr, SolverOutput *OutputPtr);

// Free the boundary line memory
void freeBoundary(Boundary *BoundaryPtr);

// Normal Boundary - simply copy the legal speed limits as boundaries
void normalSpeedBoundary(Boundary *BoundaryPtr, EnvFactor *EnvPtr);

// DP Optimization - Boundary Line (Speed)
void customSpeedBoundary(Boundary *BoundaryPtr, SolverInput *SolverInputPtr, DynParameter *ParaPtr, EnvFactor *EnvPtr,
                         real_T X0);

// Normal Boundary - simply copy the required Temperature as boundaries
void normalThermalBoundary(Boundary *BoundaryPtr, EnvFactor *EnvPtr);

// DP Optimization - Boundary Line (Thermal)
void customThermalBoundary(Boundary *BoundaryPtr, SolverInput *SolverInputPtr, DynParameter *ParaPtr, EnvFactor *EnvPtr,
                           Bridge *BridgePtr, real_T X0);

#endif