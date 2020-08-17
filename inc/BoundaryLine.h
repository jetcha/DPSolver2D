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


/*--- Public Functions ---*/
#if defined(CUSTOMBOUND)

// Allocate memory to boundary lines
void initSpeedBoundary(Boundary *BoundaryPtr);

// Copy the boundary line results to the output pointer
void copySpeedBoundary(Boundary *BoundaryPtr, SolverOutput *OutputPtr);

// Free the boundary line memory
void freeBoundary(Boundary *BoundaryPtr);

// DP Optimization - Boundary Line (Speed)
void customSpeedBoundary(Boundary *BoundaryPtr, SolverInput *SolverInputPtr, DynParameter *ParaPtr, EnvFactor *EnvPtr,
                         real_T X0);


#endif