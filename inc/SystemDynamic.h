#include "SolverStruct.h"
#include "BoundaryLine.h"

/*--- Public Functions ---*/
// Get the Parameter Setting from the main.c
void PassParameters(SolverInput *InputPtr, DynParameter *ModelParaPtr, EnvFactor *EnvFactorPtr);

// Create Line Space [min, max]
void createLineSpace(real_T *Vector, real_T min, real_T max, uint32_t N);

// Create State Plane
void createStatePlane(StateTuple (*Plane)[NT], const real_T *VectorV, const real_T *VectorT, uint16_t Nv, uint16_t Nt);

// Pass the staring plane info
void passStatePlane(StateTuple (*Plane)[NT], uint16_t Nv, uint16_t Nt);

// Calculate the system dynamics
void systemDynamics(StateTuple (*Xnext)[NT][NF][NQ], real_T (*ArcCost)[NT][NF][NQ], uint8_t (*InfFlag)[NT][NF][NQ],
                    real_T const *SpeedVec, real_T const *ForceVec, real_T const *TempVec, real_T const *InletVec,
                    uint16_t V0_index, uint16_t T0_index, Boundary *BoundaryPtr, uint16_t N);
