#include "SolverStruct.h"
#include "BoundaryLine.h"

/*--- Public Functions ---*/
// Get the Parameter Setting from the main.c
void PassParameters(SolverInput *InputPtr, DynParameter *ModelParaPtr, EnvFactor *EnvFactorPtr);

// Create Line Space [min, max]
void createLineSpace(real_T *Vector, real_T min, real_T max, uint32_t N);

// Calculate the Speed dynamics
void speedDynamics(uint16_t Nx, uint16_t Nu, real_T (*Xnext)[Nu], real_T (*ArcCost)[Nu], uint8_t (*InfFlag)[Nu],
                   real_T const *StateVec, real_T const *ControlVec, Boundary *BoundaryPtr, uint16_t N,
                   uint16_t X0_index);

// Calculate the Thermal dynamics
void thermalDynamics(uint16_t Nx, uint16_t Nu, real_T (*Xnext)[Nu], real_T (*ArcCost)[Nu], uint8_t (*InfFlag)[Nu],
                     real_T const *StateVec, real_T const *ControlVec, Boundary *BoundaryPtr, Bridge *BridgePtr,
                     uint16_t N,
                     uint16_t X0_index);

// Calculate the system dynamics
void systemDynamics(StateTuple (*Xnext)[NT][NF][NQ], real_T (*ArcCost)[NT][NF][NQ], uint8_t (*InfFlag)[NT][NF][NQ],
                    real_T const *SpeedVec, real_T const *ForceVec, real_T const *TempVec, real_T const *InletVec,
                    uint16_t V0_index, uint16_t T0_index, uint16_t N);

// Initialize bridge memory
void initBridge(Bridge *BridgePtr);

// Free bridge memory
void freeBridge(Bridge *BridgePtr);

// Calculate the Pdc[k] and tDelta[k] which are needed for the thermal solver
void bridgeConnection(Bridge *BridgePtr, SolverOutput *OutputPtr, real_T V0);