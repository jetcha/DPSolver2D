#include "SolverStruct.h"

/*--- Public Functions ---*/

// Create Box Edges Vector
void createBoxEdges(real_T *BoxEdgesVector, real_T const *StateVector, uint32_t N);

// Create Adaptive State Grid (Speed)
void createSpeedGrid(real_T (*StateGrid)[NV], real_T *StateVec, uint16_t Nx, uint16_t Nhrz);

// Create Adaptive State Grid (Thermal)
void createThermalGrid(real_T (*StateGrid)[NT], real_T *StateVec, uint16_t Nx, uint16_t Nhrz);