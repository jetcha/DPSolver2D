#include "SolverStruct.h"

#ifndef MATHFUNC
#define MATHFUNC

typedef struct {
    uint32_t Npoints;                // Number of the points in the table
    real_T *XData;                    //
    real_T *YData;                    //
    uint8_t init_success;            // Flag for successful initialization
    uint32_t lastIdx;
}
        LookupTable;

#endif //!MATHFUNC

/*--- Public Functions ---*/
// Initialize Lookup Table
void lookuptable_init(LookupTable *TableData, const real_T *XData, const real_T *YData, uint32_t length);

// Free the Lookup Table memory
void lookuptable_free(LookupTable *TableData);

// Perform interpolation
real_T interpolation(LookupTable *TableData, real_T *Xn, real_T *Yn, uint32_t length);

// Perform extrapolation
real_T extrapolation(real_T *XData, real_T *YData, real_T *Xn, real_T *Yn, uint8_t Flag);

// Find the index with the nearest value in the Vector (based on the given Value)
uint32_t findNearest(real_T *Vector, real_T Value, uint32_t length);

// Find the Maximum value in the Vector that is Less or Equal to the given Value
uint32_t findMaxLEQ(real_T *Vector, real_T Value, uint32_t length);

// Find the Minimum value in the Vector that is Greater or Equal to the given Value
uint32_t findMinGEQ(real_T *Vector, real_T Value, uint32_t length);

// Find the number of unique elements (remove those duplicated elements)
uint32_t findUnique(real_T *Vector, real_T *Output, uint32_t length);

// Sort the index in the way that X[0] to X[length] is monotonically increasing
void sortIdx(real_T *Vector, uint32_t *idx, uint32_t length);

// Reorder the Vector based on the given index - idx_Reordered
void reorderVector(real_T *Vector, uint32_t *idx_Reordered, uint32_t length);

