#include "../inc/MathFunction.h"

#ifdef INTERPOCOUNTER
uint32_t counterInterpo = 0;
#endif // INTERPOCOUNTER

/*--- Private (global) variable ---*/
real_T *CmpArray;

/*--- Private Function ---*/
// Compare Function used in sorting
int cmpFunc(const void *a, const void *b);

/*--- Public Function Definition ---*/
void lookuptable_init(LookupTable *TableData, const real_T *XData, const real_T *YData, uint32_t length) {
    // Initialize struct parameters
    TableData->init_success = 0;
    TableData->Npoints = length;
    TableData->lastIdx = length - 1;
    TableData->XData = (real_T *) malloc(length * sizeof(real_T));
    TableData->YData = (real_T *) malloc(length * sizeof(real_T));

    // Check if the given XData[] is a set of monotonically increasing data
    uint32_t i;
    for (i = 0; i < length - 1; i++) {
        if (XData[i + 1] < XData[i]) {
            printf("Given Data is not monotonically increasing.\n");
            return;
        }
    }

    // Copy XData & YData to the Table
    memcpy(TableData->XData, XData, length * sizeof(real_T));
    memcpy(TableData->YData, YData, length * sizeof(real_T));

    // Flag for successful initialization
    TableData->init_success = 1;
}

void lookuptable_free(LookupTable *TableData) {
    free(TableData->XData);
    free(TableData->YData);
    //free(TableData);
}

real_T interpolation(LookupTable *TableData, real_T *Xn, real_T *Yn, uint32_t length) {
    uint32_t i;

    // Default Yn[i] value
    for (i = 0; i < length; i++) {
#ifdef NAN
        Yn[i] = NAN;
#else
        Yn[i] = 1E60;
#endif
    }

    //
    if (TableData->init_success == 1) {
        uint32_t lastIdx = TableData->lastIdx;               // local copy of last index of the table
        uint32_t scanningIdx;
        uint32_t resultIdx;

        real_T Xcurrent;                                    // Current X value

        uint8_t found;                                      // flag for whether a match has been found
        uint32_t iLeft = 0;                                 // Left and Right index
        uint32_t iRight = lastIdx;

        for (i = 0; i < length; i++) {
            Xcurrent = Xn[i];                               // local copy of the current X value (StateVector)
            found = 0U;                                     // Reset the found flag

            // if Xcurrent exceeds the range of XData
            // Yn[i] = YData[max] or YData[min]
            if (Xcurrent >= TableData->XData[lastIdx]) {
                Yn[i] = TableData->YData[lastIdx];
            } else if (Xcurrent <= TableData->XData[0]) {
                Yn[i] = TableData->YData[0];
            }
                // Binary search to find XData[resultIdx] < Xcurrent < XData[resultIdx + 1]
            else {
                scanningIdx = lastIdx;

                while (found == 0U) {
                    if (Xcurrent < TableData->XData[scanningIdx]) {
                        iRight = scanningIdx - 1U;
                        scanningIdx = (iRight + iLeft) >> 1U;
                    } else if (Xcurrent < TableData->XData[scanningIdx + 1]) {
                        found = 1U;
                        resultIdx = scanningIdx;
                    } else {
                        iLeft = scanningIdx + 1U;
                        scanningIdx = (iRight + iLeft) >> 1U;
                    }
                }

                // Interpolate Cost-to-come value
                Yn[i] = TableData->YData[resultIdx] + \
                    (TableData->YData[resultIdx + 1] - TableData->YData[resultIdx]) /
                    (TableData->XData[resultIdx + 1] - TableData->XData[resultIdx]) * \
                    (Xcurrent - TableData->XData[resultIdx]);

#ifdef INTERPOCOUNTER
                counterInterpo++;
#endif // INTERPOCOUNTER
            }
        }
    }

    return Yn[0];
}

real_T extrapolation(real_T *XData, real_T *YData, real_T *Xn, real_T *Yn, uint8_t Flag) {
    if (Flag) {
        // Extrapolation for the upper point
        Yn[0] = YData[1] + ((YData[1] - YData[0]) / (XData[1] - XData[0])) * (Xn[0] - XData[1]);
    } else {
        // Extrapolation for the lower point
        Yn[0] = YData[0] + ((YData[1] - YData[0]) / (XData[1] - XData[0])) * (Xn[0] - XData[0]);
    }
    return Yn[0];
}

real_T multiInterp(StateTuple *States, real_T *Outputs, real_T Vk, real_T Tl) {

    // Linear equations with 3 variables - ax + by + c = P

    // Calculate the coefficients
    real_T m1 = 1 / States[0].V;
    real_T m2 = 1 / States[1].V;
    real_T m3 = 1 / States[2].V;
    real_T n1 = States[0].T / States[0].V;
    real_T n2 = States[1].T / States[1].V;
    real_T n3 = States[2].T / States[2].V;
    real_T o1 = Outputs[0] / States[0].V;
    real_T o2 = Outputs[1] / States[1].V;
    real_T o3 = Outputs[2] / States[2].V;

    // Use the above coefficients to find a, b, c
    real_T c = ((o1 - o2) / (n1 - n2) - (o2 - o3) / (n2 - n3)) / ((m1 - m2) / (n1 - n2) - (m2 - m3) / (n2 - n3));
    real_T b = (o1 - o2) / (n1 - n2) - ((m1 - m2) / (n1 - n2)) * c;
    real_T a = o1 - n1 * b - m1 * c;

#ifdef INTERPOCOUNTER
    counterInterpo++;
#endif // INTERPOCOUNTER

    // Then we can have the output = aVk + bTk + c
    return (a * Vk + b * Tl + c);

}

uint32_t findNearest(real_T *Vector, real_T Value, uint32_t length) {
    real_T minError = FLT_MAX;
    real_T currentError;
    uint32_t nearestIdx = 0;
    uint32_t i;

    for (i = 0; i < length; i++) {
        currentError = fabs(Vector[i] - Value);

        if (currentError < minError) {
            minError = currentError;
            nearestIdx = i;
        }
    }

    return nearestIdx;
}

uint32_t findMaxLEQ(real_T *Vector, real_T Value, uint32_t length) {
    real_T minError = FLT_MAX;
    real_T currentError;
    uint32_t maxIdx = 0;
    uint32_t i;

    for (i = 0; i < length; i++) {
        currentError = fabs(Vector[i] - Value);

        if ((currentError < minError) && ((Value - Vector[i]) >= 0)) {
            minError = currentError;
            maxIdx = i;
        }
    }

    // If there is no answer
    if(minError == FLT_MAX){
        maxIdx = length;
    }

    return maxIdx;
}

uint32_t findMinGEQ(real_T *Vector, real_T Value, uint32_t length) {
    real_T minError = FLT_MAX;
    real_T currentError;
    uint32_t minIdx = 0;
    uint32_t i;

    for (i = 0; i < length; i++) {
        currentError = fabs(Vector[i] - Value);

        if ((currentError < minError) && ((Value - Vector[i]) <= 0)) {
            minError = currentError;
            minIdx = i;
        }
    }

    // If there is no answer
    if(minError == FLT_MAX){
        minIdx = length;
    }

    return minIdx;
}

real_T findMaxLEQ_speed(StateTuple (*Xnext)[NT][NF][NQ], real_T Value, uint16_t Nv, uint16_t Nt, uint16_t Nf, uint16_t Nq) {
    real_T minError = FLT_MAX;
    real_T currentError;

    real_T maxSpeed;

    uint16_t i, j, k, l;

    for (i = 0; i < Nv; i++) {
        for (j = 0; j < Nt; j++) {
            for (k = 0; k < Nf; k++) {
                for (l = 0; l < Nq; l++) {
                    currentError = fabs(Xnext[i][j][k][l].V - Value);

                    if ((currentError < minError) && ((Value - Xnext[i][j][k][l].V) >= 0)) {
                        minError = currentError;
                        maxSpeed = Xnext[i][j][k][l].V;
                    }
                }
            }
        }
    }

    return maxSpeed;
}

real_T findMinGEQ_speed(StateTuple (*Xnext)[NT][NF][NQ], real_T Value, uint16_t Nv, uint16_t Nt, uint16_t Nf, uint16_t Nq) {
    real_T minError = FLT_MAX;
    real_T currentError;

    real_T minSpeed = 0;

    uint16_t i, j, k, l;

    for (i = 0; i < Nv; i++) {
        for (j = 0; j < Nt; j++) {
            for (k = 0; k < Nf; k++) {
                for (l = 0; l < Nq; l++) {
                    currentError = fabs(Xnext[i][j][k][l].V - Value);

                    if ((currentError < minError) && ((Value - Xnext[i][j][k][l].V) <= 0)) {
                        minError = currentError;
                        minSpeed = Xnext[i][j][k][l].V;
                    }
                }
            }
        }
    }

    return minSpeed;
}

uint32_t findUnique(real_T *Vector, real_T *Output, uint32_t length) {
    real_T diff = 1e-6;
    uint32_t uniqueCount = 0;

    memcpy(Output, Vector, length * sizeof(real_T));

    uint32_t i;
    uint32_t j;

    for (i = 0; i < length; i++) {
        if (!(isnan(Output[i]))) {
            for (j = i + 1; j < length; j++) {
                if (fabs(Vector[i] - Vector[j]) < diff) {
                    // Set duplicated element to be NAN
                    Output[j] = NAN;
                }
            }
            // Only store unique elements
            Output[uniqueCount] = Output[i];
            // record the number of unique elements
            uniqueCount++;
        }
    }
    return uniqueCount;
}

void sortIdx(real_T *Vector, uint32_t *idx, uint32_t length) {
    CmpArray = malloc(length * sizeof(real_T));

    // Make a Global Copy of the given vector
    memcpy(CmpArray, Vector, length * sizeof(real_T));

    // Initialize the index array [0, 1, ..., length-1]
    uint32_t i;
    for (i = 0; i < length; i++) {
        idx[i] = i;
    }

    qsort(idx, length, sizeof(uint32_t), cmpFunc);

    free(CmpArray);
}

void reorderVector(real_T *Vector, uint32_t *idx_Reordered, uint32_t length) {
    real_T *tmpPtr = (real_T *) malloc(length * sizeof(real_T));

    // Make a Local Copy of the given vector
    memcpy(tmpPtr, Vector, length * sizeof(real_T));

    // Reorder the vector based on the given idx_Reordered
    uint32_t i;
    for (i = 0; i < length; i++) {
        Vector[i] = tmpPtr[idx_Reordered[i]];
    }

    free(tmpPtr);
}

/*--- Private Function Definition ---*/
int cmpFunc(const void *a, const void *b) {
    uint32_t ia = *(uint32_t *) a;
    uint32_t ib = *(uint32_t *) b;
    return (int) (CmpArray[ia] > CmpArray[ib]);
}