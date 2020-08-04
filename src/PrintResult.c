#include "../inc/PrintResult.h"

void printInputInfo(real_T *StateVec, real_T *ControlVec, uint16_t Nx, uint16_t Nu) {

    uint16_t i;

    printf(" -- Successfully Initialized -- \n\n");

    printf("State grid:\n");
    for (i = 0; i < Nx; i++) {
        printf("%f ", StateVec[i]);

        if (i % 10 == 9) {
            printf("\n");
        }
    }
    printf("\n");

    printf("Control grid:\n");
    for (i = 0; i < Nu; i++) {
        printf("%f ", ControlVec[i]);

        if (i % 10 == 9) {
            printf("\n");
        }
    }
    printf("\n\n");
}

void printSpeedSolution(SolverInput *InputPtr, SolverOutput *OutputPtr) {
    // Make local copies of Grid sizes
    uint16_t Nhrz = HORIZON;
    uint16_t i;

    printf(" -- Output --\n\n");

    printf("Optimal Speed Trajectory: \n");
    for (i = 0; i < Nhrz; i++) {
        printf("%f ", OutputPtr->Vo[i]);
        if (i % 10 == 9) {
            printf("\n");
        }
    }
    printf("\n");

    printf("Optimal Force Trajectory: \n");
    for (i = 0; i < Nhrz; i++) {
        printf("%f ", OutputPtr->Fo[i]);
        if (i % 10 == 9) {
            printf("\n");
        }
    }
    printf("\n");
}

void printThermalSolution(SolverInput *InputPtr, SolverOutput *OutputPtr) {
    // Make local copies of Grid sizes
    uint16_t Nhrz = HORIZON;
    uint16_t i;

    printf(" -- Output --\n\n");

    printf("Optimal Temperature Trajectory: \n");
    for (i = 0; i < Nhrz; i++) {
        printf("%f ", OutputPtr->To[i]);
        if (i % 10 == 9) {
            printf("\n");
        }
    }
    printf("\n");

    printf("Optimal Temp Inlet Trajectory: \n");
    for (i = 0; i < Nhrz; i++) {
        printf("%f ", OutputPtr->Qo[i]);
        if (i % 10 == 9) {
            printf("\n");
        }
    }
    printf("\n");
}