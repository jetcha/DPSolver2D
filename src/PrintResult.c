#include "../inc/PrintResult.h"

void printInputInfo(SolverInput *InputPtr, real_T X0, real_T X0_round, uint16_t startIdx, real_T *StateVec,
                    real_T *ControlVec, uint16_t Nx, uint16_t Nu) {

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

    printf("The given initial state: %f\n\n", X0);

    printf("The rounded initial state: %f\n\n", X0_round);

    printf("Starting Index: %d\n\n", startIdx);

    printf("Cost of Infeasibility: %f\n\n\n", InputPtr->SolverLimit.infValue);
}

void printSpeedSolution(SolverInput *InputPtr, real_T X0_round, SolverOutput *OutputPtr) {
    // Make local copies of Grid sizes
    uint16_t Nhrz = HORIZON;
    uint16_t i;

    printf(" -- Output --\n\n");

    printf("Optimal Speed Trajectory: \n");
    printf("Initial Speed: %f \n", X0_round);
    for (i = 0; i < Nhrz; i++) {
        printf("%f ", OutputPtr->Vo[i]);
        if (i % 10 == 9) {
            printf("\n");
        }
    }
    printf("\n");

    printf("Optimal Control Trajectory: \n");
    for (i = 0; i < Nhrz; i++) {
        printf("%f ", OutputPtr->Fo[i]);
        if (i % 10 == 9) {
            printf("\n");
        }
    }
    printf("\n");

    printf("Minimum Total Cost: %f\n", OutputPtr->Cost);
}

void printThermalSolution(SolverInput *InputPtr, real_T X0_round, SolverOutput *OutputPtr) {
    // Make local copies of Grid sizes
    uint16_t Nhrz = RES_THERMAL;
    uint16_t i;

    printf(" -- Output --\n\n");

    printf("Optimal Temperature Trajectory: \n");
    printf("Initial Temperature: %f \n", X0_round);
    for (i = 0; i < Nhrz; i++) {
        printf("%f ", OutputPtr->To[i]);
        if (i % 10 == 9) {
            printf("\n");
        }
    }
    printf("\n");

    printf("Optimal Control Trajectory: \n");
    for (i = 0; i < Nhrz; i++) {
        printf("%f ", OutputPtr->Qo[i]);
        if (i % 10 == 9) {
            printf("\n");
        }
    }
    printf("\n");

    printf("Minimum Total Cost: %f\n", OutputPtr->Cost);
}