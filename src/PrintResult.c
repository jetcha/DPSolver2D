#include "../inc/PrintResult.h"

void fileOutput(SolverOutput *OutputPtr, uint16_t N){
    FILE *fp;
    uint16_t i;
    uint16_t Nhrz = HORIZON;

    char *nameMain;
    char *nameExtension = ".txt";

    // Speed Result
    nameMain = "../txtResult/speedResult_";

    char speedResult[strlen(nameMain)+strlen(nameExtension)+4];
    snprintf(speedResult, sizeof(speedResult), "%s%d%s", nameMain, N, nameExtension);
    fp = fopen(speedResult, "w");
    for (i = 0; i < Nhrz; i++) {
        fprintf(fp, "%f ", OutputPtr->Vo[i]);
    }
    fclose(fp);

    // Temperature Result
    nameMain = "../txtResult/tempResult_";
    char tempResult[strlen(nameMain)+strlen(nameExtension)+4];
    snprintf(tempResult, sizeof(tempResult), "%s%d%s", nameMain, N, nameExtension);

    fp = fopen(tempResult, "w");
    for (i = 0; i < Nhrz; i++) {
        fprintf(fp, "%f ", OutputPtr->To[i]);
    }
    fclose(fp);

    // Force Result
    nameMain = "../txtResult/forceResult_";
    char forceResult[strlen(nameMain)+strlen(nameExtension)+4];
    snprintf(forceResult, sizeof(forceResult), "%s%d%s", nameMain, N, nameExtension);

    fp = fopen(forceResult, "w");
    for (i = 0; i < Nhrz; i++) {
        fprintf(fp, "%f ", OutputPtr->Fo[i]);
    }
    fclose(fp);

    // Inlet Result
    nameMain = "../txtResult/inletResult_";
    char inletResult[strlen(nameMain)+strlen(nameExtension)+4];
    snprintf(inletResult, sizeof(inletResult), "%s%d%s", nameMain, N, nameExtension);

    fp = fopen(inletResult, "w");
    for (i = 0; i < Nhrz; i++) {
        fprintf(fp, "%f ", OutputPtr->Qo[i]);
    }
    fclose(fp);
}