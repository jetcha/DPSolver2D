#include "../inc/SolverAlgorithm.h"

/*--- Global Variables ---*/
uint16_t Nv;                                        // Problem sizes
uint16_t Nf;
uint16_t Nt;
uint16_t Nq;
uint16_t Nhrz;

Boundary BoundaryStruct;                            // boundary line structure

real_T dV = 0;                                      // grid resolution
real_T dT = 0;
real_T dF = 0;
real_T dQ = 0;

/*--- External Variables ---*/
real_T Vinitial;                                    // The initial state
real_T Tinitial;

/*--- Static Variables ---*/
static real_T *SpeedVec;                            // vector of all the speed state values on the grid
static real_T *ForceVec;                            // vector of all the admissible force control values
static real_T *TempVec;                             // vector of all the cabin temperature state values on the grid
static real_T *InletVec;                            // vector of all the admissible inlet heat control values

static SolverInput *SolverInputPtr;                 // local copy of the input pointer
static DynParameter *ParameterPtr;                  // local copy of the parameter pointer
static EnvFactor *EnvFactorPtr;                     // local copy of the Env factor pointer

static uint16_t V0_index;                           // the index of the rounded initial states
static uint16_t T0_index;

#ifdef ITERATIVEDP
static real_T (*SpeedGrid)[NV];                     // Speed grid for IDP
static real_T (*TempGrid)[NT];                      // Temp grid for IDP
static real_T (*ForceGrid)[NF];                     // Force Grid for IDP
static  real_T (*InletGrid)[NQ];                    // Inlet Grid for IDP
#endif

#ifdef ADAPTIVEGRID
static real_T *BoxEdges;                            // Box Edges for adaptive grid method
static real_T (*StateGrid)[NV];                     // Adaptive State Grid
#endif

/*--- Private Structure ---*/
typedef struct {
    real_T (*CostToCome)[NT];                       // The most recent cost-to-come vector [NV][NT]
    Coordinate *startIdx;                           // All the possible starting indexes [NV][NT]
    uint16_t Nstart;                                // Number of possible starting points
    Coordinate (*Xn)[NV][NT];                         // Optimal State Trajectory [Nhrz][NV][NT]
    ControlTuple (*Un)[NV][NT];                     // Optimal Control Policy [Nhrz][NV][NT]
}
        Solution;

typedef struct {
    real_T (*arcCost)[NT][NV][NT];                  // Vector with arc costs from one to another state [NV][NT][NV][NT]
    ControlTuple (*arcU)[NT][NV][NT];               // Vector with instant controls from one to another state [NV][NT][NV][NT]
#ifdef ADAPTIVEGRID
    real_T (*arcX)[NT][NV][NT];                     // (Used in Adaptive Grid) Vector that notes the possible state to state [NV][NT][NV][NT]
#endif
}
        ArcProcess;


/*--- Private Functions ---*/
/*--- Using static to restrict functions in this file ---*/

// Initialize Solution Pointer
static void solutionStruct_init(Solution *SolutionPtr);

// Free Solution Pointer
static void solutionStruct_free(Solution *SolutionPtr);

// Find the optimal trajectory based on the Cost-to-Come Matrix calculated in calculate_min_costTocome function
static void
findSolution(SolverOutput *OutputPtr, Solution *SolutionPtr, real_T Vfmin, real_T Vfmax, real_T Tfmin, real_T Tfmax);

// Provide the Initial State to the Solution Pointer
static void x0_init(Solution *SolutionPtr, real_T V0, real_T T0);

// Update all the feasible starting states at the current step
static void updateStartX(Solution *SolutionPtr, real_T (*CostToCome)[NT]);

// Initialize Arc Process Pointer
static void arcStruct_init(ArcProcess *ArcPtr);

// Free Arc Process Pointer
static void arcStruct_free(ArcProcess *ArcPtr);

// Calculate the Arc Costs
static void calculate_arc_cost(ArcProcess *ArcPtr, uint16_t N);

// Calculate the Cost-to-Come Values at each state
static void calculate_costTocome(Solution *SolutionPtr, uint16_t N);


/*--- Pubic Function Definition ---*/
void MagicBox(SolverInput *InputPtr, DynParameter *ParaPtr, EnvFactor *EnvPtr, SolverOutput *OutputPtr, real_T V0,
              real_T T0, real_T Vfmin, real_T Vfmax, real_T Tfmin, real_T Tfmax) {

    // Local copy of the input pointer
    SolverInputPtr = InputPtr;
    ParameterPtr = ParaPtr;
    EnvFactorPtr = EnvPtr;

    // Global Copy of X0 (External)
    Vinitial = V0;
    Tinitial = T0;

    // Global copy of the problem size (Speed)
    Nv = SolverInputPtr->GridSize.Nv;
    Nf = SolverInputPtr->GridSize.Nf;
    Nt = SolverInputPtr->GridSize.Nt;
    Nq = SolverInputPtr->GridSize.Nq;
    Nhrz = SolverInputPtr->GridSize.Nhrz;

    // Allocate memory to State Vector and Control Vector
    SpeedVec = malloc(Nv * sizeof(real_T));
    ForceVec = malloc(Nf * sizeof(real_T));
    TempVec = malloc(Nt * sizeof(real_T));
    InletVec = malloc(Nq * sizeof(real_T));

    uint16_t i, num;

    // Pass Parameters to SystemDynamics in order to perform calculation
    PassParameters(SolverInputPtr, ParameterPtr, EnvFactorPtr);

    // Initialize State Vector and Control Vector (based on given [Vmin, Vmax], [Fmin, Fmax], [Tmin, Tmax], [Qmin, Qmax])
    createLineSpace(SpeedVec, SolverInputPtr->Constraint.Vmin, SolverInputPtr->Constraint.Vmax, Nv);
    createLineSpace(ForceVec, SolverInputPtr->Constraint.Fmin, SolverInputPtr->Constraint.Fmax, Nf);
    createLineSpace(TempVec, SolverInputPtr->Constraint.Tmin, SolverInputPtr->Constraint.Tmax, Nt);
    createLineSpace(InletVec, SolverInputPtr->Constraint.Qmin, SolverInputPtr->Constraint.Qmax, Nq);



#ifdef ADAPTIVEGRID
    // Initialize Box Edges
    BoxEdges = malloc((Nx + 1) * sizeof(real_T));
    createBoxEdges(BoxEdges, StateVec, Nx);
    // Initialize the Adaptive State Grid
    StateGrid = malloc(sizeof(real_T[Nhrz + 1][Nx]));
    createSpeedGrid(StateGrid, StateVec, Nx, Nhrz);
#endif

    // Obtain the Boundary Line
#ifdef CUSTOMBOUND
    // Calculate the boundary lines before hand
    initSpeedBoundary(&BoundaryStruct);
    customSpeedBoundary(&BoundaryStruct, SolverInputPtr, ParameterPtr, EnvFactorPtr, V0);
#endif

#ifdef ITERATIVEDP
    // Allocate memory to State Grid
    SpeedGrid = malloc(sizeof(real_T[Nhrz + 1][Nv]));
    TempGrid = malloc(sizeof(real_T[Nhrz + 1][Nt]));
    ForceGrid = malloc(sizeof(real_T[Nhrz + 1][Nf]));
    InletGrid = malloc(sizeof(real_T[Nhrz + 1][Nf]));

    // Initialize the State Grid
    for (i = 0; i <= Nhrz; i++) {
        memcpy(SpeedGrid[i], SpeedVec, Nv * sizeof(real_T));
        memcpy(TempGrid[i], TempVec, Nt * sizeof(real_T));
        memcpy(ForceGrid[i], ForceVec, Nf * sizeof(real_T));
        memcpy(InletGrid[i], ForceVec, Nq * sizeof(real_T));
    }

    // Iterative DP with NUM_IDP loops
    for(num = 0; num < NUM_IDP; num++){

        // Initialize Solution Structure
        Solution SolutionStruct;
        solutionStruct_init(&SolutionStruct);

        // Give the initial state X0 to the solution structure
        x0_init(&SolutionStruct, V0, T0);

        V0_index = SolutionStruct.startIdx[0].X;
        T0_index = SolutionStruct.startIdx[0].Y;

        real_T V0_round = SpeedGrid[0][SolutionStruct.startIdx[0].X];
        real_T T0_round = TempGrid[0][SolutionStruct.startIdx[0].Y];

        // Print Input Info
        printf("The given initial Speed: %f\n", V0);
        printf("The starting Speed (rounded): %f\n", V0_round);
        printf("The starting Speed index: %d\n", V0_index);
        printf("The given initial Temp: %f\n", T0);
        printf("The starting Temp (rounded): %f\n", T0_round);
        printf("The starting Temp index: %d\n", T0_index);
        printf("\n");

        printf(" (Speed Part) \n");
        printInputInfo(SpeedGrid[0], ForceGrid[0], Nv, Nf);
        printf(" (Thermal Part) \n");
        printInputInfo(TempGrid[0], InletGrid[0], Nt, Nq);

        // Find the minimum Cost-to-come value step by step
        for (i = 0; i < Nhrz; i++) {
            calculate_costTocome(&SolutionStruct, i);
        }

        // Retrieve the optimal solution
        findSolution(OutputPtr, &SolutionStruct, Vfmin, Vfmax, Tfmin, Tfmax);

        // Free the solution struct memory
        solutionStruct_free(&SolutionStruct);

        // New step size of the state grid, shrinking by half
        dV = (SpeedGrid[0][1] - SpeedGrid[0][0]) / GAMMA;
        dT = (TempGrid[0][1] - TempGrid[0][0]) / GAMMA;
        dF = (ForceGrid[0][1] - ForceGrid[0][0]) / GAMMA;
        dQ = (InletGrid[0][1] - InletGrid[0][0]) / GAMMA;

        printf("\n");
        printf("dV: %f\n", dV);
        printf("dT: %f\n", dT);
        printf("dF: %f\n", dF);
        printf("dQ: %f\n", dQ);
        printf("\n");

        real_T Vmax, Vmin, Tmax, Tmin, Fmax, Fmin, Qmax, Qmin;

        // state grid at the initial step (not necessary)
        Vmax = V0 + dV * Nv / 2;
        Vmin = V0 - dV * Nv / 2;
        Tmax = T0 + dT * Nt / 2;
        Tmin = T0 - dT * Nt / 2;

        createLineSpace(SpeedGrid[0], Vmin, Vmax, Nv);
        createLineSpace(TempGrid[0], Tmin, Tmax, Nt);

        // Generate the new state grid based on the previous run
        for (i = 0; i < Nhrz; i++){

            // Calculate the ranges of interest based on the midpoints
            Vmax = OutputPtr->Vo[i] + dV * Nv / 2;
            Vmin = OutputPtr->Vo[i] - dV * Nv / 2;
            Tmax = OutputPtr->To[i] + dT * Nt / 2;
            Tmin = OutputPtr->To[i] - dT * Nt / 2;
            Fmax = OutputPtr->Fo[i] + dF * Nf / 2;
            Fmin = OutputPtr->Fo[i] - dF * Nf / 2;
            Qmax = OutputPtr->Qo[i] + dQ * Nq / 2;
            Qmin = OutputPtr->Qo[i] - dQ * Nq / 2;

            // Create the new state grids
            createLineSpace(SpeedGrid[i + 1], Vmin, Vmax, Nv);
            createLineSpace(TempGrid[i + 1], Tmin, Tmax, Nt);
            createLineSpace(ForceGrid[i], Fmin, Fmax, Nf);
            createLineSpace(InletGrid[i], Qmin, Qmax, Nq);
        }
    }

    free(SpeedGrid);
    free(TempGrid);
    free(ForceGrid);
    free(InletGrid);

#else
    // Initialize Solution Structure
    Solution SolutionStruct;
    solutionStruct_init(&SolutionStruct);

    // Give the initial state X0 to the solution structure
    x0_init(&SolutionStruct, V0, T0);

    V0_index = SolutionStruct.startIdx[0].X;
    T0_index = SolutionStruct.startIdx[0].Y;

    real_T V0_round = SpeedVec[SolutionStruct.startIdx[0].X];
    real_T T0_round = TempVec[SolutionStruct.startIdx[0].Y];

    // Print Input Info
    printf("The given initial Speed: %f\n", V0);
    printf("The starting Speed (rounded): %f\n", V0_round);
    printf("The starting Speed index: %d\n", V0_index);
    printf("The given initial Temp: %f\n", T0);
    printf("The starting Temp (rounded): %f\n", T0_round);
    printf("The starting Temp index: %d\n", T0_index);
    printf("\n");

    printf(" (Speed Part) \n");
    printInputInfo(SpeedVec, ForceVec, Nv, Nf);
    printf(" (Thermal Part) \n");
    printInputInfo(TempVec, InletVec, Nt, Nq);

    // Find the minimum Cost-to-come value step by step
    for (i = 0; i < Nhrz; i++) {
        calculate_costTocome(&SolutionStruct, i);
    }

    // Retrieve the optimal solution
    findSolution(OutputPtr, &SolutionStruct, Vfmin, Vfmax, Tfmin, Tfmax);

    // Free the solution struct memory
    solutionStruct_free(&SolutionStruct);
#endif

#if defined(CUSTOMBOUND)
    // Get the boundary line to the output pointer
    copySpeedBoundary(&BoundaryStruct, OutputPtr);
#endif

    // Print Output Solution
    printf("Initial Speed: %f\n", V0);
    printSpeedSolution(SolverInputPtr, OutputPtr);
    printf("Initial Temperature: %f\n", T0);
    printThermalSolution(SolverInputPtr, OutputPtr);
    printf("Minimum Total Cost: %f\n", OutputPtr->Cost);

#ifdef DYNCOUNTER
    printf("The number of dynamics computation: %d\n", counterDynamics);
#endif // DYNCOUNTER

#ifdef INTERPOCOUNTER
    printf("The number of interpolation computation: %d\n", counterInterpo);
#endif // INTERPOCOUNTER

#ifdef BOUNDCOUNTER
    printf("The number of boundary computation: %d\n", counterBound);
#endif // BOUNDCOUNTER

    // Free the memory
    free(SpeedVec);
    free(ForceVec);
    free(TempVec);
    free(InletVec);
#if defined(NORMALBOUND) || defined(CUSTOMBOUND)
    freeBoundary(&BoundaryStruct);
#endif

#ifdef ADAPTIVEGRID
    free(BoxEdges);
    free(StateGrid);
#endif
}

static void solutionStruct_init(Solution *SolutionPtr) {
    SolutionPtr->CostToCome = malloc(sizeof(real_T[Nv][Nt]));
    SolutionPtr->startIdx = (Coordinate *) malloc(Nv * Nt * sizeof(Coordinate));

    SolutionPtr->Xn = malloc(sizeof(Coordinate[Nhrz][Nv][Nt]));
    SolutionPtr->Un = malloc(sizeof(ControlTuple[Nhrz][Nv][Nt]));

    uint16_t i, j, l;

    // Initialize the most recent cost-to-come value and possible starting points
    for (i = 0; i < Nv; i++) {
        for (j = 0; j < Nt; j++) {
            SolutionPtr->CostToCome[i][j] = SolverInputPtr->SolverLimit.infValue;
            SolutionPtr->startIdx[i * Nt + j].X = 0;
            SolutionPtr->startIdx[i * Nt + j].Y = 0;
        }
    }

    // Initialize Optimal state trajectory and control policy
    for (l = 0; l < Nhrz; l++) {
        for (i = 0; i < Nv; i++) {
            for (j = 0; j < Nt; j++) {
                SolutionPtr->Xn[l][i][j].X = 0;
                SolutionPtr->Xn[l][i][j].Y = 0;
                SolutionPtr->Un[l][i][j].F = NAN;
                SolutionPtr->Un[l][i][j].Q = NAN;
            }
        }
    }

    // Initialize number of possible starting points
    SolutionPtr->Nstart = 0;
}

static void solutionStruct_free(Solution *SolutionPtr) {
    free(SolutionPtr->CostToCome);
    free(SolutionPtr->startIdx);
    free(SolutionPtr->Xn);
    free(SolutionPtr->Un);
}

static void x0_init(Solution *SolutionPtr, real_T V0, real_T T0) {
    // Round the initial state to the closet point in the state vector
    uint16_t x, y;

#ifdef ITERATIVEDP
    x = (uint16_t) findNearest(SpeedGrid[0], V0, Nv);
    y = (uint16_t) findNearest(TempGrid[0], T0, Nt);
#else
    x = (uint16_t) findNearest(SpeedVec, V0, Nv);
    y = (uint16_t) findNearest(TempVec, T0, Nt);
#endif

    SolutionPtr->startIdx[0].X = x;
    SolutionPtr->startIdx[0].Y = y;

    SolutionPtr->Nstart = 1;
}

static void updateStartX(Solution *SolutionPtr, real_T (*CostToCome)[NT]) {
    uint16_t counter = 0;
    uint16_t i, j;

    // Update the possible starting points at the next step
    for (i = 0; i < Nv; i++) {
        for (j = 0; j < Nt; j++) {
            if (CostToCome[i][j] < SolverInputPtr->SolverLimit.infValue) {
                SolutionPtr->startIdx[counter].X = i;
                SolutionPtr->startIdx[counter].Y = j;
                counter++;
            }
        }
    }
    SolutionPtr->Nstart = counter;
}

static void arcStruct_init(ArcProcess *ArcPtr) {
    uint16_t i, j, k, l;

    ArcPtr->arcCost = malloc(sizeof(real_T[Nv][Nt][Nv][Nt]));
    ArcPtr->arcU = malloc(sizeof(ControlTuple[Nv][Nt][Nv][Nt]));
#ifdef ADAPTIVEGRID
    ArcPtr->arcX = malloc(sizeof(real_T[Nx][Nx]));
#endif

    // Initialize the arc costs
    for (i = 0; i < Nv; i++) {
        for (j = 0; j < Nt; j++) {
            for (k = 0; k < Nv; k++) {
                for (l = 0; l < Nt; l++) {
                    ArcPtr->arcCost[i][j][k][l] = SolverInputPtr->SolverLimit.infValue;
                    //real_T a = ArcPtr->arcCost[i * Nt + j * Nv + k * Nt + l][0][0][0];
                    ArcPtr->arcU[i][j][k][l].F = 0;
                    ArcPtr->arcU[i][j][k][l].Q = 0;
#ifdef ADAPTIVEGRID
                    ArcPtr->arcX[i][j] = 0.0;
#endif
                }
            }
        }
    }
}

static void arcStruct_free(ArcProcess *ArcPtr) {
    free(ArcPtr->arcCost);
    free(ArcPtr->arcU);
#ifdef ADAPTIVEGRID
    free(ArcPtr->arcX);
#endif

}

static void calculate_costTocome(Solution *SolutionPtr, uint16_t N)        // (N is iteration index)
{
    // Initialize arc process structure
    ArcProcess ArcStruct;
    arcStruct_init(&ArcStruct);

    // Calculate Arc Costs
    calculate_arc_cost(&ArcStruct, N);

    // Initialize cost-to-come [Nv][Nt]
    real_T (*CostToCome)[Nt] = malloc(sizeof(real_T[Nv][Nt]));
    // used to store all the possible cost-to-come values to reach each state
    real_T (*CostToBeComp)[Nt] = malloc(sizeof(real_T[Nv][Nt]));

    uint16_t i, j;

    // Initialize the Cost-to-Come Value
    for (i = 0; i < Nv; i++) {
        for (j = 0; j < Nt; j++) {
            CostToCome[i][j] = SolverInputPtr->SolverLimit.infValue;
        }
    }

    for (int n = 0; n < SolutionPtr->Nstart; n++) {
        uint16_t startIdxV = SolutionPtr->startIdx[n].X;
        uint16_t startIdxT = SolutionPtr->startIdx[n].Y;

        //printf("#%d - Start index: %d %d\n", N, startIdxV, startIdxT);

        // If it is the initial point
        if (N == 0) {
            // Since it is the first step, simply use arc costs as cost-to-come (from the starting index)
            // TODO: Check if this is the correct way of copying
            memcpy(CostToBeComp, ArcStruct.arcCost[startIdxV][startIdxT][0], Nv * Nt * sizeof(real_T));
        } else {
            // We need to add the cost-to-some values to the arc costs

            for (i = 0; i < Nv; i++) {
                for (j = 0; j < Nt; j++) {
                    CostToBeComp[i][j] = SolutionPtr->CostToCome[startIdxV][startIdxT] +
                                         ArcStruct.arcCost[startIdxV][startIdxT][i][j];
                }
            }
        }

        // Pick the minimum cost to be the cost-to-come value
        for (i = 0; i < Nv; i++) {
            for (j = 0; j < Nt; j++) {
                if (CostToBeComp[i][j] < CostToCome[i][j]) {
                    SolutionPtr->Xn[N][i][j].X = startIdxV;
                    SolutionPtr->Xn[N][i][j].Y = startIdxT;
                    SolutionPtr->Un[N][i][j].F = ArcStruct.arcU[startIdxV][startIdxT][i][j].F;
                    SolutionPtr->Un[N][i][j].Q = ArcStruct.arcU[startIdxV][startIdxT][i][j].Q;
                    CostToCome[i][j] = CostToBeComp[i][j];
#ifdef ADAPTIVEGRID
                    StateGrid[N + 1][j] = ArcStruct.arcX[startIdx][j];
#endif
                }
            }
        }
    }


    // Obtain the possible starting points at the next step.
    updateStartX(SolutionPtr, CostToCome);

    // Check the number of possible steps
    if (SolutionPtr->Nstart == 0) {
        printf("Loop #%d has problem\n", N);
    } else {
        printf("Loop #%d is fine\n", N);
    }

    // Copy the output back
    memcpy(SolutionPtr->CostToCome, CostToCome, Nv * Nt * sizeof(real_T));

    // Free the intermediate memory
    free(CostToCome);
    free(CostToBeComp);
    arcStruct_free(&ArcStruct);
}

static void calculate_arc_cost(ArcProcess *ArcPtr, uint16_t N)    // N is iteration index
{
    uint16_t i, j, k, l;

    // Initialize 4-D arrays [Nv][Nt][Nf][Nq]
    StateTuple (*Xnext)[Nt][Nf][Nq] = malloc(sizeof(StateTuple[Nv][Nt][Nf][Nq]));
    ControlTuple (*Control)[Nt][Nf][Nq] = malloc(sizeof(ControlTuple[Nv][Nt][Nf][Nq]));
    real_T (*ArcCost)[Nt][Nf][Nq] = malloc(sizeof(real_T[Nv][Nt][Nf][Nq]));
    uint8_t (*InfFlag)[Nt][Nf][Nq] = malloc(sizeof(uint8_t[Nv][Nt][Nf][Nq]));

    // Initialize 2-D arrays
    uint16_t (*FeasibleCounter)[Nt] = malloc(sizeof(uint16_t[Nv][Nt]));

    // The copy of the speed vector and temp vector
    real_T *SpeedVecCopy = (real_T *) malloc(Nv * sizeof(real_T));
    real_T * TempVecCopy = (real_T *) malloc(Nt * sizeof(real_T));

#ifdef ITERATIVEDP
    memcpy(SpeedVecCopy, SpeedGrid[N + 1], Nv * sizeof(real_T));
    memcpy(TempVecCopy, TempGrid[N + 1], Nt * sizeof(real_T));
#else
    memcpy(SpeedVecCopy, SpeedVec, Nv * sizeof(real_T));
    memcpy(TempVecCopy, TempVec, Nt * sizeof(real_T));
#endif


    for (i = 0; i < Nv; i++)
    {
        for (j = 0; j < Nt; j++) {
            for (k = 0; k < Nf; k++) {
                for (l = 0; l < Nq; l++) {
                    Xnext[i][j][k][l].V = 0;
                    Xnext[i][j][k][l].T = 0;
                    Control[i][j][k][l].F = 0;
                    Control[i][j][k][l].Q = 0;
                    ArcCost[i][j][k][l] = 0;
                    InfFlag[i][j][k][l] = 1;
                }
            }
        }
    }

#ifdef ADAPTIVEGRID
    // Calculate System Dynamics
    speedDynamics(Nx, Nu, Xnext, ArcCost, InfFlag, StateGrid[N], ControlVec, &BoundaryStruct, N, X0_index);
#elif defined(ITERATIVEDP)
    // Calculate System Dynamics
    systemDynamics(Xnext, ArcCost, InfFlag, SpeedGrid[N], ForceGrid[N], TempGrid[N], InletGrid[N], V0_index, T0_index, &BoundaryStruct, N);
#else
    // Calculate System Dynamics
    systemDynamics(Xnext, ArcCost, InfFlag, SpeedVec, ForceVec, TempVec, InletVec, V0_index, T0_index, &BoundaryStruct,
                   N);
#endif

#ifdef BOUNDCALIBRATION
    // Adjust the state vector to be fit on the boundary lines
    // 1, Get the theoretical upper and lower bounds at Xnext
    real_T Vmax_end = BoundaryStruct.upperBound[N + 1];
    real_T Vmin_end = BoundaryStruct.lowerBound[N + 1];

    // 2, Get the index of state vector that just exceeds the bounds
    uint16_t minIdx = (uint16_t) findMaxLEQ(SpeedVecCopy, Vmin_end, Nv);
    uint16_t maxIdx = (uint16_t) findMinGEQ(SpeedVecCopy, Vmax_end, Nv);

    // 3, Get the actual upper and lower bound within Xnext
    // 4, Shift the upper and lower points on the state vector
    SpeedVecCopy[minIdx] = findMinGEQ_speed(Xnext, Vmin_end, Nv, Nt, Nf, Nq);
    SpeedVecCopy[maxIdx] = findMaxLEQ_speed(Xnext, Vmax_end, Nv, Nt, Nf, Nq);

    // 5, Update the boundary memory
    BoundaryStruct.boundMemo[N][0] = SpeedVecCopy[minIdx];                 // Actual lower bound at step N+1
    BoundaryStruct.boundMemo[N][1] = SpeedVecCopy[maxIdx];                 // Actual upper bound at step N+1
    BoundaryStruct.boundMemo[N][2] = minIdx;                               // The index of the lower bound in the vector
    BoundaryStruct.boundMemo[N][3] = maxIdx;                               // The index of the upper bound in the vector
#endif

#ifdef ADAPTIVEGRID
    // Local copy of box edges (for step N)
    real_T *BoxEdgesCopy = (real_T *) malloc((Nx + 1) * sizeof(real_T));
    memcpy(BoxEdgesCopy, BoxEdges, (Nx + 1) * sizeof(real_T));

    // Shift the min and max Edge to be the same as Calibrated StateVecCopy
    //uint16_t minIdxEdge = findMaxLEQ(BoxEdgesCopy, StateVecCopy[minIdx], (Nx + 1));
    //uint16_t maxIdxEdge = findMinGEQ(BoxEdgesCopy, StateVecCopy[maxIdx], (Nx + 1));
    //BoxEdgesCopy[minIdxEdge] = StateVecCopy[minIdx];
    //BoxEdgesCopy[maxIdxEdge] = StateVecCopy[maxIdx];
#endif

    // Count the number of feasible control signals per state
    uint32_t counter;

    for (i = 0; i < Nv; i++) {
        for (j = 0; j < Nt; j++) {
            counter = 0;
            for (k = 0; k < Nf; k++) {
                for (l = 0; l < Nq; l++) {
                    // If the transition is infeasible
                    if (InfFlag[i][j][k][l] == 1) {
                        ArcCost[i][j][k][l] = SolverInputPtr->SolverLimit.infValue;
                        Xnext[i][j][k][l].V = 0.0;
                        Xnext[i][j][k][l].T = 0.0;
                    } else {
                        *(Xnext[i][j][0] + counter) = Xnext[i][j][k][l];
                        *(ArcCost[i][j][0] + counter) = ArcCost[i][j][k][l];
#ifdef ITERATIVEDP
                        (*(Control[i][j][0] + counter)).F = ForceGrid[N][k];
                        (*(Control[i][j][0] + counter)).Q = InletGrid[N][l];
#else
                        (*(Control[i][j][0] + counter)).F = ForceVec[k];
                        (*(Control[i][j][0] + counter)).Q = InletVec[l];
#endif

                        // Counting the number of feasible input combos
                        counter++;
                    }
                }
            }
            FeasibleCounter[i][j] = counter;

#ifdef ADAPTIVEGRID
            if (FeasibleCounter[i] > 0) {
            // Start and End indexes of the Box of interest (use findMaxLEQ to find endIdx in this case)
            uint16_t startBox = findMaxLEQ(BoxEdgesCopy, Xnext[i][0], (Nx + 1));
            uint16_t endBox = findMaxLEQ(BoxEdgesCopy, Xnext[i][(FeasibleCounter[i] - 1)], (Nx + 1));

            uint16_t k = 0;

            for (j = startBox; j <= endBox; j++) {
                while (k < FeasibleCounter[i]) {
                    // Break when Xnext has entered the next Box
                    if (Xnext[i][k] > BoxEdgesCopy[j + 1]) {
                        break;
                    }

                    // Pick the minimum cost among all the possible costs to the box
                    if (ArcCost[i][k] < ArcPtr->arcCost[i][j]) {
                        ArcPtr->arcCost[i][j] = ArcCost[i][k];
                        ArcPtr->arcU[i][j] = Control[i][k];
                        ArcPtr->arcX[i][j] = Xnext[i][k];
                    }
                    k++;
                }
            }
        }
#endif
        }
    }


#ifndef ADAPTIVEGRID

    // Grid gaps (only initialize in the first run)
    if(dV == 0 && dT == 0){
        dV = SpeedVec[1] - SpeedVec[0];
        dT = TempVec[1] - TempVec[0];
    }

    for (i = 0; i < Nv; i++) {
        for (j = 0; j < Nt; j++) {

            // Arc Cost (after picking the nearest-neighbor) from one state to another (point-to-point)
            real_T *p2pCost = ArcPtr->arcCost[i][j][0];

            // Arc Control values (after picking the nearest-neighbor) from one state to another (point-to-point)
            ControlTuple *p2pControl = ArcPtr->arcU[i][j][0];

            // Point to each state (head)
            StateTuple *Xnext_real = Xnext[i][j][0];
            ControlTuple *Control_real = Control[i][j][0];
            real_T *ArcCost_real = ArcCost[i][j][0];

            for (k = 0; k < Nv; k++) {

#ifdef MULTILINEAR
                for (l = 0; l < Nt; l++) {

                    *(p2pCost + k * Nv + l) = SolverInputPtr->SolverLimit.infValue;
                    (*(p2pControl + k * Nv + l)).F = SolverInputPtr->SolverLimit.infValue;
                    (*(p2pControl + k * Nv + l)).Q = SolverInputPtr->SolverLimit.infValue;

                    // Only Find the arc costs if the speed is within the legal range
#ifdef CUSTOMBOUND
                    if (SpeedVecCopy[k] > BoundaryStruct.upperBound[N + 1] ||
                        SpeedVecCopy[k] < BoundaryStruct.lowerBound[N + 1]) {
                        continue;
                    }
#elif defined(NOBOUND)
                    if (SpeedVecCopy[k] > EnvFactorPtr->Vmax_env[N + 1] || SpeedVecCopy[k] < EnvFactorPtr->Vmin_env[N + 1]) {
                        continue;
                    }
#endif

                    // If there is no feasible transition
                    if (FeasibleCounter[i][j] == 0) {
                        continue;
                    }

                    real_T Distance;
                    real_T vDistance, tDistance;

                    real_T minDistance1 = FLT_MAX;
                    real_T minDistance2 = FLT_MAX;
                    real_T minDistance3 = FLT_MAX;

                    StateTuple nearestState1, nearestState2, nearestState3;
                    ControlTuple nearestControl1, nearestControl2, nearestControl3;
                    real_T nearestCost1, nearestCost2, nearestCost3;

                    for (counter = 0; counter < FeasibleCounter[i][j]; counter++) {
                        // The distance from the grid point to the real point
                        vDistance = fabs((*(Xnext_real + counter)).V - SpeedVecCopy[k]);
                        tDistance = fabs((*(Xnext_real + counter)).T - TempVecCopy[l]);

                        if (vDistance < dV && tDistance < dT) {
                            Distance = vDistance * vDistance + tDistance * tDistance;
                            if (Distance < minDistance1) {
                                minDistance1 = Distance;
                                nearestState1 = *(Xnext_real + counter);
                                nearestControl1 = *(Control_real + counter);
                                nearestCost1 = *(ArcCost_real + counter);
                            }
                        }
                    }

                    for (counter = 0; counter < FeasibleCounter[i][j]; counter++) {
                        // The distance from the grid point to the real point
                        vDistance = fabs((*(Xnext_real + counter)).V - SpeedVecCopy[k]);
                        tDistance = fabs((*(Xnext_real + counter)).T - TempVecCopy[l]);

                        if (vDistance < dV && tDistance < dT) {
                            Distance = vDistance * vDistance + tDistance * tDistance;
                            if (Distance < minDistance2 && Distance > minDistance1) {
                                minDistance2 = Distance;
                                nearestState2 = *(Xnext_real + counter);
                                nearestControl2 = *(Control_real + counter);
                                nearestCost2 = *(ArcCost_real + counter);
                            }
                        }
                    }


                    for (counter = 0; counter < FeasibleCounter[i][j]; counter++) {

                        // The distance from the grid point to the real point
                        vDistance = fabs((*(Xnext_real + counter)).V - SpeedVecCopy[k]);
                        tDistance = fabs((*(Xnext_real + counter)).T - TempVecCopy[l]);

                        if (vDistance < dV && tDistance < dT) {
                            Distance = vDistance * vDistance + tDistance * tDistance;
                            if (Distance < minDistance3 && Distance > minDistance2 &&
                                ((*(Xnext_real + counter)).V != nearestState1.V &&
                                 (*(Xnext_real + counter)).T != nearestState1.T)) {
                                minDistance3 = Distance;
                                nearestState3 = *(Xnext_real + counter);
                                nearestControl3 = *(Control_real + counter);
                                nearestCost3 = *(ArcCost_real + counter);
                            }
                        }
                    }


                    if (minDistance3 >= (dV * dT)) {

                        // If there are less than 3 usable points, pick the nearest one
                        if (minDistance1 < (dV * dT)) {
                            *(p2pCost + k * Nv + l) = nearestCost1;
                            (*(p2pControl + k * Nv + l)).F = nearestControl1.F;
                            (*(p2pControl + k * Nv + l)).Q = nearestControl1.Q;
                        }
                    } else {

                        // If there are more or equal to 3 usable points, do the multi linear interpolation
                        // Pick the nearest 3 points
                        StateTuple nearest3states[3] = {nearestState1, nearestState2, nearestState3};
                        real_T nearest3forces[3] = {nearestControl1.F, nearestControl2.F, nearestControl3.F};
                        real_T nearest3inlets[3] = {nearestControl1.Q, nearestControl2.Q, nearestControl3.Q};
                        real_T nearest3Costs[3] = {nearestCost1, nearestCost2, nearestCost3};

//                        real_T costInspect = multiInterp(nearest3states, nearest3Costs, SpeedVecCopy[k], TempVec[l]);
//                        real_T forceInspect = multiInterp(nearest3states, nearest3forces, SpeedVecCopy[k],
//                                                          TempVec[l]);
//                        real_T inletInspect = multiInterp(nearest3states, nearest3inlets, SpeedVecCopy[k],
//                                                          TempVec[l]);

                        *(p2pCost + k * Nv + l) = multiInterp(nearest3states, nearest3Costs, SpeedVecCopy[k],
                                                              TempVecCopy[l]);
                        (*(p2pControl + k * Nv + l)).F = multiInterp(nearest3states, nearest3forces, SpeedVecCopy[k],
                                                                     TempVecCopy[l]);
                        (*(p2pControl + k * Nv + l)).Q = multiInterp(nearest3states, nearest3inlets, SpeedVecCopy[k],
                                                                     TempVecCopy[l]);
                    }
                }
#endif

#ifdef NEARESTNEIGHBOR
                for (l = 0; l < Nt; l++) {
                    real_T minDistance = FLT_MAX;
                    real_T vDistance, tDistance, Distance;

                    *(p2pCost + k * Nv + l) = SolverInputPtr->SolverLimit.infValue;
                    (*(p2pControl + k * Nv + l)).F = SolverInputPtr->SolverLimit.infValue;
                    (*(p2pControl + k * Nv + l)).Q = SolverInputPtr->SolverLimit.infValue;

                    // Only Find the arc costs if the speed is within the legal range
#ifdef CUSTOMBOUND
                    if (SpeedVecCopy[k] > BoundaryStruct.upperBound[N + 1] ||
                        SpeedVecCopy[k] < BoundaryStruct.lowerBound[N + 1]) {
                        continue;
                    }
#elif defined(NOBOUND)
                    if (SpeedVecCopy[k] > EnvFactorPtr->Vmax_env[N + 1] || SpeedVecCopy[k] < EnvFactorPtr->Vmin_env[N + 1]) {
                        continue;
                    }
#endif
                    for (counter = 0; counter < FeasibleCounter[i][j]; counter++) {

                        // The distance from the grid point to the real point
                        vDistance = fabs((*(Xnext_real + counter)).V - SpeedVecCopy[k]);
                        tDistance = fabs((*(Xnext_real + counter)).T - TempVecCopy[l]);

                        if (vDistance < dV && tDistance < dT) {
                            Distance = vDistance * vDistance + tDistance * tDistance;
                            if (Distance < minDistance) {
                                minDistance = Distance;
                                *(p2pCost + k * Nv + l) = *(ArcCost_real + counter);
                                (*(p2pControl + k * Nv + l)).F = (*(Control_real + counter)).F;
                                (*(p2pControl + k * Nv + l)).Q = (*(Control_real + counter)).Q;
                            }
                        }
                    }
                }
#endif
            }
        }
    }


#endif

    free(Xnext);
    free(ArcCost);
    free(InfFlag);
    free(Control);
    free(FeasibleCounter);
    free(SpeedVecCopy);
    free(TempVecCopy);
#ifdef ADAPTIVEGRID
    free(BoxEdgesCopy);
#endif
}

static void
findSolution(SolverOutput *OutputPtr, Solution *SolutionPtr, real_T Vfmin, real_T Vfmax, real_T Tfmin, real_T Tfmax) {

    real_T minCost = SolverInputPtr->SolverLimit.infValue;

    uint16_t i, j;
    int16_t k;
    Coordinate finalIdx;

    // Compare all the possible final cost-to-come values, take the minimum one as the minimum total cost
    for (i = 0; i < Nv; i++) {
        for (j = 0; j < Nt; j++) {
            if (SolutionPtr->CostToCome[i][j] < minCost) {
                minCost = SolutionPtr->CostToCome[i][j];
                finalIdx.X = i;
                finalIdx.Y = j;
            }
        }
    }

    //Minimum total cost
    OutputPtr->Cost = minCost;

    // Memory for the optimal state trajectory
    Coordinate *optimalstateIdx = (Coordinate *) calloc(Nhrz + 1, sizeof(Coordinate));

    // Load the final state index
    optimalstateIdx[Nhrz] = finalIdx;

    // Make sure there is at least one solution
    if (minCost < SolverInputPtr->SolverLimit.infValue) {
        // Retrieve the optimal index trajectory
        k = Nhrz - 1;

        while (k >= 0) {
            optimalstateIdx[k] = SolutionPtr->Xn[k][optimalstateIdx[k + 1].X][optimalstateIdx[k + 1].Y];
            k--;
        }

        // Find the optimal control policy and state trajectory
        for (i = 0; i < Nhrz; i++) {
#ifdef ADAPTIVEGRID
            OutputPtr->Vo[i] = StateGrid[i + 1][optimalstateIdx[i + 1]];
#elif defined(ITERATIVEDP)
            OutputPtr->Vo[i] = SpeedGrid[i + 1][optimalstateIdx[i + 1].X];
            OutputPtr->To[i] = TempGrid[i + 1][optimalstateIdx[i + 1].Y];

//            if (optimalstateIdx[i + 1].X == BoundaryStruct.boundMemo[i][2]) {
//                OutputPtr->Vo[i] = BoundaryStruct.boundMemo[i][0];
//                printf("Hit!!!!\n");
//                printf("%d\n", i);
//            } else if (optimalstateIdx[i + 1].X == BoundaryStruct.boundMemo[i][3]) {
//                OutputPtr->Vo[i] = BoundaryStruct.boundMemo[i][1];
//                printf("Hit!!!!\n");
//                printf("%d\n", i);
//            } else {
//                OutputPtr->Vo[i] = SpeedVec[optimalstateIdx[i + 1].X];
//            }
//
//            OutputPtr->To[i] = TempVec[optimalstateIdx[i + 1].Y];
#else
            OutputPtr->Vo[i] = SpeedVec[optimalstateIdx[i + 1].X];
            OutputPtr->To[i] = TempVec[optimalstateIdx[i + 1].Y];
#endif

#ifdef BOUNDCALIBRATION
            if (optimalstateIdx[i + 1].X == BoundaryStruct.boundMemo[i][2]) {
                OutputPtr->Vo[i] = BoundaryStruct.boundMemo[i][0];
                printf("Hit!!!!\n");
                printf("%d\n", i);
            } else if (optimalstateIdx[i + 1].X == BoundaryStruct.boundMemo[i][3]) {
                OutputPtr->Vo[i] = BoundaryStruct.boundMemo[i][1];
                printf("Hit!!!!\n");
                printf("%d\n", i);
            }
#endif
            OutputPtr->Fo[i] = SolutionPtr->Un[i][optimalstateIdx[i + 1].X][optimalstateIdx[i + 1].Y].F;
            OutputPtr->Qo[i] = SolutionPtr->Un[i][optimalstateIdx[i + 1].X][optimalstateIdx[i + 1].Y].Q;
        }
    }

    printf("%f\n", OutputPtr->Cost);
    // Free the memory
    free(optimalstateIdx);
}