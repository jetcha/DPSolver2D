#include "../inc/SpeedSolver.h"

/*--- Global Variables ---*/
uint16_t Nv;                                        // Problem sizes
uint16_t Nf;
uint16_t Nt;
uint16_t Nq;
uint16_t Nhrz;
Boundary BoundaryStruct;                            // boundary line structure

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

#ifdef ADAPTIVEGRID
static real_T *BoxEdges;                            // Box Edges for adaptive grid method
static real_T (*StateGrid)[NV];                     // Adaptive State Grid
#endif

/*--- Private Structure ---*/
typedef struct {
    real_T (*CostToCome)[NT];                       // The most recent cost-to-come vector [NV][NT]
    Coordinate (*startIdx)[NT];                     // All the possible starting indexes [NV][NT]
    uint16_t Nstart;                                // Number of possible starting points
    uint16_t (*Xn)[NV][NT];                         // Optimal State Trajectory [Nhrz][NV][NT]
    real_T (*Un)[NV][NT];                           // Optimal Control Policy [Nhrz][NV][NT]
}
        Solution;

typedef struct {
    real_T (*arcCost)[NT][NV][NT];                  // Vector with arc costs from one to another state [NV][NT][NV][NT]
    real_T (*arcU)[NT][NV][NT];                     // Vector with instant controls from one to another state [NV][NT][NV][NT]
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
static void findSolution(SolverOutput *OutputPtr, Solution *SolutionPtr, real_T Xfmin, real_T Xfmax);

// Provide the Initial State to the Solution Pointer
static void x0_init(Solution *SolutionPtr, real_T V0, real_T T0);

// Update all the feasible starting states at the current step
static void updateStartX(Solution *SolutionPtr, real_T *CostToCome);

// Initialize Arc Process Pointer
static void arcStruct_init(ArcProcess *ArcPtr);

// Free Arc Process Pointer
static void arcStruct_free(ArcProcess *ArcPtr);

// Calculate the Arc Costs
static void calculate_arc_cost(ArcProcess *ArcPtr, uint16_t N);

// Calculate the Cost-to-Come Values at each state
static void calculate_costTocome(Solution *SolutionPtr, uint16_t N);


/*--- Pubic Function Definition ---*/
void speedSolver(SolverInput *InputPtr, DynParameter *ParaPtr, EnvFactor *EnvPtr, SolverOutput *OutputPtr, real_T V0,
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

    // Initialize Solution Structure
    Solution SolutionStruct;
    solutionStruct_init(&SolutionStruct);

    // Give the initial state X0 to the solution structure
    x0_init(&SolutionStruct, V0, T0);

    V0_index = SolutionStruct.startIdx[0][0].X;
    T0_index = SolutionStruct.startIdx[0][0].Y;

    real_T V0_round = SpeedVec[SolutionStruct.startIdx[0][0].X];
    real_T T0_round = TempVec[SolutionStruct.startIdx[0][0].Y];

    // Print Input Info
    printf("The given initial Speed: %f\n", V0);
    printf("The starting Speed (rounded): %f\n", V0_round);
    printf("The starting Speed index: %d\n", V0_index);
    printf("The given initial Temp: %f\n", T0);
    printf("The starting Temp (rounded): %f\n", T0_round);
    printf("The starting Temp index: %d\n", T0_index);

    printf(" (Speed Part) \n");
    printInputInfo(SpeedVec, ForceVec, Nv, Nf);
    printf(" (Thermal Part) \n");
    printInputInfo(TempVec, InletVec, Nt, Nq);

    // Obtain the Boundary Line
#ifdef CUSTOMBOUND
    initSpeedBoundary(&BoundaryStruct);
    customSpeedBoundary(&BoundaryStruct, SolverInputPtr, ParameterPtr, EnvFactorPtr, X0);
#elif defined NORMALBOUND
    initSpeedBoundary(&BoundaryStruct);
    normalSpeedBoundary(&BoundaryStruct, EnvFactorPtr);
#endif


    // Find the minimum Cost-to-come value step by step
    uint16_t i;
    uint16_t j;
    for (i = 0; i < Nhrz; i++) {
        calculate_costTocome(&SolutionStruct, i);
    }

    // Retrieve the optimal solution
    findSolution(OutputPtr, &SolutionStruct, Xfmin, Xfmax);

#if defined(NORMALBOUND) || defined(CUSTOMBOUND)
    // Get the boundary line to the output pointer
    copySpeedBoundary(&BoundaryStruct, OutputPtr);
#endif

    // Print Output Solution
    printSpeedSolution(SolverInputPtr, X0_round, OutputPtr);

#ifdef DYNCOUNTER
    printf("(Speed) The number of dynamics computation: %d\n", counterDynamics);
#endif // DYNCOUNTER

#ifdef INTERPOCOUNTER
    printf("(Speed) The number of interpolation computation: %d\n", counterInterpo);
#endif // INTERPOCOUNTER

#ifdef BOUNDCOUNTER
    printf("(Speed) The number of boundary computation: %d\n", counterBound);
#endif // BOUNDCOUNTER

    // Free the memory
    solutionStruct_free(&SolutionStruct);
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
    SolutionPtr->startIdx = malloc(sizeof(Coordinate[Nv][Nt]));

    SolutionPtr->Xn = malloc(sizeof(uint16_t[Nhrz][Nv][Nt]));
    SolutionPtr->Un = malloc(sizeof(real_T[Nhrz][Nv][Nt]));

    uint16_t i, j, l;

    // Initialize the most recent cost-to-come value and possible starting points
    for (i = 0; i < Nv; i++) {
        for (j = 0; j < Nt; j++) {
            SolutionPtr->CostToCome[i][j] = SolverInputPtr->SolverLimit.infValue;
            SolutionPtr->startIdx[i][j].X = 0;
            SolutionPtr->startIdx[i][j].Y = 0;
        }
    }

    // Initialize Optimal state trajectory and control policy
    for (l = 0; l < Nhrz; l++) {
        for (i = 0; i < Nv; i++) {
            for (j = 0; j < Nt; j++) {
                SolutionPtr->Xn[l][i][j] = 0;
                SolutionPtr->Un[l][i][j] = NAN;
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

    x = (uint16_t) findNearest(SpeedVec, V0, Nv);
    y = (uint16_t) findNearest(TempVec, T0, Nt);

    SolutionPtr->startIdx[0][0].X = x;
    SolutionPtr->startIdx[0][0].Y = y;

    SolutionPtr->Nstart = 1;
}

static void updateStartX(Solution *SolutionPtr, real_T *CostToCome) {
    uint16_t counter = 0;
    uint16_t i;

    for (i = 0; i < Nx; i++) {
        if (CostToCome[i] < SolverInputPtr->SolverLimit.infValue) {
            SolutionPtr->startIdx[counter] = i;
            counter++;
        }
    }
    // Update the number of possible starting points at the next step
    SolutionPtr->Nstart = counter;
}

static void arcStruct_init(ArcProcess *ArcPtr) {
    uint16_t i, j, k, l;

    ArcPtr->arcCost = malloc(sizeof(real_T[Nv][Nt][Nv][Nt]));
    ArcPtr->arcU = malloc(sizeof(real_T[Nv][Nt][Nv][Nt]));
#ifdef ADAPTIVEGRID
    ArcPtr->arcX = malloc(sizeof(real_T[Nx][Nx]));
#endif

    // Initialize the arc costs
    for (i = 0; i < Nv; i++) {
        for (j = 0; j < Nt; j++) {
            for (k = 0; k < Nv; k++) {
                for (l = 0; l < Nt; l++) {
                    ArcPtr->arcCost[i][j][k][l] = SolverInputPtr->SolverLimit.infValue;
                    ArcPtr->arcU[i][j][k][l] = 0;
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

    // Initialize cost-to-come vector
    real_T *CostToCome = malloc(Nx * sizeof(real_T));
    // used to store all the possible cost-to-come values to reach each state
    real_T *CostToBeComp = malloc(Nx * sizeof(real_T));

    uint16_t i;
    uint16_t j;

    // Initialize the Cost-to-Come Value at the current State Vector
    for (i = 0; i < Nx; i++) {
        CostToCome[i] = SolverInputPtr->SolverLimit.infValue;
    }

    for (i = 0; i < SolutionPtr->Nstart; i++) {
        uint16_t startIdx = SolutionPtr->startIdx[i];

        // If it is at the initial point
        if (N == 0) {
            // Since it is the first step, simply use arc costs as cost-to-come (from the starting index)
            memcpy(CostToBeComp, ArcStruct.arcCost[startIdx], Nx * sizeof(real_T));
        }
            // If it is not the initial state
        else {
            // we need to add the cost-to-come values to the arc costs
            for (j = 0; j < Nx; j++) {
                CostToBeComp[j] = SolutionPtr->CostToCome[startIdx] + ArcStruct.arcCost[startIdx][j];
            }
        }

        // Pick the minimum cost to be the cost-to-come value
        for (j = 0; j < Nx; j++) {
            if (CostToBeComp[j] < CostToCome[j]) {
                SolutionPtr->Xn[N][j] = startIdx;
                SolutionPtr->Un[N][j] = ArcStruct.arcU[startIdx][j];
                CostToCome[j] = CostToBeComp[j];
#ifdef ADAPTIVEGRID
                StateGrid[N + 1][j] = ArcStruct.arcX[startIdx][j];
#endif
            }
        }
    }

    // Obtain the number of possible starting points at the next step.
    updateStartX(SolutionPtr, CostToCome);

    // Copy the output back
    memcpy(SolutionPtr->CostToCome, CostToCome, Nx * sizeof(real_T));

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
    real_T (*ArcCost)[Nt][Nf][Nq]= malloc(sizeof(real_T[Nv][Nt][Nf][Nq]));
    uint8_t (*InfFlag)[Nt][Nf][Nq] = malloc(sizeof(uint8_t[Nv][Nt][Nf][Nq]));

    // Initialize 2-D arrays
    uint16_t (*FeasibleCounter)[Nt] = malloc(sizeof(uint16_t[Nv][Nt]));

    for (i = 0; i < Nv; i++) {
        for (j = 0; j < Nt; j++) {
            for (i = 0; i < Nf; i++) {
                for (j = 0; j < Nq; j++) {
                    Xnext[i][j][k][l].V = 0;
                    Xnext[i][j][k][l].T = 0;
                    Control[i][j][k][l].F = 0;
                    Control[i][j][k][l].Q = 0;
                    ArcCost[i][j][k][l] = 0;
                    InfFlag[i][j][k][l] = 0;
                }
            }
        }
    }

#ifdef ADAPTIVEGRID
    // Calculate System Dynamics
    speedDynamics(Nx, Nu, Xnext, ArcCost, InfFlag, StateGrid[N], ControlVec, &BoundaryStruct, N, X0_index);
#else
    // Calculate System Dynamics
    systemDynamics(Xnext, ArcCost, InfFlag, SpeedVec, ForceVec, TempVec, InletVec, V0_index, T0_index, N);
#endif

    // Local Copy the State and Control grids
    real_T *StateVecCopy = (real_T *) malloc(Nx * sizeof(real_T));
    memcpy(StateVecCopy, StateVec, Nx * sizeof(real_T));

#ifdef BOUNDCALIBRATION
    // Adjust the state vector to be fit on the boundary lines
    // 1, Get the theoretical upper and lower bounds at Xnext
    real_T Vmax_end = BoundaryStruct.upperBound[N + 1];
    real_T Vmin_end = BoundaryStruct.lowerBound[N + 1];

    // 2, Get the index of state vector that just exceeds the bounds
    uint16_t minIdx = (uint16_t) findMaxLEQ(StateVecCopy, Vmin_end, Nx);
    uint16_t maxIdx = (uint16_t) findMinGEQ(StateVecCopy, Vmax_end, Nx);

    // 3, Get the actual upper and lower bound within Xnext
    uint16_t minIdxReal = (uint16_t) findMinGEQ(Xnext[0], Vmin_end, Nx * Nu);
    uint16_t maxIdxReal = (uint16_t) findMaxLEQ(Xnext[0], Vmax_end, Nx * Nu);

    // 4, Shift the upper and lower points on the state vector
    StateVecCopy[minIdx] = *(Xnext[0] + minIdxReal);
    StateVecCopy[maxIdx] = *(Xnext[0] + maxIdxReal);

    // 5, Update the boundary memory
    BoundaryStruct.boundMemo[N][0] = StateVecCopy[minIdx];                 // Actual lower bound at step N+1
    BoundaryStruct.boundMemo[N][1] = StateVecCopy[maxIdx];                 // Actual upper bound at step N+1
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

    for (i = 0; i < Nx; i++) {
        counter = 0;
        for (j = 0; j < Nu; j++) {
            if (InfFlag[i][j] == 1) {
                ArcCost[i][j] = SolverInputPtr->SolverLimit.infValue;
                Xnext[i][j] = 0.0;
            } else {
                Xnext[i][counter] = Xnext[i][j];
                ArcCost[i][counter] = ArcCost[i][j];
                InfFlag[i][counter] = InfFlag[i][j];
                Control[i][counter] = ControlVec[j];

                // Counting the number of feasible inputs
                counter++;
            }
        }
        // Store the number of feasible control inputs per starting state
        FeasibleCounter[i] = counter;

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

#ifndef ADAPTIVEGRID

    // Interpolate the cost-to-come for states that can be reached
    for (i = 0; i < Nx; i++) {

        /*--- Point to the head of each row ---*/

        // Arc Cost (after interpolation) from one state to another (point-to-point)
        real_T *p2pCost = ArcPtr->arcCost[i];

        // Arc Control values (after interpolation) from one state to another (point-to-point)
        real_T *p2pControl = ArcPtr->arcU[i];

        // Per starting state (on the grid)
        // All these values are before interpolation
        // Xnext are not necessarily on the grid, we want them to be on the grid so that we can calculate cost-to-come grid
        real_T *Xnext_real = Xnext[i];
        real_T *ArcCost_real = ArcCost[i];
        real_T *Control_real = Control[i];

        // If there is only one control input possible
        if (FeasibleCounter[i] == 1) {

            // Find the nearest possible state on the state grid (Since it is not possible to interpolate with only one state)
            uint16_t idx = (uint16_t) findNearest(StateVecCopy, Xnext_real[0], Nx);

            // Store the cost and control to come to this point
            p2pCost[idx] = ArcCost_real[0];
            p2pControl[idx] = Control_real[0];
        } else if (FeasibleCounter[i] > 1) {
            // Sort the idx in the way the Xnext_real[0] to Xnext_real[FeasibleCounter[i]] is monotonically increasing
            sortIdx(Xnext_real, idxSort, FeasibleCounter[i]);

            // Reorder the actual vector based on the calculated idxSort
            reorderVector(Xnext_real, idxSort, FeasibleCounter[i]);
            reorderVector(ArcCost_real, idxSort, FeasibleCounter[i]);
            reorderVector(Control_real, idxSort, FeasibleCounter[i]);

            // Find the range of the states (in the StateVector) that can be reached
            uint16_t idxMax = (uint16_t) findMaxLEQ(StateVecCopy, Xnext_real[FeasibleCounter[i] - 1], Nx);
            uint16_t idxMin = (uint16_t) findMinGEQ(StateVecCopy, Xnext_real[0], Nx);


            // Linear Interpolation
            if (idxMax >= idxMin) {
                // Calculate all the possible arc costs to each state (in the State Vector)
                LookupTable CostComeTable;
                lookuptable_init(&CostComeTable, Xnext_real, ArcCost_real, FeasibleCounter[i]);
                interpolation(&CostComeTable, StateVecCopy + idxMin, p2pCost + idxMin, idxMax - idxMin + 1);
                lookuptable_free(&CostComeTable);

                // Calculate all the possible arc control signals to each state (in the State Vector)
                LookupTable ControlComeTable;
                lookuptable_init(&ControlComeTable, Xnext_real, Control_real, FeasibleCounter[i]);
                interpolation(&ControlComeTable, StateVecCopy + idxMin, p2pControl + idxMin, idxMax - idxMin + 1);
                lookuptable_free(&ControlComeTable);
            }

        }
    }
#endif

    free(Xnext);
    free(ArcCost);
    free(InfFlag);
    free(Control);
    free(FeasibleCounter);
    free(idxSort);
    free(StateVecCopy);
#ifdef ADAPTIVEGRID
    free(BoxEdgesCopy);
#endif
}

static void findSolution(SolverOutput *OutputPtr, Solution *SolutionPtr, real_T Xfmin, real_T Xfmax) {

    // Find the index range within the final state constraints
    uint16_t minIdx = (uint16_t) findMinGEQ(StateVec, Xfmin, Nx);
    uint16_t maxIdx = (uint16_t) findMaxLEQ(StateVec, Xfmax, Nx);

    real_T minCost = SolverInputPtr->SolverLimit.infValue;

    uint16_t i;
    int16_t k;
    uint16_t finalIdx;

    // Compare all the possible final cost-to-come values, take the minimum one as the minimum total cost
    for (i = minIdx; i <= maxIdx; i++) {
        if (SolutionPtr->CostToCome[i] < minCost) {
            minCost = SolutionPtr->CostToCome[i];
            finalIdx = i;
        }
    }

    //Minimum total cost
    OutputPtr->Cost = minCost;

    // Memory for the optimal state trajectory
    uint16_t *optimalstateIdx = (uint16_t *) calloc(Nhrz + 1, sizeof(uint16_t));

    // Load the final state index
    optimalstateIdx[Nhrz] = finalIdx;

    // Make sure there is at least one solution
    if (minCost < SolverInputPtr->SolverLimit.infValue) {
        // Retrieve the optimal index trajectory
        k = Nhrz - 1;

        while (k >= 0) {
            optimalstateIdx[k] = SolutionPtr->Xn[k][optimalstateIdx[k + 1]];
            k--;
        }

        // Find the optimal control policy and state trajectory
        for (i = 0; i < Nhrz; i++) {
#ifdef ADAPTIVEGRID
            OutputPtr->Vo[i] = StateGrid[i + 1][optimalstateIdx[i + 1]];
#elif defined(BOUNDCALIBRATION)
            if (optimalstateIdx[i + 1] == BoundaryStruct.boundMemo[i][2]) {
                OutputPtr->Vo[i] = BoundaryStruct.boundMemo[i][0];
                printf("Hit!!!!\n");
                printf("%d\n", i);
            } else if (optimalstateIdx[i + 1] == BoundaryStruct.boundMemo[i][3]) {
                OutputPtr->Vo[i] = BoundaryStruct.boundMemo[i][1];
                printf("Hit!!!!\n");
                printf("%d\n", i);
            } else {
                OutputPtr->Vo[i] = StateVec[optimalstateIdx[i + 1]];
            }
#else
            OutputPtr->Vo[i] = StateVec[optimalstateIdx[i + 1]];
#endif
            OutputPtr->Fo[i] = SolutionPtr->Un[i][optimalstateIdx[i + 1]];
        }
    }

    printf("%f\n", OutputPtr->Cost);
    // Free the memory
    free(optimalstateIdx);
}