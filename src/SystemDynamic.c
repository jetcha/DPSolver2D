#include "../inc/SystemDynamic.h"
#include "../inc/MathFunction.h"

/*--- Static Variables ---*/
static const DynParameter *ModelParameter;
static EnvFactor *EnvironmentalFactor;
static SolverInput *SolverInputPtr;

#ifdef DYNCOUNTER
uint32_t counterDynamics = 0;
#endif // DYNCOUNTER

/*--- Public Function Definition ---*/
void PassParameters(SolverInput *InputPtr, DynParameter *ModelParaPtr, EnvFactor *EnvFactorPtr) {
    ModelParameter = ModelParaPtr;
    EnvironmentalFactor = EnvFactorPtr;
    SolverInputPtr = InputPtr;
}

void createLineSpace(real_T *Vector, real_T min, real_T max, uint32_t N) {
    real_T delta = (max - min) / (N - 1);
    uint32_t i;
    Vector[0] = min;

    for (i = 1; i < N; i++) {
        Vector[i] = Vector[i - 1] + delta;
    }
}

/*
void
speedDynamics(uint16_t Nx, uint16_t Nu, real_T (*Xnext)[Nu], real_T (*ArcCost)[Nu], uint8_t (*InfFlag)[Nu],
              real_T const *StateVec, real_T const *ControlVec, Boundary *BoundaryPtr, uint16_t N, uint16_t X0_index) {

    // Local Copy the State and Control grids
    real_T *Xin = (real_T *) malloc(Nx * sizeof(real_T));
    real_T *Uin = (real_T *) malloc(Nu * sizeof(real_T));
    memcpy(Xin, StateVec, Nx * sizeof(real_T));
    memcpy(Uin, ControlVec, Nu * sizeof(real_T));

    // Environmental Factors
#if defined(NORMALBOUND) || defined(CUSTOMBOUND)
    // Copy the boundary lines
    real_T Vmax_start = BoundaryPtr->upperBound[N];
    real_T Vmin_start = BoundaryPtr->lowerBound[N];
    real_T Vmax_end = BoundaryPtr->upperBound[N + 1];
    real_T Vmin_end = BoundaryPtr->lowerBound[N + 1];
#elif defined(NOBOUND)
    real_T Vmax_end = EnvironmentalFactor->Vmax_env[N + 1];
    real_T Vmin_end = EnvironmentalFactor->Vmin_env[N + 1];
#endif

#ifdef BOUNDCALIBRATION
    // Besides the points on the state grid, also consider the points on the boundary
    if (N > 0) {
        uint16_t minIdx = (uint16_t) findMaxLEQ(Xin, BoundaryPtr->boundMemo[N - 1][0], Nx);
        uint16_t maxIdx = (uint16_t) findMinGEQ(Xin, BoundaryPtr->boundMemo[N - 1][1], Nx);
        Xin[minIdx] = BoundaryPtr->boundMemo[N - 1][0];
        Xin[maxIdx] = BoundaryPtr->boundMemo[N - 1][1];
    }
#endif

    // Read the slope angle
    real_T angle = EnvironmentalFactor->Angle_env[N + 1];

    // Intermediate Variables
    real_T Pwh;
    real_T Pm;
    real_T Pinv;
    real_T Pdc;
    real_T Pbatt;
    real_T dt;

    // Hard Constraints
    real_T PDmax = SolverInputPtr->Constraint.PDmax;
    real_T PAmax = SolverInputPtr->Constraint.PAmax;

    // Parameters
    real_T m = ModelParameter->m;
    real_T g = ModelParameter->g;
    real_T crr = ModelParameter->crr;
    real_T CdA = ModelParameter->CdA;
    real_T ds = ModelParameter->ds;

    real_T penalty = ModelParameter->speedPenalty;

    real_T eta_trans = ModelParameter->eta_trans;
    real_T eta_dc = ModelParameter->eta_dc;
    real_T alpha0 = ModelParameter->alpha0;
    real_T alpha1 = ModelParameter->alpha1;
    real_T alpha2 = ModelParameter->alpha2;
    real_T beta0 = ModelParameter->beta0;

    // Calculations (all the possibilities)
    uint16_t i;
    uint16_t j;

    // Preserve the initial state accuracy
    if (N == 0) {
        Xin[X0_index] = Xinitial;
    }

    for (i = 0; i < Nx; i++) {

        for (j = 0; j < Nu; j++) {
#if defined(NORMALBOUND) || defined(CUSTOMBOUND)
            // Only calculate the states within the boundaries
            if (Xin[i] > Vmax_start || Xin[i] > SolverInputPtr->Constraint.Vmax || Xin[i] < Vmin_start ||
                Xin[i] < SolverInputPtr->Constraint.Vmin) {
                InfFlag[i][j] = 1;
                continue;
            }
#endif

            Pwh = Uin[j] * Xin[i];

#ifdef DYNCOUNTER
            counterDynamics++;
#endif // DYNCOUNTER


            // First determine if the Pwh has already exceeded the limit
            // Acceleration
            if (Pwh > 0) {
                if (Pwh > PAmax)                    // if the wheel power exceeds the limit...
                {
                    InfFlag[i][j] = 1;                // mark it as 'infeasible'
                    continue;
                }

                Pm = Pwh / eta_trans;
                Pinv = ((1 - alpha1) - sqrt((alpha1 - 1) * (alpha1 - 1) - 4 * alpha2 * (alpha0 + Pm))) / (2 * alpha2);
                Pdc = Pinv / eta_dc;
                Pbatt = (1 - sqrt(1 - 4 * beta0 * Pdc)) / (2 * beta0);

            }
                // Deceleration
            else {
                if (Pwh < PDmax) {
                    InfFlag[i][j] = 1;
                    continue;
                }
                Pm = Pwh * eta_trans;
                Pinv = ((1 - alpha1) - sqrt((alpha1 - 1) * (alpha1 - 1) - 4 * alpha2 * (alpha0 + Pm))) / (2 * alpha2);
                Pdc = Pinv * eta_dc;
                Pbatt = (1 - sqrt(1 - 4 * beta0 * Pdc)) / (2 * beta0);

            }

            // Calculate the speed at the next step -> Check the derivation in the paper
            real_T X_squared = (2 * ds / m) * Uin[j] + (1 - 2 * ds * CdA / m) * Xin[i] * Xin[i] -
                               2 * ds * g * (sin(angle) + crr * cos(angle));


            // Check if the speed squared becomes smaller than 0
            if (X_squared < 0) {
                InfFlag[i][j] = 1;
                continue;
            }

            // The speed at the next step
            Xnext[i][j] = sqrt(X_squared);


            // Check if the speed result is inside the legal speed range and physical speed limits
            if (Xnext[i][j] > Vmax_end || Xnext[i][j] > SolverInputPtr->Constraint.Vmax || Xnext[i][j] < Vmin_end ||
                Xnext[i][j] < SolverInputPtr->Constraint.Vmin) {
                InfFlag[i][j] = 1;
                continue;
            }

            // Calculate dt
            dt = 2 * ds / (Xnext[i][j] + Xin[i]);

            // Arc Cost of this combination: Xin[i] and Uin[j]
            ArcCost[i][j] = (Pbatt + penalty) * dt;

        }
    }

    // Free the memory of local state and control vectors
    free(Xin);
    free(Uin);
}

void thermalDynamics(uint16_t Nx, uint16_t Nu, real_T (*Xnext)[Nu], real_T (*ArcCost)[Nu], uint8_t (*InfFlag)[Nu],
                     real_T const *StateVec, real_T const *ControlVec, Boundary *BoundaryPtr, Bridge *BridgePtr,
                     uint16_t N, uint16_t X0_index) {

    // The real index in the horizon
    uint16_t start = N * HORIZON / RES_THERMAL;
    uint16_t end = (N + 1) * HORIZON / RES_THERMAL;

    // Local Copy the State and Control grids
    real_T *Xin = (real_T *) malloc(Nx * sizeof(real_T));
    real_T *Uin = (real_T *) malloc(Nu * sizeof(real_T));
    memcpy(Xin, StateVec, Nx * sizeof(real_T));
    memcpy(Uin, ControlVec, Nu * sizeof(real_T));

    // Required temperature by the driver
#if defined(NORMALBOUND) || defined(CUSTOMBOUND)
    // Copy the boundary lines
    real_T Tmax_start = BoundaryPtr->upperBound[N];
    real_T Tmin_start = BoundaryPtr->lowerBound[N];
    real_T Tmax_end = BoundaryPtr->upperBound[N+1];
    real_T Tmin_end = BoundaryPtr->lowerBound[N+1];
#elif defined(NOBOUND)
    real_T Tmax_end = EnvironmentalFactor->T_required[end] + 5;
    real_T Tmin_end = EnvironmentalFactor->T_required[end] - 5;
#endif

#ifdef BOUNDCALIBRATION
    // Besides the points on the state grid, also consider the points on the boundary
    if (N > 0) {
        uint16_t minIdx = (uint16_t) findMaxLEQ(Xin, BoundaryPtr->boundMemo[N - 1][0], Nx);
        uint16_t maxIdx = (uint16_t) findMinGEQ(Xin, BoundaryPtr->boundMemo[N - 1][1], Nx);
        Xin[minIdx] = BoundaryPtr->boundMemo[N - 1][0];
        Xin[maxIdx] = BoundaryPtr->boundMemo[N - 1][1];
    }
#endif

    // Intermediate Variables
    real_T Pdc = BridgePtr->Pdc[start];
    real_T Phvac;
    real_T Ps;
    real_T Pbatt;
    real_T Tinlet;

    // Hard Constraints
    real_T PACmax = SolverInputPtr->Constraint.PACmax;
    real_T Tmax_inlet = SolverInputPtr->Constraint.Tmax_inlet;
    real_T Tmin_inlet = SolverInputPtr->Constraint.Tmin_inlet;

    //Parameters
    real_T Cth = ModelParameter->Cth;
    real_T Rth = ModelParameter->Rth;
    real_T Qsun = ModelParameter->Qsun;
    real_T Qpas = ModelParameter->Qpas;
    real_T Cp = ModelParameter->Cp;
    real_T rho = ModelParameter->rho;
    real_T mDot = ModelParameter->mDot;
    real_T CoP_pos = ModelParameter->CoP_pos;
    real_T CoP_neg = ModelParameter->CoP_neg;
    real_T Tamb = ModelParameter->Tamb;

    real_T beta0 = ModelParameter->beta0;

    real_T T_required = EnvironmentalFactor->T_required[end];
    real_T penalty = ModelParameter->thermalPenalty;

    uint16_t i;
    uint16_t j;

    // To calculate delta t between a step
    real_T tDelta = 0;
    for (i = start; i < end; i++) {
        tDelta += BridgePtr->tDelta[i];
    }

    printf("Duration: %f\n", tDelta);

    // Preserve the initial state accuracy
    if (N == 0) {
        Xin[X0_index] = Xinitial;
    }

    for (i = 0; i < Nx; i++) {

        for (j = 0; j < Nu; j++) {

#if defined(NORMALBOUND) || defined(CUSTOMBOUND)
            // Only calculate the states with in the boundary
            if (Xin[i] > Tmax_start || Xin[i] > SolverInputPtr->Constraint.Tmax || Xin[i] < Tmin_start ||
                Xin[i] < SolverInputPtr->Constraint.Tmin) {
                InfFlag[i][j] = 1;
                continue;
            }
#endif

            // Heating or Cooling
            if (Uin[j] > 0) {
                Phvac = Uin[j] / CoP_pos;
            } else if (Uin[j] == 0) {
                Phvac = 0;
            } else {
                Phvac = Uin[j] / CoP_neg;
            }

#ifdef DYNCOUNTER
            counterDynamics++;
#endif // DYNCOUNTER

            // Check if Phvac exceeds the limits
            if (Phvac > PACmax) {
                InfFlag[i][j] = 1;
                continue;
            }

            Ps = Pdc + Phvac;
            Pbatt = (1 - sqrt(1 - 4 * beta0 * Ps)) / (2 * beta0);

            // Check if Tinlet exceeds the limits
            Tinlet = Xin[i] + Uin[j] / (Cp * rho * mDot);
            if (Tinlet > Tmax_inlet || Tinlet < Tmin_inlet) {
                InfFlag[i][j] = 1;
                continue;
            }

            // Temperature at the next step
            Xnext[i][j] = Xin[i] + (tDelta / Cth) * (Uin[j] + Qsun + Qpas + (Tamb - Xin[i]) / Rth);

            // Check if it stays in the cabin temperature limit
            if (Xnext[i][j] > Tmax_end || Xnext[i][j] < Tmin_end) {
                InfFlag[i][j] = 1;
                continue;
            }

            // ArcCost - added a L2 norm as the penalization
            ArcCost[i][j] = Pbatt * tDelta + penalty * (Xnext[i][j] - T_required) * (Xnext[i][j] - T_required);
        }
    }

    // Free the memory of local state and control vectors
    free(Xin);
    free(Uin);
}

void initBridge(Bridge *BridgePtr) {
    BridgePtr->Pdc = (real_T *) malloc(HORIZON * sizeof(real_T));
    BridgePtr->tDelta = (real_T *) malloc(HORIZON * sizeof(real_T));
}

void freeBridge(Bridge *BridgePtr) {
    free(BridgePtr->Pdc);
    free(BridgePtr->tDelta);
}

void bridgeConnection(Bridge *BridgePtr, SolverOutput *OutputPtr, real_T V0) {

    // Intermediate Variables
    real_T Pwh;
    real_T Pm;
    real_T Pinv;

    // Parameters
    real_T ds = ModelParameter->ds;

    real_T eta_trans = ModelParameter->eta_trans;
    real_T eta_dc = ModelParameter->eta_dc;
    real_T alpha0 = ModelParameter->alpha0;
    real_T alpha1 = ModelParameter->alpha1;
    real_T alpha2 = ModelParameter->alpha2;
    real_T beta0 = ModelParameter->beta0;

    uint16_t i;

    real_T duration = 0;
    for (i = 0; i < HORIZON; i++) {

        // Calculate delta t
        // Whether it's the initial state or not
        if (i == 0) {
            BridgePtr->tDelta[i] = 2 * ds / (OutputPtr->Vo[i] + V0);
            Pwh = OutputPtr->Fo[i] * V0;
        } else {
            BridgePtr->tDelta[i] = 2 * ds / (OutputPtr->Vo[i] + OutputPtr->Vo[i - 1]);
            Pwh = OutputPtr->Fo[i] * OutputPtr->Vo[i - 1];
        }

        duration += BridgePtr->tDelta[i];

        // Calculate Pdc
        // Acceleration
        if (Pwh > 0) {
            Pm = Pwh / eta_trans;
            Pinv = ((1 - alpha1) - sqrt((alpha1 - 1) * (alpha1 - 1) - 4 * alpha2 * (alpha0 + Pm))) / (2 * alpha2);
            BridgePtr->Pdc[i] = Pinv / eta_dc;
        }
            //Deceleration
        else {
            Pm = Pwh * eta_trans;
            Pinv = ((1 - alpha1) - sqrt((alpha1 - 1) * (alpha1 - 1) - 4 * alpha2 * (alpha0 + Pm))) / (2 * alpha2);
            BridgePtr->Pdc[i] = Pinv * eta_dc;
        }
    }

    real_T Pbatt;
    real_T Edem = 0;

    for (i = 0; i < HORIZON; i++) {
        Pbatt = (1 - sqrt(1 - 4 * beta0 * BridgePtr->Pdc[i])) / (2 * beta0);
        Edem += Pbatt * BridgePtr->tDelta[i];
    }

    printf("\nPower demanding: %f\n", Edem);
    printf("\nTotal duration:%f\n\n", duration);

}
*/

void systemDynamics(StateTuple (*Xnext)[NT][NF][NQ], real_T (*ArcCost)[NT][NF][NQ], uint8_t (*InfFlag)[NT][NF][NQ],
                    real_T const *SpeedVec, real_T const *ForceVec, real_T const *TempVec, real_T const *InletVec,
                    uint16_t V0_index, uint16_t T0_index, Boundary *BoundaryPtr, uint16_t N) {

    // Grid sizes
    uint16_t Nv = SolverInputPtr->GridSize.Nv;
    uint16_t Nf = SolverInputPtr->GridSize.Nf;
    uint16_t Nt = SolverInputPtr->GridSize.Nt;
    uint16_t Nq = SolverInputPtr->GridSize.Nq;

    // Local Copy of the state and control vectors
    real_T *Vin = (real_T *) malloc(Nv * sizeof(real_T));
    real_T *Fin = (real_T *) malloc(Nf * sizeof(real_T));
    real_T *Tin = (real_T *) malloc(Nt * sizeof(real_T));
    real_T *Qin = (real_T *) malloc(Nq * sizeof(real_T));

    memcpy(Vin, SpeedVec, Nv * sizeof(real_T));
    memcpy(Fin, ForceVec, Nf * sizeof(real_T));
    memcpy(Tin, TempVec, Nt * sizeof(real_T));
    memcpy(Qin, InletVec, Nq * sizeof(real_T));

    // Environmental Factors
#if defined(CUSTOMBOUND)
    // Copy the boundary lines
    real_T Vmax_start = BoundaryPtr->upperBound[N];
    real_T Vmin_start = BoundaryPtr->lowerBound[N];
    real_T Vmax_end = BoundaryPtr->upperBound[N + 1];
    real_T Vmin_end = BoundaryPtr->lowerBound[N + 1];
#elif defined(NOBOUND)
    real_T Vmax_end = EnvironmentalFactor->Vmax_env[N + 1];
    real_T Vmin_end = EnvironmentalFactor->Vmin_env[N + 1];
#endif

#ifdef BOUNDCALIBRATION
    // Besides the points on the state grid, also consider the points on the boundary
    if (N > 0) {
        uint16_t minIdx = (uint16_t) findMaxLEQ(Vin, BoundaryPtr->boundMemo[N - 1][0], Nv);
        uint16_t maxIdx = (uint16_t) findMinGEQ(Vin, BoundaryPtr->boundMemo[N - 1][1], Nv);
        Vin[minIdx] = BoundaryPtr->boundMemo[N - 1][0];
        Vin[maxIdx] = BoundaryPtr->boundMemo[N - 1][1];
    }
#endif

    real_T Tmax_end = EnvironmentalFactor->T_required[N] + 5;
    real_T Tmin_end = EnvironmentalFactor->T_required[N] - 5;

    real_T angle = EnvironmentalFactor->Angle_env[N + 1];

    // Intermediate Variables
    real_T Pwh;
    real_T Pm;
    real_T Pinv;
    real_T Pdc;
    real_T Ps;
    real_T Pbatt;
    real_T Phvac;

    real_T dt;

    real_T Tinlet;

    // Hard Constraints
    real_T PDmax = SolverInputPtr->Constraint.PDmax;
    real_T PAmax = SolverInputPtr->Constraint.PAmax;

    real_T PACmax = SolverInputPtr->Constraint.PACmax;
    real_T Tmax_inlet = SolverInputPtr->Constraint.Tmax_inlet;
    real_T Tmin_inlet = SolverInputPtr->Constraint.Tmin_inlet;

    // Parameters
    real_T m = ModelParameter->m;
    real_T g = ModelParameter->g;
    real_T crr = ModelParameter->crr;
    real_T CdA = ModelParameter->CdA;

    real_T ds = ModelParameter->ds;
    real_T speedPenalty = ModelParameter->speedPenalty;

    real_T eta_trans = ModelParameter->eta_trans;
    real_T eta_dc = ModelParameter->eta_dc;
    real_T alpha0 = ModelParameter->alpha0;
    real_T alpha1 = ModelParameter->alpha1;
    real_T alpha2 = ModelParameter->alpha2;
    real_T beta0 = ModelParameter->beta0;

    real_T Cth = ModelParameter->Cth;
    real_T Rth = ModelParameter->Rth;
    real_T Qsun = ModelParameter->Qsun;
    real_T Qpas = ModelParameter->Qpas;
    real_T Cp = ModelParameter->Cp;
    real_T rho = ModelParameter->rho;
    real_T mDot = ModelParameter->mDot;
    real_T CoP_pos = ModelParameter->CoP_pos;
    real_T CoP_neg = ModelParameter->CoP_neg;
    real_T Tamb = ModelParameter->Tamb;

    real_T T_required = EnvironmentalFactor->T_required[N];
    real_T thermalPenalty = ModelParameter->thermalPenalty;

    // Preserve the initial state accuracy
    if (N == 0) {
        Vin[V0_index] = Vinitial;
        Tin[T0_index] = Tinitial;
    }

    // Calculations (all the possibilities)
    uint16_t i, j, k, l;

    for (i = 0; i < Nv; i++) {
        for (j = 0; j < Nt; j++) {

            // Try all the possible force inputs
            for (k = 0; k < Nf; k++) {

#if defined(CUSTOMBOUND)
                // Only calculate the states within the boundaries
                if (Vin[i] > Vmax_start || Vin[i] > SolverInputPtr->Constraint.Vmax || Vin[i] < Vmin_start ||
                    Vin[i] < SolverInputPtr->Constraint.Vmin) {
                    continue;
                }
#endif

                /// Powertrain ///
                // Wheel power
                Pwh = Vin[i] * Fin[k];

#ifdef DYNCOUNTER
                counterDynamics++;
#endif // DYNCOUNTER

                // First determine if the Pwh has already exceeded the limit
                // Acceleration
                if (Pwh > 0) {
                    // if the wheel power exceeds the limit, leave it as 'infeasible'
                    if (Pwh > PAmax) {
                        continue;
                    }

                    Pm = Pwh / eta_trans;
                    Pinv = ((1 - alpha1) - sqrt((alpha1 - 1) * (alpha1 - 1) - 4 * alpha2 * (alpha0 + Pm))) /
                           (2 * alpha2);
                    Pdc = Pinv / eta_dc;
                }
                    // Deceleration
                else {
                    if (Pwh < PDmax) {
                        continue;
                    }
                    Pm = Pwh * eta_trans;
                    Pinv = ((1 - alpha1) - sqrt((alpha1 - 1) * (alpha1 - 1) - 4 * alpha2 * (alpha0 + Pm))) /
                           (2 * alpha2);
                    Pdc = Pinv * eta_dc;
                }

                /// Speed Dynamics ///
                // Calculate the speed dynamics
                real_T V_squared = (2 * ds / m) * Fin[k] + (1 - 2 * ds * CdA / m) * Vin[i] * Vin[i] -
                                   2 * ds * g * (sin(angle) + crr * cos(angle));

                // Check if the speed squared becomes smaller than 0
                if (V_squared < 0) {
                    continue;
                }

                // Speed at the next step
                real_T Vnext = sqrt(V_squared);

                // Check if the speed result is inside the legal speed range and physical speed limits
                if (Vnext > Vmax_end || Vnext > SolverInputPtr->Constraint.Vmax || Vnext < Vmin_end ||
                    Vnext < SolverInputPtr->Constraint.Vmin) {
                    continue;
                }

                // Try all the possible heat inlet inputs
                for (l = 0; l < Nq; l++) {

                    // HVAC power
                    // Heating or Cooling
                    if (Qin[l] > 0) {
                        Phvac = Qin[l] / CoP_pos;
                    } else if (Qin[l] == 0) {
                        Phvac = 0;
                    } else {
                        Phvac = Qin[l] / CoP_neg;
                    }

#ifdef DYNCOUNTER
                    counterDynamics++;
#endif // DYNCOUNTER

                    // Check if Phvac exceeds the limits
                    if (Phvac > PACmax) {
                        continue;
                    }

                    // Check if Tinlet exceeds the limits
                    Tinlet = Tin[j] + Qin[l] / (Cp * rho * mDot);
                    if (Tinlet > Tmax_inlet || Tinlet < Tmin_inlet) {
                        continue;
                    }

                    // Calculate the power of the battery
                    Ps = Pdc + Phvac;
                    Pbatt = (1 - sqrt(1 - 4 * beta0 * Ps)) / (2 * beta0);

                    // Speed at the next step
                    Xnext[i][j][k][l].V = Vnext;

                    // Calculate dt
                    dt = 2 * ds / (Xnext[i][j][k][l].V + Vin[i]);

                    /// Thermal Dynamics ///
                    // Temperature at the next step
                    Xnext[i][j][k][l].T = Tin[j] + (dt / Cth) * (Qin[l] + Qsun + Qpas + (Tamb - Tin[j]) / Rth);

                    // Check if it stays in the cabin temperature limit
                    if (Xnext[i][j][k][l].T > Tmax_end || Xnext[i][j][k][l].T < Tmin_end) {
                        continue;
                    }

                    /// Arc Cost ///
                    // ArcCost - added a L2 norm as the penalization
                    ArcCost[i][j][k][l] = (Pbatt + speedPenalty) * dt +
                                          thermalPenalty * (Xnext[i][j][k][l].T - T_required) *
                                          (Xnext[i][j][k][l].T - T_required);

                    // When all the constraints are required, mark it as feasible
                    InfFlag[i][j][k][l] = 0;
                }
            }
        }
    }

    // Free the local copies
    free(Vin);
    free(Tin);
    free(Fin);
    free(Qin);
}