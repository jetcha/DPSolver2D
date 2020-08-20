#include "../inc/SystemDynamic.h"
#include "../inc/MathFunction.h"

/*--- Static Variables ---*/
static const DynParameter *ModelParameter;
static EnvFactor *EnvironmentalFactor;
static SolverInput *SolverInputPtr;

/*--- Global Variables ---*/
#ifdef ADAPTIVEGRID
real_T (*AdaptiveSpeedPlane)[NT];               // Adaptive Grid, Used to store the starting point
real_T (*AdaptiveTempPlane)[NT];
#endif

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

void createStatePlane(StateTuple (*Plane)[NT], const real_T *VectorV, const real_T *VectorT, uint16_t Nv, uint16_t Nt){
    uint16_t i, j;

    for(i = 0; i < Nv; i++){
        for(j = 0; j < Nt; j++){
            Plane[i][j].V = VectorV[i];
            Plane[i][j].T = VectorT[j];
        }
    }
}

#ifdef ADAPTIVEGRID
void passStatePlane(StateTuple (*Plane)[NT], uint16_t Nv, uint16_t Nt){
    uint16_t i, j;

    AdaptiveSpeedPlane = malloc(sizeof(real_T[Nv][Nt]));
    AdaptiveTempPlane = malloc(sizeof(real_T[Nv][Nt]));

    for(i = 0; i < Nv; i++){
        for(j = 0; j < Nt; j++) {
            AdaptiveSpeedPlane[i][j] = Plane[i][j].V;
            AdaptiveTempPlane[i][j] = Plane[i][j].T;
        }
    }
}
#endif

void systemDynamics(StateTuple (*Xnext)[NT][NF][NQ], real_T (*ArcCost)[NT][NF][NQ], uint8_t (*InfFlag)[NT][NF][NQ],
                    real_T const *SpeedVec, real_T const *ForceVec, real_T const *TempVec, real_T const *InletVec,
                    uint16_t V0_index, uint16_t T0_index, Boundary *BoundaryPtr, uint16_t N) {

    // Grid sizes
    uint16_t Nv = SolverInputPtr->GridSize.Nv;
    uint16_t Nf = SolverInputPtr->GridSize.Nf;
    uint16_t Nt = SolverInputPtr->GridSize.Nt;
    uint16_t Nq = SolverInputPtr->GridSize.Nq;

    // Local Copy of the state and control vectors
#ifdef ADAPTIVEGRID
    real_T (*Vin)[Nt] = malloc(sizeof(real_T[Nv][Nt]));
    real_T (*Tin)[Nt] = malloc(sizeof(real_T[Nv][Nt]));
    memcpy(Vin, AdaptiveSpeedPlane, Nv * Nt * sizeof(real_T));
    memcpy(Tin, AdaptiveTempPlane, Nv * Nt * sizeof(real_T));
    free(AdaptiveSpeedPlane);
    free(AdaptiveTempPlane);
#else
    real_T *Vin = (real_T *) malloc(Nv * sizeof(real_T));
    real_T *Tin = (real_T *) malloc(Nt * sizeof(real_T));
    memcpy(Vin, SpeedVec, Nv * sizeof(real_T));
    memcpy(Tin, TempVec, Nt * sizeof(real_T));
#endif

    real_T *Fin = (real_T *) malloc(Nf * sizeof(real_T));
    real_T *Qin = (real_T *) malloc(Nq * sizeof(real_T));
    memcpy(Fin, ForceVec, Nf * sizeof(real_T));
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
#ifdef ADAPTIVEGRID
        Vin[V0_index][T0_index] = Vinitial;
        Tin[V0_index][T0_index] = Tinitial;
#else
        Vin[V0_index] = Vinitial;
        Tin[T0_index] = Tinitial;
#endif
    }

    // Calculations (all the possibilities)
    uint16_t i, j, k, l;

    for (i = 0; i < Nv; i++) {

#if defined(CUSTOMBOUND)
        // Only start from the given initial state
        if(N == 0){
            if(i != V0_index) continue;
        }
#endif


        for (j = 0; j < Nt; j++) {

#if defined(CUSTOMBOUND)
            // Only start from the given initial state
            if(N == 0){
                if(j != T0_index) continue;
            }
#endif

            // Try all the possible force inputs
            for (k = 0; k < Nf; k++) {

#if defined(CUSTOMBOUND)
#ifdef ADAPTIVEGRID
                // Only calculate the states within the boundaries
                if (Vin[i][j] > Vmax_start || Vin[i][j] > SolverInputPtr->Constraint.Vmax || Vin[i][j] < Vmin_start ||
                    Vin[i][j] < SolverInputPtr->Constraint.Vmin) {
                    continue;
                }
#else
                // Only calculate the states within the boundaries
                if (Vin[i] > Vmax_start || Vin[i] > SolverInputPtr->Constraint.Vmax || Vin[i] < Vmin_start ||
                    Vin[i] < SolverInputPtr->Constraint.Vmin) {
                    continue;
                }
#endif
#endif

                /// Powertrain ///
                // Wheel power
#ifdef ADAPTIVEGRID
                Pwh = Vin[i][j] * Fin[k];
#else
                Pwh = Vin[i] * Fin[k];
#endif

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
#ifdef ADAPTIVEGRID
                real_T V_squared = (2 * ds / m) * Fin[k] + (1 - 2 * ds * CdA / m) * Vin[i][j] * Vin[i][j] -
                                   2 * ds * g * (sin(angle) + crr * cos(angle));
#else
                real_T V_squared = (2 * ds / m) * Fin[k] + (1 - 2 * ds * CdA / m) * Vin[i] * Vin[i] -
                                   2 * ds * g * (sin(angle) + crr * cos(angle));
#endif


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
#ifdef ADAPTIVEGRID
                    Tinlet = Tin[i][j] + Qin[l] / (Cp * rho * mDot);
#else
                    Tinlet = Tin[j] + Qin[l] / (Cp * rho * mDot);
#endif

                    if (Tinlet > Tmax_inlet || Tinlet < Tmin_inlet) {
                        continue;
                    }

                    // Calculate the power of the battery
                    Ps = Pdc + Phvac;
                    Pbatt = (1 - sqrt(1 - 4 * beta0 * Ps)) / (2 * beta0);

                    // Speed at the next step
                    Xnext[i][j][k][l].V = Vnext;

                    // Calculate dt
#ifdef ADAPTIVEGRID
                    dt = 2 * ds / (Xnext[i][j][k][l].V + Vin[i][j]);
#else
                    dt = 2 * ds / (Xnext[i][j][k][l].V + Vin[i]);
#endif

                    /// Thermal Dynamics ///
                    // Temperature at the next step
#ifdef ADAPTIVEGRID
                    Xnext[i][j][k][l].T = Tin[i][j] + (dt / Cth) * (Qin[l] + Qsun + Qpas + (Tamb - Tin[i][j]) / Rth);
#else
                    Xnext[i][j][k][l].T = Tin[j] + (dt / Cth) * (Qin[l] + Qsun + Qpas + (Tamb - Tin[j]) / Rth);
#endif


                    // Check if it stays in the cabin temperature limit
                    if (Xnext[i][j][k][l].T > Tmax_end || Xnext[i][j][k][l].T < Tmin_end) {
                        continue;
                    }

                    /// Arc Cost ///
                    // ArcCost - added a L2 norm as the penalization
                    ArcCost[i][j][k][l] = (0*Pbatt + speedPenalty) * dt +
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