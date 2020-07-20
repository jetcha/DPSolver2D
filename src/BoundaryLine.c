#include "../inc/BoundaryLine.h"

#ifdef BOUNDCOUNTER
uint32_t counterBound = 0;
#endif // BOUNDCOUNTER

#if defined(NORMALBOUND) || defined(CUSTOMBOUND)
void initSpeedBoundary(Boundary *BoundaryPtr) {
    BoundaryPtr->lowerBound = (real_T *) malloc((HORIZON + 1) * sizeof(real_T));
    BoundaryPtr->upperBound = (real_T *) malloc((HORIZON + 1) * sizeof(real_T));
#ifdef BOUNDCALIBRATION
    BoundaryPtr->boundMemo = malloc(sizeof(real_T[HORIZON][4]));
#endif
}

void initThermalBoundary(Boundary *BoundaryPtr) {
    BoundaryPtr->lowerBound = (real_T *) malloc((RES_THERMAL + 1) * sizeof(real_T));
    BoundaryPtr->upperBound = (real_T *) malloc((RES_THERMAL + 1) * sizeof(real_T));
#ifdef BOUNDCALIBRATION
    BoundaryPtr->boundMemo = malloc(sizeof(real_T[RES_THERMAL][4]));
#endif
}

void copySpeedBoundary(Boundary *BoundaryPtr, SolverOutput *OutputPtr) {
    memcpy(OutputPtr->lowerSpeedBound, BoundaryPtr->lowerBound, (HORIZON + 1) * sizeof(real_T));
    memcpy(OutputPtr->upperSpeedBound, BoundaryPtr->upperBound, (HORIZON + 1) * sizeof(real_T));
#ifdef BOUNDCALIBRATION
    for (uint16_t i = 0; i < HORIZON; i++) {
        OutputPtr->lowerSpeedActual[i] = BoundaryPtr->boundMemo[i][0];
        OutputPtr->upperSpeedActual[i] = BoundaryPtr->boundMemo[i][1];
    }
#endif
}

void copyThermalBoundary(Boundary *BoundaryPtr, SolverOutput *OutputPtr) {
    memcpy(OutputPtr->lowerTempBound, BoundaryPtr->lowerBound, (RES_THERMAL + 1) * sizeof(real_T));
    memcpy(OutputPtr->upperTempBound, BoundaryPtr->upperBound, (RES_THERMAL + 1) * sizeof(real_T));
#ifdef BOUNDCALIBRATION
    for (uint16_t i = 0; i < RES_THERMAL; i++) {
        OutputPtr->lowerSpeedActual[i] = BoundaryPtr->boundMemo[i][0];
        OutputPtr->upperSpeedActual[i] = BoundaryPtr->boundMemo[i][1];
    }
#endif
}

void freeBoundary(Boundary *BoundaryPtr) {
    free(BoundaryPtr->lowerBound);
    free(BoundaryPtr->upperBound);
#ifdef BOUNDCALIBRATION
    free(BoundaryPtr->boundMemo);
#endif
}

void normalSpeedBoundary(Boundary *BoundaryPtr, EnvFactor *EnvPtr) {
    for (uint16_t i = 0; i <= HORIZON; i++) {
        BoundaryPtr->upperBound[i] = EnvPtr->Vmax_env[i];
        BoundaryPtr->lowerBound[i] = EnvPtr->Vmin_env[i];
    }
}


void customSpeedBoundary(Boundary *BoundaryPtr, SolverInput *SolverInputPtr, DynParameter *ParaPtr, EnvFactor *EnvPtr,
                         real_T X0) {
    // Iteration Index
    uint16_t i, j, k;

    // Number of range blocks
    uint16_t numBlock = sizeof(EnvPtr->endBlock) / sizeof(EnvPtr->endBlock[0]);

    // Note down the change points
    uint16_t *changePoint = (uint16_t *) malloc((numBlock + 1) * sizeof(EnvPtr->endBlock[0]));

    // For example, when endBlock = [50, 100, 150], changePoint = [0, 50, 100, 150] (add zero to the beginning)
    changePoint[0] = 0;
    for (i = 1; i <= numBlock; i++) {
        changePoint[i] = EnvPtr->endBlock[i - 1];
    }

    // Local Copy of Horizon Length
    uint16_t Nhrz = SolverInputPtr->GridSize.Nhrz;

    // Local Copy of Constraints
    real_T Fmax = SolverInputPtr->Constraint.Fmax;
    real_T Fmin = SolverInputPtr->Constraint.Fmin;
    real_T PAmax = SolverInputPtr->Constraint.PAmax;
    real_T PDmax = SolverInputPtr->Constraint.PDmax;

    // Local Copy of Parameters
    real_T m = ParaPtr->m;
    real_T g = ParaPtr->g;
    real_T crr = ParaPtr->crr;
    real_T CdA = ParaPtr->CdA;
    real_T ds = ParaPtr->ds;
    real_T penalty = ParaPtr->speedPenalty;
    real_T eta_trans = ParaPtr->eta_trans;
    real_T eta_dc = ParaPtr->eta_dc;
    real_T alpha0 = ParaPtr->alpha0;
    real_T alpha1 = ParaPtr->alpha1;
    real_T alpha2 = ParaPtr->alpha2;
    real_T beta0 = ParaPtr->beta0;

    real_T maxForce;
    real_T minForce;
    real_T upper_squared;
    real_T lower_squared;

    // First, Copy the known legal speed ranges to upperBound and lower Bound
    for (i = 0; i <= Nhrz; i++) {
        BoundaryPtr->upperBound[i] = EnvPtr->Vmax_env[i];
        BoundaryPtr->lowerBound[i] = EnvPtr->Vmin_env[i];
    }

    // Let upper and lower bound both start from the known X0
    BoundaryPtr->upperBound[0] = X0;
    BoundaryPtr->lowerBound[0] = X0;


    // If there is only one range block
    if (numBlock == 1) {
        // Draw Upper Bound
        for (i = 0; i < Nhrz; i++) {
            // Maximum possible Force
            if (Fmax > (PAmax / BoundaryPtr->upperBound[i])) {
                maxForce = PAmax / BoundaryPtr->upperBound[i];
            } else {
                maxForce = Fmax;
            }

            upper_squared = (2 * ds / m) * maxForce +
                            (1 - 2 * ds * CdA / m) * BoundaryPtr->upperBound[i] * BoundaryPtr->upperBound[i] -
                            2 * ds * g * (sin(EnvPtr->Angle_env[i]) + crr * cos(EnvPtr->Angle_env[i]));

            // If it has exceeded the legal range, break
            if ((sqrt(upper_squared)) >= BoundaryPtr->upperBound[i + 1]) {
                break;
            }
                // If not, we give it a new (lower) value to it
            else {
                BoundaryPtr->upperBound[i + 1] = sqrt(upper_squared);
            }
        }

        // Draw Lower Bound
        for (i = 0; i < Nhrz; i++) {
            // Minimum possible Force
            if (Fmin < (PDmax / BoundaryPtr->lowerBound[i])) {
                minForce = PDmax / BoundaryPtr->lowerBound[i];
            } else {
                minForce = Fmin;
            }

            lower_squared = (2 * ds / m) * minForce +
                            (1 - 2 * ds * CdA / m) * BoundaryPtr->lowerBound[i] * BoundaryPtr->lowerBound[i] -
                            2 * ds * g * (sin(EnvPtr->Angle_env[i]) + crr * cos(EnvPtr->Angle_env[i]));

            // Make sure the lower bound value is not negative
            if (lower_squared < 0) {
                break;
            }
                // If it has exceeded the legal range, break
            else if ((sqrt(lower_squared)) <= BoundaryPtr->lowerBound[i + 1]) {
                break;
            }
                // If not, we give it a new (higher) value to it
            else {
                BoundaryPtr->lowerBound[i + 1] = sqrt(lower_squared);
            }
        }
    }
        // If there are multiple range blocks
    else {
        // Draw boundary lines block by block
        for (k = 0; k < numBlock; k++) {
            // #1 Block
            if (k == 0) {
                // Draw Upper Bound
                for (i = 0; i < (changePoint[k + 1] - 1); i++) {
                    // Maximum possible Force
                    if (Fmax > (PAmax / BoundaryPtr->upperBound[i])) {
                        maxForce = PAmax / BoundaryPtr->upperBound[i];
                    } else {
                        maxForce = Fmax;
                    }

                    upper_squared = (2 * ds / m) * maxForce +
                                    (1 - 2 * ds * CdA / m) * BoundaryPtr->upperBound[i] * BoundaryPtr->upperBound[i] -
                                    2 * ds * g * (sin(EnvPtr->Angle_env[i]) + crr * cos(EnvPtr->Angle_env[i]));
                    counterBound++;

                    // Same thing as we did when numBlock == 1
                    if ((sqrt(upper_squared)) >= BoundaryPtr->upperBound[i + 1]) {
                        break;
                    } else {
                        BoundaryPtr->upperBound[i + 1] = sqrt(upper_squared);
                    }
                }

                // Draw Lower Bound
                for (i = 0; i < (changePoint[k + 1] - 1); i++) {
                    // Minimum possible Force
                    if (Fmin < (PDmax / BoundaryPtr->lowerBound[i])) {
                        minForce = PDmax / BoundaryPtr->lowerBound[i];
                    } else {
                        minForce = Fmin;
                    }

                    lower_squared = (2 * ds / m) * minForce +
                                    (1 - 2 * ds * CdA / m) * BoundaryPtr->lowerBound[i] * BoundaryPtr->lowerBound[i] -
                                    2 * ds * g * (sin(EnvPtr->Angle_env[i]) + crr * cos(EnvPtr->Angle_env[i]));
                    counterBound++;

                    // Same thing as we did when numBlock == 1
                    if (lower_squared < 0) {
                        break;
                    } else if ((sqrt(lower_squared)) <= BoundaryPtr->lowerBound[i + 1]) {
                        break;
                    } else {
                        BoundaryPtr->lowerBound[i + 1] = sqrt(lower_squared);
                    }
                }
            }
                // From the #2 Block (4 conditions to consider)
            else {
                // If Vmax[Block_k] > Vmax[Block_k-1], draw the Left-Top upper bound forwards
                if (EnvPtr->Vmax_env[changePoint[k]] > EnvPtr->Vmax_env[changePoint[k - 1]]) {
                    for (i = changePoint[k] - 1; i < (changePoint[k + 1] - 1); i++) {
                        // Maximum possible Force
                        if (Fmax > (PAmax / BoundaryPtr->upperBound[i])) {
                            maxForce = PAmax / BoundaryPtr->upperBound[i];
                        } else {
                            maxForce = Fmax;
                        }

                        upper_squared = (2 * ds / m) * maxForce + (1 - 2 * ds * CdA / m) * BoundaryPtr->upperBound[i] *
                                                                  BoundaryPtr->upperBound[i] -
                                        2 * ds * g * (sin(EnvPtr->Angle_env[i]) + crr * cos(EnvPtr->Angle_env[i]));
                        counterBound++;

                        // Again, the same exact thing we have done
                        if ((sqrt(upper_squared)) >= BoundaryPtr->upperBound[i + 1]) {
                            break;
                        } else {
                            BoundaryPtr->upperBound[i + 1] = sqrt(upper_squared);
                        }
                    }
                }
                    // If Vmax[Block_k] < Vmax[Block_k-1], draw the Right-Top upper bound (backwards)
                else {
                    for (i = changePoint[k] - 1; i >= changePoint[k - 1]; i--) {
                        // Predictive minimum force
                        minForce = Fmin;
                        //upper_squared = ((upperBound[i + 1])*(upperBound[i + 1]) - 2 * ds*minForce / m + 2 * ds*g*(sin(EnvPtr->Angle_env[i]) + crr*cos(EnvPtr->Angle_env[i]))) / (1 - 2 * ds*CdA / m);

                        do {
                            upper_squared = ((BoundaryPtr->upperBound[i + 1]) * (BoundaryPtr->upperBound[i + 1]) -
                                             2 * ds * minForce / m + 2 * ds * g * (sin(EnvPtr->Angle_env[i]) +
                                                                                   crr * cos(EnvPtr->Angle_env[i]))) /
                                            (1 - 2 * ds * CdA / m);
                            minForce = minForce * 0.95;
                            counterBound++;
                        } while ((sqrt(upper_squared)) * minForce < PDmax);

                        if ((sqrt(upper_squared)) >= BoundaryPtr->upperBound[i]) {
                            break;
                        } else {
                            BoundaryPtr->upperBound[i] = sqrt(upper_squared);
                        }
                    }
                }

                // If Vmin[Block_k] < Vmin[Block_k-1], draw the Left-Bottom lower bound (forwards)
                if (EnvPtr->Vmin_env[changePoint[k]] < EnvPtr->Vmin_env[changePoint[k - 1]]) {
                    for (i = changePoint[k] - 1; i < (changePoint[k + 1] - 1); i++) {
                        // Minimum possible Force
                        if (Fmin < (PDmax / BoundaryPtr->lowerBound[i])) {
                            minForce = PDmax / BoundaryPtr->lowerBound[i];
                        } else {
                            minForce = Fmin;
                        }

                        lower_squared = (2 * ds / m) * minForce + (1 - 2 * ds * CdA / m) * BoundaryPtr->lowerBound[i] *
                                                                  BoundaryPtr->lowerBound[i] -
                                        2 * ds * g * (sin(EnvPtr->Angle_env[i]) + crr * cos(EnvPtr->Angle_env[i]));
                        counterBound++;

                        // Again, the same exact thing
                        if (lower_squared < 0) {
                            break;
                        } else if ((sqrt(lower_squared)) <= BoundaryPtr->lowerBound[i + 1]) {
                            break;
                        } else {
                            BoundaryPtr->lowerBound[i + 1] = sqrt(lower_squared);
                        }
                    }
                }
                    // If Vmin[Block_k] > Vmin[Block_k-1], draw the Right-Bottom lower bound (backwards)
                else {
                    for (i = changePoint[k] - 1; i >= changePoint[k - 1]; i--) {
                        // Predictive maximum force
                        maxForce = Fmax;
                        //lower_squared = ((lowerBound[i + 1])*(lowerBound[i + 1]) - 2 * ds*maxForce / m + 2 * ds*g*(sin(EnvPtr->Angle_env[i]) + crr*cos(EnvPtr->Angle_env[i]))) / (1 - 2 * ds*CdA / m);

                        do {
                            lower_squared = ((BoundaryPtr->lowerBound[i + 1]) * (BoundaryPtr->lowerBound[i + 1]) -
                                             2 * ds * maxForce / m + 2 * ds * g * (sin(EnvPtr->Angle_env[i]) +
                                                                                   crr * cos(EnvPtr->Angle_env[i]))) /
                                            (1 - 2 * ds * CdA / m);
                            maxForce = maxForce * 0.95;
                            counterBound++;
                            if (lower_squared < 0) {
                                break;
                            }
                        } while ((sqrt(lower_squared)) * maxForce > PAmax);

                        if (lower_squared < 0) {
                            break;
                        } else if ((sqrt(lower_squared)) <= BoundaryPtr->lowerBound[i]) {
                            break;
                        } else {
                            BoundaryPtr->lowerBound[i] = sqrt(lower_squared);
                        }
                    }
                }
            }
        }
    }

    free(changePoint);
}

void normalThermalBoundary(Boundary *BoundaryPtr, EnvFactor *EnvPtr) {
    for (uint16_t i = 0; i <= RES_THERMAL; i++) {
        BoundaryPtr->upperBound[i] = EnvPtr->T_required[i * HORIZON / RES_THERMAL] + 1;
        BoundaryPtr->lowerBound[i] = EnvPtr->T_required[i * HORIZON / RES_THERMAL] - 1;
    }
}

void customThermalBoundary(Boundary *BoundaryPtr, SolverInput *SolverInputPtr, DynParameter *ParaPtr, EnvFactor *EnvPtr,
                           Bridge *BridgePtr, real_T X0) {
    // Iteration Index
    uint16_t i, j, k;

    // Local Copy of Horizon Length
    uint16_t Nhrz = SolverInputPtr->GridSize.ResThermal;

    uint16_t *tmpMemory = (uint16_t *) malloc(Nhrz * sizeof(uint16_t));

    uint16_t counter = 0;
    // Detect the change point and number of blocks
    for (i = 0; i < Nhrz; i++) {
        if (EnvPtr->T_required[i * HORIZON / RES_THERMAL] != EnvPtr->T_required[(i + 1) * HORIZON / RES_THERMAL]) {
            tmpMemory[counter] = i + 1;
            counter++;
        }
    }

    // For instance, if the T_required changes at 3, 10. Then changPoint[] = [0, 3, 10];
    uint16_t numBlock = counter + 1;
    uint16_t *changePoint = (uint16_t *) malloc((numBlock + 1) * sizeof(uint16_t));
    changePoint[0] = 0;
    memcpy(changePoint + 1, tmpMemory, counter * sizeof(uint16_t));
    free(tmpMemory);

//    for (i = 0; i < numBlock; i++) {
//        printf("Change point: %d\n", changePoint[i]);
//    }

    // Local copy of constraints
    real_T Qmax = SolverInputPtr->Constraint.Qmax;
    real_T Qmin = SolverInputPtr->Constraint.Qmin;
    real_T PACmax = SolverInputPtr->Constraint.PACmax;
    real_T Tmax_inlet = SolverInputPtr->Constraint.Tmax_inlet;
    real_T Tmin_inlet = SolverInputPtr->Constraint.Tmin_inlet;

    // Local copy of parameters
    real_T Cth = ParaPtr->Cth;
    real_T Rth = ParaPtr->Rth;
    real_T Qsun = ParaPtr->Qsun;
    real_T Qpas = ParaPtr->Qpas;
    real_T Cp = ParaPtr->Cp;
    real_T rho = ParaPtr->rho;
    real_T mDot = ParaPtr->mDot;
    real_T CoP_pos = ParaPtr->CoP_pos;
    real_T CoP_neg = ParaPtr->CoP_neg;
    real_T Tamb = ParaPtr->Tamb;

    real_T beta0 = ParaPtr->beta0;

    real_T maxQ;
    real_T minQ;
    real_T upperTmp;
    real_T lowerTmp;

    // First, Copy the known required temperature to the upper and lower bounds
    for (i = 0; i <= Nhrz; i++) {
        BoundaryPtr->upperBound[i] = EnvPtr->T_required[i * HORIZON / RES_THERMAL] + 5;
        BoundaryPtr->lowerBound[i] = EnvPtr->T_required[i * HORIZON / RES_THERMAL] - 5;
    }

    // Let upper and lower bound both start from the known X0
    BoundaryPtr->upperBound[0] = X0;
    BoundaryPtr->lowerBound[0] = X0;

    // Calculate a sequence of tDelta
    real_T *tDelta = (real_T *) malloc(Nhrz * sizeof(real_T));
    for (i = 0; i < Nhrz; i++) {
        tDelta[i] = 0;
        for (j = 0; j < HORIZON / RES_THERMAL; j++) {
            tDelta[i] += BridgePtr->tDelta[i * HORIZON / RES_THERMAL + j];
        }
        printf("Duration_boundary: %f\n", tDelta[i]);
    }


    // If there is only one range block
    if (numBlock == 1) {
        // Draw Upper Bound
        for (i = 0; i < Nhrz; i++) {
            // Maximum possible Qinlet
            if (Qmax > Cp * rho * mDot * (Tmax_inlet - BoundaryPtr->upperBound[i])) {
                maxQ = Cp * rho * mDot * (Tmax_inlet - BoundaryPtr->upperBound[i]);
            } else {
                maxQ = Qmax;
            }

            upperTmp = BoundaryPtr->upperBound[i] +
                       (tDelta[i] / Cth) * (maxQ + Qsun + Qpas + (Tamb - BoundaryPtr->upperBound[i]) / Rth);

            counterBound++;

            // If it has exceeded the upper bound (required + 1), break
            if (upperTmp >= BoundaryPtr->upperBound[i + 1]) {
                break;
            }
                // If not, we give it a new (lower) value to it
            else {
                BoundaryPtr->upperBound[i + 1] = upperTmp;
            }
        }

        // Draw Lower Bound
        for (i = 0; i < Nhrz; i++) {
            // Minimum possible Qinlet
            if (Qmin < Cp * rho * mDot * (Tmin_inlet - BoundaryPtr->lowerBound[i])) {
                minQ = Cp * rho * mDot * (Tmin_inlet - BoundaryPtr->lowerBound[i]);
            } else {
                minQ = Qmin;
            }

            lowerTmp = BoundaryPtr->lowerBound[i] +
                       (tDelta[i] / Cth) * (minQ + Qsun + Qpas + (Tamb - BoundaryPtr->lowerBound[i]) / Rth);

            counterBound++;

            // If it has exceeded the lower bound (required - 1), break
            if (lowerTmp <= BoundaryPtr->lowerBound[i + 1]) {
                break;
            }
                // If not, we give it a new (higher) value to it
            else {
                BoundaryPtr->lowerBound[i + 1] = lowerTmp;
            }
        }
    }
        // If there are multiple range blocks
    else {
        for (k = 0; k < numBlock; k++) {
            // #1 Block
            if (k == 0) {
                // Draw Upper Bound
                for (i = 0; i < (changePoint[k + 1] - 1); i++) {
                    // Maximum possible Qinlet
                    if (Qmax > Cp * rho * mDot * (Tmax_inlet - BoundaryPtr->upperBound[i])) {
                        maxQ = Cp * rho * mDot * (Tmax_inlet - BoundaryPtr->upperBound[i]);
                    } else {
                        maxQ = Qmax;
                    }

                    upperTmp = BoundaryPtr->upperBound[i] +
                               (tDelta[i] / Cth) * (maxQ + Qsun + Qpas + (Tamb - BoundaryPtr->upperBound[i]) / Rth);

                    counterBound++;

                    // Same as when numBlock == 1
                    if (upperTmp >= BoundaryPtr->upperBound[i + 1]) {
                        break;
                    } else {
                        BoundaryPtr->upperBound[i + 1] = upperTmp;
                    }
                }

                // Draw Lower Bound
                for (i = 0; i < (changePoint[k + 1] - 1); i++) {
                    // Minimum possible Qinlet
                    if (Qmin < Cp * rho * mDot * (Tmin_inlet - BoundaryPtr->lowerBound[i])) {
                        minQ = Cp * rho * mDot * (Tmin_inlet - BoundaryPtr->lowerBound[i]);
                    } else {
                        minQ = Qmin;
                    }

                    lowerTmp = BoundaryPtr->lowerBound[i] +
                               (tDelta[i] / Cth) * (minQ + Qsun + Qpas + (Tamb - BoundaryPtr->lowerBound[i]) / Rth);

                    counterBound++;

                    // Same as when numBlock == 1
                    if (lowerTmp <= BoundaryPtr->lowerBound[i + 1]) {
                        break;
                    } else {
                        BoundaryPtr->lowerBound[i + 1] = lowerTmp;
                    }
                }
            }
                // From the #2 Block (2 conditions)
            else {
                // If T_required[Block_k] > T_required[Block_k-1], draw the Left-Top upper bound (forwards) and Right-Bottom lower bound (backwards)
                // (If required Temperature increases in the next range)
                if (EnvPtr->T_required[changePoint[k] * HORIZON / RES_THERMAL] >
                    EnvPtr->T_required[changePoint[k - 1] * HORIZON / RES_THERMAL]) {

                    // Upper Bound
                    for (i = changePoint[k] - 1; i < (changePoint[k + 1] - 1); i++) {
                        // Maximum possible Qinlet
                        if (Qmax > Cp * rho * mDot * (Tmax_inlet - BoundaryPtr->upperBound[i])) {
                            maxQ = Cp * rho * mDot * (Tmax_inlet - BoundaryPtr->upperBound[i]);
                        } else {
                            maxQ = Qmax;
                        }

                        upperTmp = BoundaryPtr->upperBound[i] +
                                   (tDelta[i] / Cth) * (maxQ + Qsun + Qpas + (Tamb - BoundaryPtr->upperBound[i]) / Rth);

                        counterBound++;

                        // Same
                        if (upperTmp >= BoundaryPtr->upperBound[i + 1]) {
                            break;
                        } else {
                            BoundaryPtr->upperBound[i + 1] = upperTmp;
                        }
                    }

                    // Lower Bound
                    for (i = changePoint[k] - 1; i >= changePoint[k - 1]; i--) {
                        maxQ = Qmax;

                        do {
                            lowerTmp = ((Cth / tDelta[i]) * BoundaryPtr->lowerBound[i + 1] - maxQ - Qsun - Qpas -
                                        Tamb / Rth) / (Cth / tDelta[i] - 1 / Rth);
                            maxQ = maxQ * 0.95;
                            counterBound++;
                        } while (lowerTmp + maxQ / (Cp * rho * mDot) > Tmax_inlet);

                        if (lowerTmp <= BoundaryPtr->lowerBound[i]) {
                            break;
                        } else {
                            BoundaryPtr->lowerBound[i] = lowerTmp;
                        }
                    }
                }
                    // If T_required[Block_k] > T_required[Block_k-1], draw the Right-Top upper bound (backwards) and Left-Bottom lower bound (forwards)
                else {

                    // Upper Bound
                    for (i = changePoint[k] - 1; i >= changePoint[k - 1]; i--) {
                        // Predictive minimum Qinlet
                        minQ = Qmin;

                        // Calculate the upper bound backwards
                        do {
                            upperTmp = ((Cth / tDelta[i]) * BoundaryPtr->upperBound[i + 1] - minQ - Qsun - Qpas -
                                        Tamb / Rth) / (Cth / tDelta[i] - 1 / Rth);
                            minQ = minQ * 0.95;
                            counterBound++;
                        } while (upperTmp + minQ / (Cp * rho * mDot) < Tmin_inlet);

                        if (upperTmp >= BoundaryPtr->upperBound[i]) {
                            break;
                        } else {
                            BoundaryPtr->upperBound[i] = upperTmp;
                        }
                    }

                    // Lower Bound
                    for (i = changePoint[k] - 1; i < (changePoint[k + 1] - 1); i++) {
                        // Minimum possible Qinlet
                        if (Qmin < Cp * rho * mDot * (Tmin_inlet - BoundaryPtr->lowerBound[i])) {
                            minQ = Cp * rho * mDot * (Tmin_inlet - BoundaryPtr->lowerBound[i]);
                        } else {
                            minQ = Qmin;
                        }

                        lowerTmp = BoundaryPtr->lowerBound[i] +
                                   (tDelta[i] / Cth) * (minQ + Qsun + Qpas + (Tamb - BoundaryPtr->lowerBound[i]) / Rth);

                        counterBound++;

                        // Same as when numBlock == 1
                        if (lowerTmp <= BoundaryPtr->lowerBound[i + 1]) {
                            break;
                        } else {
                            BoundaryPtr->lowerBound[i + 1] = lowerTmp;
                        }
                    }
                }
            }
        }
    }

    free(tDelta);
    free(changePoint);
}

#endif