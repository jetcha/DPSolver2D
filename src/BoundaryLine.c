#include "../inc/BoundaryLine.h"


#if defined(CUSTOMBOUND)
void initSpeedBoundary(Boundary *BoundaryPtr) {
    BoundaryPtr->lowerBound = (real_T *) malloc((HORIZON + 1) * sizeof(real_T));
    BoundaryPtr->upperBound = (real_T *) malloc((HORIZON + 1) * sizeof(real_T));
}

void copySpeedBoundary(Boundary *BoundaryPtr, SolverOutput *OutputPtr) {
    memcpy(OutputPtr->lowerSpeedBound, BoundaryPtr->lowerBound, (HORIZON + 1) * sizeof(real_T));
    memcpy(OutputPtr->upperSpeedBound, BoundaryPtr->upperBound, (HORIZON + 1) * sizeof(real_T));
}

void freeBoundary(Boundary *BoundaryPtr) {
    free(BoundaryPtr->lowerBound);
    free(BoundaryPtr->upperBound);
}

void customSpeedBoundary(Boundary *BoundaryPtr, SolverInput *SolverInputPtr, DynParameter *ParaPtr, EnvFactor *EnvPtr,
                         real_T V0) {
    // Iteration Index
    uint16_t i, k;

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
    BoundaryPtr->upperBound[0] = V0;
    BoundaryPtr->lowerBound[0] = V0;

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


#endif