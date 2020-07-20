#include "./inc/SolverStruct.h"
#include "./inc/BasicAlgorithm.h"

int main() {
    /*----------------------*/
    /*--- Initialization ---*/
    /*----------------------*/

    // Input Variable declaration (will be used in the BasicAlgorithm)
    SolverInput SolverInputPtr;                        // solver inputs: Constraints, Grid Sizes, Road Profile, etc.
    DynParameter ModelParaPtr;                        // solver inputs: Model parameters
    EnvFactor EnvFactorPtr;                            // solver inputs: Environmental Factors (Legal speed limits, Angle of slopes)
    SolverOutput SolverOutputPtr;                    // solver outputs: Minimum Total Cost, Optimal Speed Trajectory, Optimal Control Policy
    real_T V0;                                        // initial speed
    real_T T0;                                        // initial temperature
    real_T Vfmin;                                    // final speed constraints
    real_T Vfmax;
    real_T Tfmin;                                    // final thermal constraints
    real_T Tfmax;

    // Final states constraint
    Vfmin = 0 / 3.6;
    Vfmax = 400 / 3.6;
    Tfmin = 18;
    Tfmax = 30;

    // Input Settings
    SolverInputPtr.GridSize.Nv = NV;
    SolverInputPtr.GridSize.Nf = NF;
    SolverInputPtr.GridSize.Nt = NT;
    SolverInputPtr.GridSize.Nq = NQ;
    SolverInputPtr.GridSize.Nhrz = HORIZON;
    SolverInputPtr.GridSize.ResThermal = RES_THERMAL;

    SolverInputPtr.Constraint.Vmax = 200 / 3.6;        // Physical Speed limits
    SolverInputPtr.Constraint.Vmin = 0.0;
    SolverInputPtr.Constraint.Fmax = 3e3;
    SolverInputPtr.Constraint.Fmin = -3e3;
    SolverInputPtr.Constraint.PAmax = 6e4;
    SolverInputPtr.Constraint.PDmax = -6e4;
    SolverInputPtr.Constraint.Tmax = 30;
    SolverInputPtr.Constraint.Tmin = 10;
    SolverInputPtr.Constraint.Tmax_inlet = 30;
    SolverInputPtr.Constraint.Tmin_inlet = 10;
    SolverInputPtr.Constraint.Qmax = 2000;
    SolverInputPtr.Constraint.Qmin = -2000;
    SolverInputPtr.Constraint.PACmax = 2000 / 2.1379;

    SolverInputPtr.SolverLimit.infValue = FLT_MAX;

    // Speed Parameters
    ModelParaPtr.m = 2000;
    ModelParaPtr.g = 9.81;
    ModelParaPtr.CdA = 0.6;
    ModelParaPtr.crr = 0.006;

    ModelParaPtr.alpha0 = 785.0;
    ModelParaPtr.alpha1 = -43e-4;
    ModelParaPtr.alpha2 = 5.6e-7;
    ModelParaPtr.beta0 = 41e-7;
    ModelParaPtr.eta_dc = 0.99;
    ModelParaPtr.eta_trans = 0.98;

    // Thermal Parameters
    ModelParaPtr.Cth = 1.1347e5;
    ModelParaPtr.Rth = 0.015;
    ModelParaPtr.Qsun = 150 * 1.7725;
    ModelParaPtr.Qpas = 4 * 104.2975;
    ModelParaPtr.Cp = 1.0035e3;
    ModelParaPtr.rho = 1.1839;
    ModelParaPtr.mDot = 0.0842;
    ModelParaPtr.CoP_pos = 2.1379;
    ModelParaPtr.CoP_neg = -2.1379;
    ModelParaPtr.Tamb = 30;

    // Tuning Parameter
    ModelParaPtr.ds = 30;
    ModelParaPtr.speedPenalty = 11e4;
    ModelParaPtr.thermalPenalty = 10e4;

    // Environmental Information
    uint8_t numFactors = 3;
    uint16_t Nhrz = HORIZON;
    uint16_t i;

#ifdef SCENE1
    // Initial Speed
    V0 = 0 / 3.6;
    // Initial Temperature
    T0 = 30;

    real_T Vmax_GPS_1 = 50 / 3.6;
    real_T Vmin_GPS_1 = 0 / 3.6;
    real_T T_required_1 = 25;
    uint16_t endBlock_1 = 100;

    real_T Vmax_GPS_2 = 80 / 3.6;
    real_T Vmin_GPS_2 = 30 / 3.6;
    real_T T_required_2 = 25;
    uint16_t endBlock_2 = 250;

    real_T Vmax_GPS_3 = 130 / 3.6;
    real_T Vmin_GPS_3 = 60 / 3.6;
    real_T T_required_3 = 25;
    uint16_t endBlock_3 = 750;

    real_T Vmax_GPS_4 = 80 / 3.6;
    real_T Vmin_GPS_4 = 30 / 3.6;
    real_T T_required_4 = 25;
    uint16_t endBlock_4 = 850;

    real_T Vmax_GPS_5 = 50 / 3.6;
    real_T Vmin_GPS_5 = 0 / 3.6;
    real_T T_required_5 = 25;
    uint16_t endBlock_5 = Nhrz;

    EnvFactorPtr.endBlock[0] = endBlock_1;
    EnvFactorPtr.endBlock[1] = endBlock_2;
    EnvFactorPtr.endBlock[2] = endBlock_3;
    EnvFactorPtr.endBlock[3] = endBlock_4;
    EnvFactorPtr.endBlock[4] = endBlock_5;

    for (i = 0; i <= Nhrz; i++) {
        if (i < endBlock_1) {
            EnvFactorPtr.Vmax_env[i] = Vmax_GPS_1;                    // Legal Vmax
            EnvFactorPtr.Vmin_env[i] = Vmin_GPS_1;                    // Legal Vmin
            EnvFactorPtr.Angle_env[i] = 0.0;                        // Road slops
            EnvFactorPtr.T_required[i] = T_required_1;              // Required Temp
        } else if (i < endBlock_2) {
            EnvFactorPtr.Vmax_env[i] = Vmax_GPS_2;                    // Legal Vmax
            EnvFactorPtr.Vmin_env[i] = Vmin_GPS_2;                    // Legal Vmin
            EnvFactorPtr.Angle_env[i] = 0.0;                        // Road slops
            EnvFactorPtr.T_required[i] = T_required_2;              // Required Temp
        } else if (i < endBlock_3) {
            EnvFactorPtr.Vmax_env[i] = Vmax_GPS_3;                    // Legal Vmax
            EnvFactorPtr.Vmin_env[i] = Vmin_GPS_3;                    // Legal Vmin
            EnvFactorPtr.Angle_env[i] = 0.0;                        // Road slops
            EnvFactorPtr.T_required[i] = T_required_3;              // Required Temp
        } else if (i < endBlock_4) {
            EnvFactorPtr.Vmax_env[i] = Vmax_GPS_4;                    // Legal Vmax
            EnvFactorPtr.Vmin_env[i] = Vmin_GPS_4;                    // Legal Vmin
            EnvFactorPtr.Angle_env[i] = 0.0;                        // Road slops
            EnvFactorPtr.T_required[i] = T_required_4;              // Required Temp
        } else {
            EnvFactorPtr.Vmax_env[i] = Vmax_GPS_5;                    // Legal Vmax
            EnvFactorPtr.Vmin_env[i] = Vmin_GPS_5;                    // Legal Vmin
            EnvFactorPtr.Angle_env[i] = 0.0;                        // Road slops
            EnvFactorPtr.T_required[i] = T_required_5;              // Required Temp
        }
    }
#elif defined(SCENE2)
    // Initial Speed
    X0 = 0/3.6;
    T0 = 28;

    real_T Vmax_GPS_1 = 130 / 3.6;
    real_T Vmin_GPS_1 = 0 / 3.6;

    for (i = 0; i <= Nhrz; i++)
    {
        EnvFactorPtr.Vmax_env[i] = Vmax_GPS_1;					// Legal Vmax
        EnvFactorPtr.Vmin_env[i] = Vmin_GPS_1;					// Legal Vmin
        EnvFactorPtr.Angle_env[i] = 0.0;						// Road slops

        EnvFactorPtr.T_required[i] = 24;
    }
#endif

    /*-------------------------*/
    /*--- Run the Algorithm ---*/
    /*-------------------------*/
    MagicBox(&SolverInputPtr, &ModelParaPtr, &EnvFactorPtr, &SolverOutputPtr, V0, T0, Vfmin, Vfmax, Tfmin, Tfmax);



//    for(i = 0; i<=RES_THERMAL;i++)
//    {
//        printf("Upper Thermal bound: %f, Lower Thermal bound: %f\n", SolverOutputPtr.upperTempBound[i], SolverOutputPtr.lowerTempBound[i]);
//    }

    // Fix
    return 0;
}