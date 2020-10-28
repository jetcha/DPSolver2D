#include "./inc/SolverStruct.h"
#include "./inc/SolverAlgorithm.h"

int main() {
    /*----------------------*/
    /*--- Initialization ---*/
    /*----------------------*/

    // Input Variable declaration (will be used in the BasicAlgorithm)
    SolverInput SolverInputPtr;                      // solver inputs: Constraints, Grid Sizes, Road Profile, etc.
    DynParameter ModelParaPtr;                       // solver inputs: Model parameters
    EnvFactor EnvFactorPtr;                          // solver inputs: Environmental Factors (Legal speed limits, Angle of slopes)
    SolverOutput SolverOutputPtr;                    // solver outputs: Minimum Total Cost, Optimal Speed Trajectory, Optimal Control Policy
    real_T V0;                                       // initial speed
    real_T T0;                                       // initial temperature

    // Input Settings
    SolverInputPtr.GridSize.Nv = NV;
    SolverInputPtr.GridSize.Nf = NF;
    SolverInputPtr.GridSize.Nt = NT;
    SolverInputPtr.GridSize.Nq = NQ;
    SolverInputPtr.GridSize.Nhrz = HORIZON;

    SolverInputPtr.Constraint.Vmax = 130 / 3.6;        // Physical Speed limits
    SolverInputPtr.Constraint.Vmin = 0.0;
    SolverInputPtr.Constraint.Fmax = 4e3;
    SolverInputPtr.Constraint.Fmin = -4e3;
    SolverInputPtr.Constraint.PAmax = 6.5e4;
    SolverInputPtr.Constraint.PDmax = -6.5e4;
    SolverInputPtr.Constraint.Tmax = 26;
    SolverInputPtr.Constraint.Tmin = 24;
    SolverInputPtr.Constraint.Tmax_inlet = 40;
    SolverInputPtr.Constraint.Tmin_inlet = 5;
    SolverInputPtr.Constraint.Qmax = 2000;
    SolverInputPtr.Constraint.Qmin = -2000;
    SolverInputPtr.Constraint.PACmax = 2000 / 1.8;

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
    ModelParaPtr.speedPenalty = 13e4;
    ModelParaPtr.thermalPenalty = 4e4;

    // Environmental Information
    uint16_t Nhrz = HORIZON;
    uint16_t i;

#ifdef SCENE1
    // Initial Speed
    V0 = 0 / 3.6;
    // Initial Temperature
    T0 = 26;
    // Block Length
    real_T blockLength = 750;
    // Block Steps
    real_T blockStep = blockLength / ModelParaPtr.ds;

    real_T Vmax_GPS_1 = 50 / 3.6;
    real_T Vmin_GPS_1 = 0 / 3.6;
    real_T T_required_1 = 25;
    uint16_t endBlock_1 = (uint16_t) 1 * (blockStep - 0);
    //uint16_t endBlock_1 = 18;

    real_T Vmax_GPS_2 = 80 / 3.6;
    real_T Vmin_GPS_2 = 30 / 3.6;
    real_T T_required_2 = 25;
    uint16_t endBlock_2 = (uint16_t) 2 * blockStep;
    //uint16_t endBlock_2 = 37;

    real_T Vmax_GPS_3 = 100 / 3.6;
    real_T Vmin_GPS_3 = 50 / 3.6;
    real_T T_required_3 = 25;
    uint16_t endBlock_3 = (uint16_t) 3 * (blockStep - 0);
    //uint16_t endBlock_3 = 56;

    real_T Vmax_GPS_4 = 130 / 3.6;
    real_T Vmin_GPS_4 = 50 / 3.6;
    real_T T_required_4 = 25;
    uint16_t endBlock_4 = (uint16_t) 4 * blockStep;
    //uint16_t endBlock_4 = 75;

    real_T Vmax_GPS_5 = 50 / 3.6;
    real_T Vmin_GPS_5 = 0 / 3.6;
    real_T T_required_5 = 25;
    uint16_t endBlock_5 = (uint16_t) 5 * (blockStep - 0);
    //uint16_t endBlock_5 = 93;

    real_T Vmax_GPS_6 = 130 / 3.6;
    real_T Vmin_GPS_6 = 50 / 3.6;
    real_T T_required_6 = 25;
    uint16_t endBlock_6 = (uint16_t) 6 * blockStep;
    //uint16_t endBlock_6 = 112;

    real_T Vmax_GPS_7 = 80 / 3.6;
    real_T Vmin_GPS_7 = 30 / 3.6;
    real_T T_required_7 = 25;
    uint16_t endBlock_7 = (uint16_t) 7 * (blockStep - 0);
    //uint16_t endBlock_7 = 131;

    real_T Vmax_GPS_8 = 100 / 3.6;
    real_T Vmin_GPS_8 = 50 / 3.6;
    real_T T_required_8 = 25;
    uint16_t endBlock_8 = Nhrz;

    EnvFactorPtr.endBlock[0] = endBlock_1;
    EnvFactorPtr.endBlock[1] = endBlock_2;
    EnvFactorPtr.endBlock[2] = endBlock_3;
    EnvFactorPtr.endBlock[3] = endBlock_4;
    EnvFactorPtr.endBlock[4] = endBlock_5;
    EnvFactorPtr.endBlock[5] = endBlock_6;
    EnvFactorPtr.endBlock[6] = endBlock_7;
    EnvFactorPtr.endBlock[7] = endBlock_8;


    for (i = 0; i <= Nhrz; i++) {
        if (i < endBlock_1) {
            EnvFactorPtr.Vmax_env[i] = Vmax_GPS_1;                  // Legal Vmax
            EnvFactorPtr.Vmin_env[i] = Vmin_GPS_1;                  // Legal Vmin
            EnvFactorPtr.Angle_env[i] = 0.0;                        // Road slops
            EnvFactorPtr.T_required[i] = T_required_1;              // Required Temp
        } else if (i < endBlock_2) {
            EnvFactorPtr.Vmax_env[i] = Vmax_GPS_2;                  // Legal Vmax
            EnvFactorPtr.Vmin_env[i] = Vmin_GPS_2;                  // Legal Vmin
            EnvFactorPtr.Angle_env[i] = 0.0;                        // Road slops
            EnvFactorPtr.T_required[i] = T_required_2;              // Required Temp
        } else if (i < endBlock_3) {
            EnvFactorPtr.Vmax_env[i] = Vmax_GPS_3;                  // Legal Vmax
            EnvFactorPtr.Vmin_env[i] = Vmin_GPS_3;                  // Legal Vmin
            EnvFactorPtr.Angle_env[i] = 0.0;                        // Road slops
            EnvFactorPtr.T_required[i] = T_required_3;              // Required Temp
        } else if (i < endBlock_4) {
            EnvFactorPtr.Vmax_env[i] = Vmax_GPS_4;                  // Legal Vmax
            EnvFactorPtr.Vmin_env[i] = Vmin_GPS_4;                  // Legal Vmin
            EnvFactorPtr.Angle_env[i] = 0.0;                        // Road slops
            EnvFactorPtr.T_required[i] = T_required_4;              // Required Temp
        } else if (i < endBlock_5) {
            EnvFactorPtr.Vmax_env[i] = Vmax_GPS_5;                  // Legal Vmax
            EnvFactorPtr.Vmin_env[i] = Vmin_GPS_5;                  // Legal Vmin
            EnvFactorPtr.Angle_env[i] = 0.0;                        // Road slops
            EnvFactorPtr.T_required[i] = T_required_5;              // Required Temp
        } else if(i < endBlock_6) {
            EnvFactorPtr.Vmax_env[i] = Vmax_GPS_6;                  // Legal Vmax
            EnvFactorPtr.Vmin_env[i] = Vmin_GPS_6;                  // Legal Vmin
            EnvFactorPtr.Angle_env[i] = 0.0;                        // Road slops
            EnvFactorPtr.T_required[i] = T_required_6;              // Required Temp
        } else if(i < endBlock_7) {
            EnvFactorPtr.Vmax_env[i] = Vmax_GPS_7;                  // Legal Vmax
            EnvFactorPtr.Vmin_env[i] = Vmin_GPS_7;                  // Legal Vmin
            EnvFactorPtr.Angle_env[i] = 0.0;                        // Road slops
            EnvFactorPtr.T_required[i] = T_required_7;              // Required Temp
        } else {
            EnvFactorPtr.Vmax_env[i] = Vmax_GPS_8;                  // Legal Vmax
            EnvFactorPtr.Vmin_env[i] = Vmin_GPS_8;                  // Legal Vmin
            EnvFactorPtr.Angle_env[i] = 0.0;                        // Road slops
            EnvFactorPtr.T_required[i] = T_required_8;              // Required Temp
        }
    }
#elif defined(SCENE2)
    // Initial Speed
    V0 = 0/3.6;
    // Initial Temperature
    T0 = 26;
    // Block Lengths
    real_T blockLength1 = 300;
    real_T blockLength2 = 2700;
    real_T blockLength3 = 3300;

    real_T Vmax_GPS_1 = 100 / 3.6;
    real_T Vmin_GPS_1 = 0 / 3.6;
    real_T T_required_1 = 25;
    uint16_t endBlock_1 = (uint16_t) blockLength1/ ModelParaPtr.ds;
    //uint16_t endBlock_1 = 7;

    real_T Vmax_GPS_2 = 100 / 3.6;
    real_T Vmin_GPS_2 = 50 / 3.6;
    real_T T_required_2 = 25;
    uint16_t endBlock_2 = (uint16_t) blockLength2/ ModelParaPtr.ds;
    //uint16_t endBlock_2 = 68;

    real_T Vmax_GPS_3 = 50 / 3.6;
    real_T Vmin_GPS_3 = 0 / 3.6;
    real_T T_required_3 = 25;
    uint16_t endBlock_3 = (uint16_t) blockLength3/ ModelParaPtr.ds;
    //uint16_t endBlock_3 = 83;

    real_T Vmax_GPS_4 = 100 / 3.6;
    real_T Vmin_GPS_4 = 50 / 3.6;
    real_T T_required_4 = 25;
    uint16_t endBlock_4 = Nhrz;

    EnvFactorPtr.endBlock[0] = endBlock_1;
    EnvFactorPtr.endBlock[1] = endBlock_2;
    EnvFactorPtr.endBlock[2] = endBlock_3;
    EnvFactorPtr.endBlock[3] = endBlock_4;

    for (i = 0; i <= Nhrz; i++)
    {
        if (i < endBlock_1) {
            EnvFactorPtr.Vmax_env[i] = Vmax_GPS_1;                  // Legal Vmax
            EnvFactorPtr.Vmin_env[i] = Vmin_GPS_1;                  // Legal Vmin
            EnvFactorPtr.Angle_env[i] = 0.0;                        // Road slops
            EnvFactorPtr.T_required[i] = T_required_1;              // Required Temp
        } else if (i < endBlock_2) {
            EnvFactorPtr.Vmax_env[i] = Vmax_GPS_2;                  // Legal Vmax
            EnvFactorPtr.Vmin_env[i] = Vmin_GPS_2;                  // Legal Vmin
            EnvFactorPtr.Angle_env[i] = 0.0;                        // Road slops
            EnvFactorPtr.T_required[i] = T_required_2;              // Required Temp
        } else if (i < endBlock_3) {
            EnvFactorPtr.Vmax_env[i] = Vmax_GPS_3;                  // Legal Vmax
            EnvFactorPtr.Vmin_env[i] = Vmin_GPS_3;                  // Legal Vmin
            EnvFactorPtr.Angle_env[i] = 0.0;                        // Road slops
            EnvFactorPtr.T_required[i] = T_required_3;              // Required Temp
        } else {
            EnvFactorPtr.Vmax_env[i] = Vmax_GPS_4;                  // Legal Vmax
            EnvFactorPtr.Vmin_env[i] = Vmin_GPS_4;                  // Legal Vmin
            EnvFactorPtr.Angle_env[i] = 0.0;                        // Road slops
            EnvFactorPtr.T_required[i] = T_required_4;              // Required Temp
        }
    }
#endif

    /*-------------------------*/
    /*--- Run the Algorithm ---*/
    /*-------------------------*/
    MagicBox(&SolverInputPtr, &ModelParaPtr, &EnvFactorPtr, &SolverOutputPtr, V0, T0);

    return 0;
}