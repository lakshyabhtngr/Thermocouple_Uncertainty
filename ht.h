/**********************************************************************\
File:         ht.h

Purpose:      Tables for HT Loss Calculations
\***********************************************************************/

  /* Temperatures */
  #define T_Max                           10
  #define MN_Max                          6
  
  /* Physical Constants */
  #define Rankine                         459.67
  #define StefanBoltzman                  4.76E-13
  #define Emissivity                      0.85
  #define U_Emissivity                    0.05
  #define AirDensityConstant              2.7

  /* Inlet Pressures (PSIG) */
  #define CIT_P  14.7
  #define FEX_P  35.0
  #define CDT_P  300.0
  #define TIT_P  300.0

  /* Mach Calculation Constants */
  #define g_C0               1.3930336 
  #define g_C1               6.81374E-5
  #define g_C2              -1.11831E-7
  #define g_C3               3.16776E-11
  #define GravityConstat     32.174       /* FT/SEC^2 */
  #define GasConstant        1545.33      /* LBM/LBMOL */
  #define AirMolecularWeight 28.9664      /* LBF/(LBMOL*DEGR) */

  /* Dynamic Air Viscosity (u) */
  /* u = C0 + T*(C1 + T*C2) */
  #define DAV_C0    4.141794E-06 /* Dynamic Air Viscosity C0 */
  #define DAV_C1    1.709688E-08 /* Dynamic Air Viscosity C1 */
  #define DAV_C2   -2.127130E-12 /* Dynamic Air Viscosity C2 */
  #define TDAV_MIN   559.67      /* DEGR */
  #define TDAV_MAX  2459.67      /* DEGR */
  #define DAV_MIN   1.274722E-05 /* LB/(FT*SEC) */
  #define DAV_MAX   3.351389E-05 /* LB/(FT*SEC) */
  #define U_DAV     0.50         /* % */

  /* Air Thermo Conductivity (kf) */
  /* kf = C0 + T*(C1 + T*C2) */
  #define ATC_C0    8.72732655E-07   /* Air Thermo Conductivity C0 */
  #define ATC_C1    6.84123412E-09   /* Air Thermo Conductivity C1 */
  #define ATC_C2   -7.18448218E-13   /* Air Thermo Conductivity C2 */
  #define TATC_MIN   559.67      /* DEGR */
  #define TATC_MAX  2459.67      /* DEGR */
  #define ATC_MIN   4.3639E-06   /* BTU/(SEC*FT*DEGR) */
  #define ATC_MAX   1.3725E-05   /* BTU/(SEC*FT*DEGR) */
  #define U_ATC     0.50         /* % */

  /* Type K Thermo Conductivity (ks) */
  /* ks(K) = C0 + T*(C1 + T*C2) */
  #define KTC_C0    2.122727E-02 /* K Thermo Conductivity C0 */
  #define KTC_C1   -1.472331E-05 /* K Thermo Conductivity C1 */
  #define KTC_C2    4.983878E-09 /* K Thermo Conductivity C2 */
  #define TKTC_MIN   360.27      /* DEGR */
  #define TKTC_MAX  2160.27      /* DEGR */
  #define KTC_MIN   1.656978E-02 /* BTU/(SEC*FT*DEGR) */
  #define KTC_MAX   1.267954E-02 /* BTU/(SEC*FT*DEGR) */
  #define U_KTC     0.50         /* % */

  /* Type J Thermo Conductivity (ks) */
  /* ks(J) = C0 + T*(C1 + T*C2) */
  #define JTC_C0    1.007941E-02 /* J Thermo Conductivity C0 */
  #define JTC_C1   -4.839988E-06 /* J Thermo Conductivity C1 */
  #define JTC_C2    9.680179E-10 /* J Thermo Conductivity C2 */
  #define TJTC_MIN   360.27      /* DEGR */
  #define TJTC_MAX  2160.27      /* DEGR */
  #define JTC_MIN   1.656978E-02 /* BTU/(SEC*FT*DEGR) */
  #define JTC_MAX   4.141243E-03 /* BTU/(SEC*FT*DEGR) */
  #define U_JTC     0.50         /* % */

  /* Type E Thermo Conductivity (ks) */
  /* ks(E) = C0 + T*(C1 + T*C2) */
  #define ETC_C0    1.105492E-02  /* E Thermo Conductivity C0 */
  #define ETC_C1    5.928694E-06  /* E Thermo Conductivity C1 */
  #define ETC_C2    2.133430E-09 /* E Thermo Conductivity C2 */
  #define TETC_MIN   360.27      /* DEGR */
  #define TETC_MAX  2160.27      /* DEGR */
  #define ETC_MIN   1.346776E-02 /* BTU/(SEC*FT*DEGR) */
  #define ETC_MAX   2.864256E-02 /* BTU/(SEC*FT*DEGR) */
  #define U_ETC     0.50         /* % */

  /* Type N Thermo Conductivity (ks) */
  /* ks(N) = C0 + T*(C1 + T*C2) */
  #define NTC_C0    2.430266E-03 /* N Thermo Conductivity C0 */
  #define NTC_C1    8.074214E-07 /* N Thermo Conductivity C1 */
  #define NTC_C2   -2.373496E-10 /* N Thermo Conductivity C2 */
  #define TNTC_MIN   360.27      /* DEGR */
  #define TNTC_MAX  2160.27      /* DEGR */
  #define NTC_MIN   2.690349E-03 /* BTU/(SEC*FT*DEGR) */
  #define NTC_MAX   3.066859E-03 /* BTU/(SEC*FT*DEGR) */
  #define U_NTC     0.50         /* % */

  /* Assumed Wall Conditions */
  #define WALL_SLOPE                     0.08      /* DEGR/DEGR */
  #define WALL_OFFSET                   -21.5736   /* DEGR */
  #define ParallelMachRecovery            0.96
  #define ParallelMachRecoverySpread      0.03
  #define NormalMachRecovery              0.68
  #define NormalMachRecoverySpread        0.07

