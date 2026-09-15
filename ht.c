/**********************************************************************\
File:         ht.c

Purpose:      Compute Heat Transfer Losses for Temperature Probes
\***********************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ht.h"

#define L_Max 6
static double L_Ratio[L_Max] = {10, 15, 20, 25, 30, 40};
static double T_CIT[T_Max] = { -80, -40,    0,   20,   40,   60,   80,  100,  150,  200};
static double T_FEX[T_Max] = {  32,  50,   75,  100,  150,  200,  250,  300,  350,  400};
static double T_CDT[T_Max] = { 200, 300,  400 , 500,  600,  800, 1000, 1200, 1350, 1500};
static double T_TIT[T_Max] = { 500, 750, 1000, 1250, 1500, 1750, 2000, 2200, 2400, 2500};
static double T_FIL[T_Max];
static double T_FIL[T_Max];
static double FIL_P;
static double MN_Table[MN_Max] = {0.2,0.3,0.4,0.5,0.6,0.7};

/***********************************************************************
 * Function lstsqr.c 
* 
*  Description:
*
*      LEAST SQUARES SUBROUTINE
*      REF: NUMERICAL METHODS BY ROBERT HORNBECK, QUANTUM PUBLISHERS
*           1975, PP 121-125 FOR THE GUASS-JORDAN ELIMINATION METHOD
*           FOR SOLVING THE SIMULTANOUS EQUATIONS
*           DESIGN AND ANALYSIS OF ENGINEERING EXEPERIMENTS, LIPSON AND
*           SHETH, MCGRAW-HILL,     PP   400   FOR THE  CORELATION
*           COEFFICIENT.
*
*  Parameters:
*
*     X(NPTS) - INDEPENDENT VARIABLE ARRAY OF SIZE NPTS
*     Y(NPTS) - DEPENDENT VARIABLE ARRAY OF SIZE NPTS
*     NPTS    - NUMBER OF PTS IN THE VARIABLE ARRAYS
*     ORDER   - ORDER OF THE POLYNOMIAL TO BE FITTED (1-6)
*     C( , )  - 'STIFFNESS' ARRAY OF THE SIMULTANEOUS EQUATIONS
*     R( )    - 'FORCE' ARRAY OF THE SIMULTAINOUS EQUATIONS
*                (THE SUMS OF (X**?)*Y)
* 
*************************************************************************/
#define  MAXCOEF   7

static void lstsqr(double *x, double *y, int npts, int order, double *ans)
{
  double  c[MAXCOEF+1][MAXCOEF+1], r[MAXCOEF+1];
  double  sumy, sumyy, yic, diff;
  double  temp, sy, syx;
  int i, j, k, l;

  sumy = sumyy = yic = diff = 0.0;

  /* initialize arrays to zero */

  for (j = 0; j < MAXCOEF+1; j++)
  {
    r[j] = 0.0;
    for (i = 0; i < MAXCOEF; i++)
      c[j][i] = 0.0;
  }

  /* generate the matricies c and r */

  for (k = 0; k < npts; k++)
  {
    for (j = 0; j < order +1; j++)
    {
      if (j == 0)
        r[j] = r[j] + y[k];
      else
        r[j] = r[j] + y[k] * pow(x[k],j);
      for (i = 0; i < order+1; i++)
      {
        if (i == 0 && j ==0)
          c[j][i] = c[j][i] + 1;
        else
          c[j][i] = c[j][i] + pow(x[k],i+j);
      }
    }
  }
  /*
    USE GAUSS-JORDAN ELIMINATION TO SOLVE THE EQUATIONS
    THE SOLUTION IS STORED IN THE R ARRAY (THE INITIAL VALUES
    OF THE TWO ARRAYS ARE DESTROYED)
  */

  for (k = 0; k < order + 1; k++)
  {
    temp = 1.0 / c[k][k];
    for (j = 0; j < order + 1; j++)
      c[k][j] = c[k][j] * temp;
    r[k] = r[k] * temp;
    for (j = 0; j < order + 1; j++)
    {
      if (j != k)
      {
        temp = c[j][k];
        for (l = k; l < order +1; l++)
          c[j][l] = c[j][l] - c[k][l] * temp;
        r[j] = r[j] - r[k] * temp;
      }
    }
  }
  /*
    CALCULATE THE CORRELATION COEFFICIENT SQUARED AND STORE IN
    THE LAST ELEMENT IN R. REFERENCE LIPSON AND SHETH
    NOTE THAT IF THE SAMPLE SIZE EQUALS OR IS LESS THAN THE
    NUMBER OF COEFFICIENTS IN THE EQUATION THE CORRELLATION
    COEFFICIENT IS EQUAL TO ONE. NEED THIS BECAUSE THE LINE
    'SYX=DIFF/(NPTS-ORDER-1)'
    WOULD GIVE A DIVIDE BY ZERO ERROR.
  */

  if (npts - order <= 0)
    r[order + 1] = 1.0;
  else
  {
    for  (i = 0; i < npts; i++)
    {
      sumy = sumy + y[i];
      sumyy = sumyy + pow(y[i], 2);
      yic = 0.0;
      for (j = 0; j < order + 1; j++)
        if (j == 0)
          yic = yic + r[j];
        else
          yic = yic + r[j] * pow(x[i], j);
      diff = diff + pow(y[i]-yic, 2);
     }
     sy = (npts * sumyy - pow(sumy,2)) / npts / (npts-1);
     syx = diff / (npts - order);
     r[order + 1] = fabs(1.0 - syx/sy);
  }   
  /*
    EQUATE THE SOLUTIONS FOUND THAT ARE STORED IN THE DBLE
    PRECISION R ARRAY TO THE SINGLE PRECISION ANS ARRAY AND
    RETURN
  */
  for  (i = 0; i < order+2; i++)
    ans[i] = r[i];
}  

/***********************************************************************/
static void GetTempAndPres()
/*
Purpose:      Get Temperature Range and Pressure From File
              File name data.cfg and located where program is run
\***********************************************************************/
{
     int   i;
    char   DataName[50];
    char   LineBuf[128];
    FILE  *DataFile;

  /* Open the data File */
  sprintf(DataName, "./data.cfg");
  if ((DataFile = fopen(DataName,"rt")) == NULL)
  {
    printf("File Not Found\n");
    exit(0);
  }
  /* Get the Temperature Must be 10 */
  for (i = 0; i < 10; i++)
  {
    if (fgets(LineBuf, 128, DataFile) == NULL)
    {
      fclose(DataFile);
      printf("All Data Not There\n");
      exit(0);
    }
    else
    {
      LineBuf[strlen(LineBuf)-1] = 0;
      T_FIL[i] = atof(LineBuf);
    }
  }
  /* Get the Inlet Pressure Must be 1 */
  if (fgets(LineBuf, 128, DataFile) == NULL)
  {
    fclose(DataFile);
    printf("All Data Not There\n");
    exit(0);
  }
  else
  {
    LineBuf[strlen(LineBuf)-1] = 0;
    FIL_P = atof(LineBuf);
  }
  /* Get the 6 L/d values 6 L/d values */
  for (i = 0; i < 6; i++)
  {
    /* Get the L/d value be 6 */
    if (fgets(LineBuf, 128, DataFile) == NULL)
    {
      if (i == 0) return;
      fclose(DataFile);
      printf("All Data Not There\n");
      exit(0);
    }
    else
    {
      LineBuf[strlen(LineBuf)-1] = 0;
      L_Ratio[i] = atof(LineBuf);
    }
  }
  fclose(DataFile);
}

/***********************************************************************/
static double DEGF_to_DEGR(double T)
/*
Purpose:      Convert Fahrneheit to Rankin
\***********************************************************************/
{
  return(T + Rankine);
}

/***********************************************************************/
static double DEGR_to_DEGF(double T)
/*
Purpose:      Convert Rankin to Fahrneheit
\***********************************************************************/
{
  return(T - Rankine);
}

/***********************************************************************/
static double FlowPressure(char *Loc)
/*
Purpose:      Flow Pressure based on Location.
\***********************************************************************/
{
  if (strcmp(Loc, "CIT") == 0)
    return(CIT_P);
  else if (strcmp(Loc, "FEX") == 0)
    return(FEX_P);
  else if (strcmp(Loc, "CDT") == 0)
    return(CDT_P);
  else if (strcmp(Loc, "TIT") == 0)
    return(TIT_P);
  else if (strcmp(Loc, "FIL") == 0)
    return(FIL_P);
  else
    return(0);
}

/***********************************************************************/
static double GetTemp(int  i,
                      char *Loc)
/*
Purpose:      Temperature Value based on Location and index.
\***********************************************************************/
{
  if (strcmp(Loc, "CIT") == 0)
    return(T_CIT[i]);
  else if (strcmp(Loc, "FEX") == 0)
    return(T_FEX[i]);
  else if (strcmp(Loc, "CDT") == 0)
    return(T_CDT[i]);
  else if (strcmp(Loc, "TIT") == 0)
    return(T_TIT[i]);
  else if (strcmp(Loc, "FIL") == 0)
    return(T_FIL[i]);
  else
    return(0);
}

/***********************************************************************/
static double AirDensity(double  Temp,
                         double  Pres)
/*
Purpose:      Compute the Density of Air
\***********************************************************************/
{
  return(AirDensityConstant*Pres/Temp);
}


/***********************************************************************/
static double ReynoldsNumber(double  Vg,
                             double  rho,
                             double  d,
                             double  dav,
                             double  u_dav,
                             double *u)
/*
Purpose:      Compute Reynolds Number
\***********************************************************************/
{
  double  Th1;

 /* Reynolds Number */
  Th1 = -rho*Vg*d/(dav*dav);
  *u = sqrt(Th1*Th1*u_dav*u_dav);
  return(rho*Vg*d/dav);
}

/***********************************************************************/
static double NusseltNumber(   int  Dir,
                            double  Nre,
                            double *u)
/*
Purpose:      Compute Reynolds Number
\***********************************************************************/
{
  /*
  DIR = 0 -> Parallel
  DIR = 1 -> Perpendicular
  */

  /* Perpendicular Orientation? */
  if (Dir)
  {
    *u = 0.06*sqrt(Nre);
    return(0.44*sqrt(Nre));
  }
  /* Parallel Orientation */
  else
  {
    *u = 0.009*pow(Nre,0.674);
    return(0.085*pow(Nre,0.674));
  }
}

/***********************************************************************/
static double SetMachNumber(double  Pt,
                            double  Tj_DEGR,
                            double  MnSet,
                            double *g_set,
                            double *Vg,
                            double *Dp)
/*
Purpose:      Set the Mach Number using inlet pressure and delta pressure
\***********************************************************************/
{
  double Ps;
  double DP;
  double K;
  double g;
  double Tsa;
  double Tta;
  double Mn;
  double MnDiff;
  double Vel;

  Tta = Tj_DEGR;
  g = g_C0 + Tta*(g_C1 + Tta*(g_C2 + Tta*g_C3));
  DP = 0;
  do
  {
    DP += 0.005;
    Ps = Pt - DP;
    K = Pt/Ps;
    Mn = sqrt((pow(K,((g-1)/g))-1)*(2/(g-1)));
    //printf("MnSet = %f Mn = %f\n", MnSet, Mn);
    MnDiff = fabs(MnSet-Mn);
  } while (MnDiff > 0.001);

  Tsa = Tta/(1 + ((g-1)/2)*Mn*Mn);
  Vel = sqrt(GravityConstat*GasConstant*g*Tsa/AirMolecularWeight);
  *Vg = Vel*Mn/2;
  *Dp = DP;
  *g_set = g;
  return(Mn);
}

/***********************************************************************/
static double AirDynamicViscosity(double  T,
                                  double *u)
/*
Purpose:      Compute Dynamic Air Viscosity
\***********************************************************************/
{
  double dav;

 /* Dynamic Air Viscosity */
  if (T < TDAV_MIN)
  {
    *u = U_DAV*DAV_MIN/100;
    return(DAV_MIN);
  }
  if (T > TDAV_MAX)
  {
    *u = U_DAV*DAV_MAX/100;
    return(DAV_MAX);
  }
  dav = DAV_C0 + T*(DAV_C1 + T*DAV_C2);
  *u = U_DAV*dav/100;
  return(dav);
}

/***********************************************************************/
static double AirThermoConductivity(double  T,
                                    double *u)
/*
Purpose:      Compute Air Thermo Conductivity
\***********************************************************************/
{
  double atc;

  /* Air Thermo Conductivity */
  if (T < TATC_MIN)
  {
    *u = U_ATC*ATC_MIN/100;
    return(ATC_MIN);
  }
  if (T > TATC_MAX)
  {
    *u = U_ATC*ATC_MAX/100;
    return(ATC_MAX);
  }
  atc = ATC_C0 + T*(ATC_C1 + T*ATC_C2);
  *u = U_ATC*atc/100;
  return(atc);
}

/***********************************************************************/
static double WireThermoConductivity(double  T,
                                       char  Type,
                                     double *u)
/*
Purpose:      Compute Air Thermo Conductivity
\***********************************************************************/
{
  double wtc;

  /* Wire Thermo Conductivity */
  switch (Type)
  {
    case 'E':
      if (T < TETC_MIN)
      {
        *u = U_ETC*ETC_MIN/100;
        return(ETC_MIN);
      }
      if (T > TETC_MAX)
      {
        *u = U_ETC*ETC_MAX/100;
        return(ETC_MAX);
      }
      wtc = ETC_C0 + T*(ETC_C1 + T*ETC_C2);
      *u = U_ETC*wtc/100;
      return(wtc);
      break;
    case 'J':
      if (T < TKTC_MIN)
      {
        *u = U_KTC*KTC_MIN/100;
        return(JTC_MIN);
      }
      if (T > TJTC_MAX)
      {
        *u = U_JTC*JTC_MAX/100;
        return(JTC_MAX);
      }
      wtc = JTC_C0 + T*(JTC_C1 + T*JTC_C2);
      *u = U_JTC*wtc/100;
      return(wtc);
      break;
    case 'K':
      if (T < TKTC_MIN)
      {
        *u = U_KTC*KTC_MIN/100;
        return(ATC_MIN);
      }
      if (T > TKTC_MAX)
      {
        *u = U_KTC*KTC_MAX/100;
        return(KTC_MAX);
      }
      wtc = KTC_C0 + T*(KTC_C1 + T*KTC_C2);
      *u = U_KTC*wtc/100;
      return(wtc);
      break;
    case 'N':
      if (T < TNTC_MIN)
      {
        *u = U_NTC*NTC_MIN/100;
        return(NTC_MIN);
      }
      if (T > TNTC_MAX)
      {
        *u = U_NTC*NTC_MAX/100;
        return(NTC_MAX);
      }
      wtc = NTC_C0 + T*(NTC_C1 + T*NTC_C2);
      *u = U_NTC*wtc/100;
      return(wtc);
      break;
  }
}

/***********************************************************************/
static double ConvectiveHeatTransferCoefficient(double  Kf,
                                                double  u_Kf,
                                                double  Nnu,
                                                double  u_Nnu,
                                                double  d,
                                                double *u)
/*
Purpose:      Compute Convective Heat Transfer Coefficient
\***********************************************************************/
{
  double Th1;
  double Th2;

 /* Convective Heat Transfer Coefficient */
  Th1 = Kf/d;
  Th2 = Nnu/d;
  *u = sqrt(Th1*Th1*u_Nnu*u_Nnu + Th2*Th2*u_Kf*u_Kf);
  return(Nnu*Kf/d);
}

/***********************************************************************/
static double ConductiveDeltaTemp(double  Tj,
                                  double  w)
/*
Purpose:      Conductive Delta Temperature
\***********************************************************************/
{
  double T_Delta;

  T_Delta =  DEGR_to_DEGF(Tj)*(1 - w);

  return(T_Delta);
}

/***********************************************************************/
static double RadiativeInnerTemp(   int  shield,
                                   char  *Loc,
                                  double  Tj,
                                  double  x,
                                  double  hc,
                                  double  u_hc,
                                  double *DT,
                                  double *u)
/*
Purpose:      Compute Radiative Inner Temperature
\***********************************************************************/
{
    int  i;
  double u_w;
  double A;
  double Th1;
  double Th2;
  double Th3;
  double Th4;
  double T;
  double Tinc;
  double Tend;
  double Twall;
  double DTwall;
  double Tstart;
  double Ti;
  double Ti_In[31];
  double Ti_Eq[31];
  double Ans[4];


  /* Compute Wall Temperature */
  //Twall = x*DEGR_to_DEGF(Tj);
  Twall = x*Tj;
  
  /* Compute Delta Wall Temperature */
  //DTwall = (1-x)*DEGR_to_DEGF(Tj);
  DTwall = (1-x)*Tj;
  
  /* Sheilded Probe? */
  if (!shield)
  {
     *u = 0.0;
     *DT = 0.0;
     Ti = DEGF_to_DEGR(Twall);
     return(Ti);
  }
  if (strcmp(Loc, "CIT") == 0)
  {
    Tstart = T_CIT[0]-200;
    Tend = T_CIT[T_Max-1]+200;
  }
  else if (strcmp(Loc, "FEX") == 0)
  {
    Tstart = T_FEX[0]-200;
    Tend = T_FEX[T_Max-1]+200;
  }
  else if (strcmp(Loc, "CDT") == 0)
  {
    Tstart = T_CIT[0]-200;
    Tend = T_CIT[T_Max-1]+200;
  }
  else if (strcmp(Loc, "TIT") == 0)
  {
    Tstart = T_TIT[0]-200;
    Tend = T_TIT[T_Max-1]+200;
  }
  else if (strcmp(Loc, "FIL") == 0)
  {
    Tstart = T_FIL[0]-200;
    Tend = T_FIL[T_Max-1]+200;
  }
  A = hc/(4.4*StefanBoltzman*Emissivity*DTwall);

  /* Compute Inlet Equation */
  //printf("Ti_In, Ti_Eq\n");
  T = DEGF_to_DEGR(Tstart);
  Tinc = (Tend - Tstart)/30;
  for (i = 0; i < 31; i++)
  {
    T = DEGF_to_DEGR(Tstart);
    Ti_In[i] = T;
    Ti_Eq[i] = T*T*T + A*T - A*Tj;
    //printf("%f\n", Ti_Eq[i]);
    Tstart += Tinc;
  }
  lstsqr(Ti_Eq, Ti_In, 31, 1, Ans);
  Ti = Ans[0];

  Th1 = 1/(4.4*StefanBoltzman*Emissivity*Twall);               /*dA/dhc*/
  Th2 = -hc/(4.4*StefanBoltzman*Emissivity*Emissivity*Twall);  /*dA/de*/
  Th3 = (Tj - Ti)/(3*Ti*Ti + A);                               /*dTi/dA*/

  *u = sqrt(Th1*Th3*u_hc*Th1*Th3*u_hc + Th3*Th2*U_Emissivity*Th3*Th2*U_Emissivity);
  *DT = DTwall; 
  return(Ti);
}

/***********************************************************************/
static double RadiativeInnerTemp1(    int  shield,
                                   char  *Loc,
                                  double  Tj,
                                  double  x,
                                  double  hc,
                                  double  u_hc,
                                  double *u)
/*
Purpose:      Compute Radiative Inner Temperature
\***********************************************************************/
{
    int  i;
  double W_Delta;
  double u_w;
  double A;
  double Th1;
  double Th2;
  double Th3;
  double Th4;
  double T;
  double Tinc;
  double Tend;
  double Tstart;
  double Ti;
  double Ti_In[51];
  double Ti_Eq[51];
  double Ans[4];


  /* Compute Delta Wall Temperature */
  W_Delta = x*Tj;
  
  /* Sheilded Probe? */
  if (!shield)
  {
     *u = 0.0;
     Ti = Tj-W_Delta;
     return(Ti);
  }
  /* Cover full range of Temperature */
  Tstart = -400;
  Tend = 3400;

  A = hc/(4.4*StefanBoltzman*Emissivity*W_Delta);

  /* Compute Inlet Equation */
  Tinc = (Tend - Tstart)/50;
  for (i = 0; i < 51; i++)
  {
    T = DEGF_to_DEGR(Tstart);
    Ti_In[i] = T;
    Ti_Eq[i] = T*T*T + A*T - A*Tj;
    Tstart += Tinc;
  }
  lstsqr(Ti_Eq, Ti_In, 51, 1, Ans);
  Ti = Ans[0];

  Th1 = 1/(4.4*StefanBoltzman*Emissivity*W_Delta);               /*dA/dhc*/
  Th2 = -hc/(4.4*StefanBoltzman*Emissivity*Emissivity*W_Delta);  /*dA/de*/
  Th3 = (Tj - Ti)/(3*Ti*Ti + A);                                 /*dTi/dA*/

  *u = sqrt(Th1*Th3*u_hc*Th1*Th3*u_hc + Th3*Th2*U_Emissivity*Th3*Th2*U_Emissivity);
 
  return(Ti);
}

/***********************************************************************/
static double ConductiveLoss(double  DT,
                             double  L,	
                             double  d,	
                             double  hc,
                             double  ks,
                             double  u_hc,
                             double  u_ks,
                             double *u)
/*
Purpose:      Conductive Temperature Loss
\***********************************************************************/
{
  double X;
  double ThX;
  double ThXhc;
  double ThXks;
  double ThTChc;
  double ThTCks;
  double ThTCDT;
  double rf;
  double Loss;

  /* Conductive Loss */
  Loss = DT/cosh(L*sqrt((4*hc)/(d*ks)));

  /* Uncertainty */
  X = sqrt(4*hc/(d*ks));
  ThXhc = (4/(d*ks))/(2*sqrt(4*hc/(d*ks)));
  ThXks = (-4*hc/(d*ks*ks))/(2*sqrt(4*hc/(d*ks)));

  ThX = -(DT*L*sinh(L*X))/(cosh(L*X)*cosh(L*X));
  ThTChc = ThX*ThXhc;
  ThTCks = ThX*ThXks;

  *u = sqrt((ThTChc*u_hc)*(ThTChc*u_hc) + 
            (ThTCks*u_ks)*(ThTCks*u_ks));

  return(Loss);
}

/***********************************************************************/
static double RadiativeLoss(double  T,
                            double  Ti,	
                            double  hc,
                            double  u_hc,
                            double  u_Ti,
                            double *u)
/*
Purpose:      Radiative Temperature Loss
\***********************************************************************/
{
  double ThTiA;
  double ThTRTi;
  double ThTRe;
  double ThTRhc;
  double Loss;

  /* Radiative Loss */
  Loss = StefanBoltzman*Emissivity*(T*T*T*T - Ti*Ti*Ti*Ti)/hc;

  /* Uncertainty */
  ThTRTi = -4*(StefanBoltzman*Emissivity*Ti*Ti*Ti)/hc;
  ThTRe = StefanBoltzman*(T*T*T*T - Ti*Ti*Ti*Ti)/hc;
  ThTRhc = -StefanBoltzman*Emissivity*(T*T*T*T - Ti*Ti*Ti*Ti)/(hc*hc);

  *u = sqrt((ThTRTi*u_Ti)*(ThTRTi*u_Ti) +
            (ThTRe*U_Emissivity)*(ThTRe*U_Emissivity) +
            (ThTRhc*u_hc)*(ThTRhc*u_hc));

  return(Loss);
}

/***********************************************************************/
static void HTLoss(double  w,
                   double  x,
                   double  d_in,
                      int  Dir,
                     char  Type,
                      int  Shield,
                     char  *Location)
/*
Purpose:      Compute the Heat Losses
\***********************************************************************/
{
     int  i;
     int  j;
     int  l;
    char  s_dir[20];
    char  s_shield[20];
  double  g;
  double  d;
  double  L;
  double  PT;
  double  Vg;
  double  Pd;
  double  DTm;
  double  DTwall;
  double  AD;
  double  ADV;
  double  u_ADV;
  double  Nre;
  double  u_Nre;
  double  Nnu;
  double  u_Nnu;
  double  kf;
  double  u_kf;
  double  ks;
  double  u_ks;
  double  hc;
  double  u_hc;
  double  Ti;
  double  u_Ti;
  double  T;
  double  MN_Set;
  double  MN_Act;
  double  C_Loss;
  double  u_CL;
  double  R_Loss;
  double  u_RL;

  /* Get Inlet Condition */
  d = d_in/12;
  PT = FlowPressure(Location);

  /* Convert to String */
  if (Dir)
    strcpy(s_dir, "Perpendicular");
  else
    strcpy(s_dir, "Parallel");

  if (Shield)
    strcpy(s_shield, "On");
  else
    strcpy(s_shield, "Off");

  for (l = 0; l < L_Max; l++)
  {
    L = L_Ratio[l]*d;
    printf("Condition, Type, Dir, Shield,d,L/d,w,x,Pflow\n");
    printf("%s, %c,%s,%s,%f,%f,%f,%f,%f\n",Location, Type, s_dir, s_shield, d_in, L_Ratio[l], w, x, PT);
    printf("Mach Number, Tj, Ti, DTwall, DTm, Conductive Loss, U_CL, Radiative Loss, U_RL\n");
//    printf("T, MN, DT, AD, ADV, Nre, Nnu, kf, hc, ks, Ti, C_Loss, u_CL, R_Loss, u_RL\n");
    /* Get Temperatures based on location */
    for (i = 0; i < T_Max; i++)
    {
      T = GetTemp(i, Location);
      T = DEGF_to_DEGR(T);
  
      /* Compute Conductive loss at Different MN */
      for (j = 0; j < MN_Max; j++)
      {
        MN_Set = MN_Table[j];

        MN_Act = SetMachNumber(PT, T, MN_Set, &g, &Vg, &Pd);
  
        AD = AirDensity(T, PT);
  
        ADV = AirDynamicViscosity(T, &u_ADV);
  
        Nre = ReynoldsNumber(Vg, AD, d, ADV, u_ADV, &u_Nre);
  
        Nnu = NusseltNumber(Dir, Nre, &u_Nnu);
  
        kf = AirThermoConductivity(T, &u_kf);
  
        hc = ConvectiveHeatTransferCoefficient(kf, u_kf, Nnu, u_Nnu, d, &u_hc);
  
        ks = WireThermoConductivity(T, Type, &u_ks);
  
        DTm = ConductiveDeltaTemp(T, w);

        Ti = RadiativeInnerTemp(Shield, Location, T, x, hc, u_hc, &DTwall, &u_Ti);
  
        C_Loss = ConductiveLoss(DTm, L, d, hc, ks, u_hc, u_ks, &u_CL);
  
        R_Loss = RadiativeLoss(T, Ti, hc, u_hc, u_Ti, &u_RL);
  
//        printf("%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n", DEGR_to_DEGF(T), MN_Act, DT, AD, ADV, Nre, Nnu, kf, hc, ks, DEGR_to_DEGF(Ti), C_Loss, u_CL, R_Loss, u_RL);
        printf("%f,%f,%f,%f,%f,%f,%f,%f,%f\n", MN_Act, DEGR_to_DEGF(T), DEGR_to_DEGF(Ti), DTwall, DTm, C_Loss, u_CL, R_Loss, u_RL);
      }
      printf("\n");
    }
    printf("\n");
  }
}

/***********************************************************************/
int main( int   argc,
         char  *argv[])
/*
Purpose:      Main Program.
\***********************************************************************/
{
     int  i;
     int  Dir;
     int  Shield;
    char  Type;
    char  Location[4];
  double  W;
  double  X;
  double  d;

  if (argc < 7)
  {
    printf("Usage: ht w x d Dir Type Shield Location\n");
    printf("       Output is csv text best to redirect into file and open in spreadsheet\n");
    printf("       w is leakage % value of Tj 0 to 1\n");
    printf("       x is Wall DT % value of Tj 0 to 1\n");
    printf("       d is diameter of T/C wire\n");
    printf("       Dir is orientation of probe in flow - Parallel: 0 Perpendicular: 1, if not forced to 0 or 1\n");
    printf("       Type is T/C type - Only Types allowed are E, J, K, or N\n");
    printf("       Shield is ether ON: 1 or OFF: 0, if not forced to 0 or 1\n");
    printf("       Location is probe location in flow - Only CIT, FEX, CDT, or TIT\n") ;
    exit(1);
  }

  /* Check W is between 0 and 1 */
  W = atof(argv[1]);
  if (W < 0 || W > 1)
  {
    printf("w must be between 0 and 1\n");
    exit(1);
  }

  /* Check x is between 0 and 1 */
  X = atof(argv[2]);
  if (X < 0 || X > 1)
  {
    printf("x must be between 0 and 1\n");
    exit(1);
  }

  /* Check that d is not real large */
  d = atof(argv[3]);
  if (d > 1)
    printf("d seems to be a large number\n");

  /* Make sure Type is usuable */
  Type = argv[5][0];
  switch (Type)
  {
    case 'E':
      break;
    case 'J':
      break;
    case 'K':
      break;
    case 'N':
      break;
    default:
      printf("Type must be E, J, K, or N\n");
      exit(1);
  }

  /* Check Location */
  strcpy(Location, argv[7]);
  Location[4] = 0;
  if (strcmp(Location, "CIT") == 0)
    goto OK;
  else if (strcmp(Location, "FEX") == 0)
    goto OK;
  else if (strcmp(Location, "CDT") == 0)
    goto OK;
  else if (strcmp(Location, "TIT") == 0)
    goto OK;
  else if (strcmp(Location, "FIL") == 0)
    goto OK;
  else
  {
    printf("Location must be CIT, FEX, CDT, TIT, or FIL\n");
    exit(1);
  }
  OK: 
  /* Limit the Boolean inputs */
  Dir = atoi(argv[4]);
  if (Dir < 0)
    Dir = 0;
  if (Dir > 1)
    Dir = 1;
  Shield = atoi(argv[6]);
  if (Shield < 0)
    Shield = 0;
  if (Shield > 1)
    Shield = 1;

  if (strcmp(Location, "FIL") == 0)
    GetTempAndPres();

  /* do 3 ranges of w */
  for (i = 0; i < 3; i++)
  {
    HTLoss(W, X, d, Dir, Type, Shield, Location);
    W += 0.05;
    X += 0.05;
  }
}
