using System;
using UnityEngine;

/*
Copyright (c) 2008, Kiran Pashikanti

Permission to use, copy, modify, and/or distribute this software for any
purpose with or without fee is hereby granted, provided that the above
copyright notice and this permission notice appear in all copies.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR
IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

public static class IAPWS95
{

//public static double IAPWS95_IG_phi(double delta, double tau);
//public static double IAPWS95_IG_phi_delta(double delta, double tau);
//public static double IAPWS95_IG_phi_delta_delta(double delta,double tau);
//public static double IAPWS95_IG_phi_tau(double delta, double tau);
//public static double IAPWS95_IG_phi_tau_tau(double delta,double tau);
//public static double IAPWS95_IG_phi_delta_tau(double delta, double tau);
//public static double IAPWS95_RES_phi(double delta, double tau);
//public static double IAPWS95_RES_phi_delta(double delta, double tau);
//public static double IAPWS95_RES_phi_delta_delta(double delta,double tau);
//public static double IAPWS95_RES_phi_tau(double delta,double tau);
//public static double IAPWS95_RES_phi_tau_tau(double delta,double tau);
//public static double IAPWS95_RES_phi_delta_tau(double delta,double tau);
//
//public static double IAPWS95_pressure(double rho, double T);                                   /*Input: rho in kg/m3, T in K, Output: Pa*/
//public static double IAPWS95_internal_energy(double rho, double T);                            /*Input: rho in kg/m3, T in K, Output: Pa*/
//public static double IAPWS95_entropy(double rho, double T);                                    /*Input: rho in kg/m3, T in K, Output: kJ/kg-K*/
//public static double IAPWS95_enthalpy(double rho, double T);                                   /*Input: rho in kg/m3, T in K, Output: kJ/kg*/
//public static double IAPWS95_isochoric_heat_capacity(double rho, double T);                    /*Input: rho in kg/m3, T in K, Output: kJ/kg-K*/
//public static double IAPWS95_isobaric_heat_capacity(double rho, double T);                     /*Input: rho in kg/m3, T in K, Output: kJ/kg-K*/
//public static double IAPWS95_speed_of_sound(double rho, double T);                             /*Input: rho in kg/m3, T in K, Output: m/s*/
//public static double IAPWS95_joule_thompson_coefficient(double rho, double T);                 /*Input: rho in kg/m3, T in K*/
//public static double IAPWS95_isothermal_throttling_coefficient(double rho, double T);          /*Input: rho in kg/m3, T in K*/
//public static double IAPWS95_isentropic_temperature_pressure_coefficent(double rho, double T); /*Input: rho in kg/m3, T in K*/

/*
Copyright (c) 2008, Kiran Pashikanti

Permission to use, copy, modify, and/or distribute this software for any
purpose with or without fee is hereby granted, provided that the above
copyright notice and this permission notice appear in all copies.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR
IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

const double IAPWS95_CRITICAL_TEMPERATURE  = 647.096;   //#K
const double IAPWS95_CRITICAL_DENSITY    = 322.0;     //#kg m^-3
const double IAPWS95_SPECIFIC_GAS_CONSTANT = 0.46151805;  //#kJ kg^-1 K^-1

static double [] IG_n = {0,-8.32044648201,6.6832105268,3.00632,0.012436,0.97315,1.27950,0.96956,0.24873};
static double [] IG_gamma = {0,0,0,0,1.28728967,3.53734222,7.74073708,9.24437796,27.5075105};

static double [] RES_n =  {0,
      0.12533547935523e-1,  0.78957634722828e+1, -0.87803203303561e+1,
      0.31802509345418e+0, -0.26145533859358e+0, -0.78199751687981e-2,
      0.88089493102134e-2, -0.66856572307965e+0,  0.20433810950965e+0,
     -0.66212605039687e-4, -0.19232721156002e+0, -0.25709043003438e+0,
      0.16074868486251e+0, -0.40092828925807e-1,  0.39343422603254e-6,
     -0.75941377088144e-5,  0.56250979351888e-3, -0.15608652257135e-4,
      0.11537996422951e-8,  0.36582165144204e-6, -0.13251180074668e-11,
     -0.62639586912454e-9, -0.10793600908932e+0,  0.17611491008752e-1,
      0.22132295167546e+0, -0.40247669763528e+0,  0.58083399985759e+0,
      0.49969146990806e-2, -0.31358700712549e-1, -0.74315929710341e+0,
      0.47807329915480e+0,  0.20527940895948e-1, -0.13636435110343e+0,
      0.14180634400617e-1,  0.83326504880713e-2, -0.29052336009585e-1,
      0.38615085574206e-1, -0.20393486513704e-1, -0.16554050063734e-2,
      0.19955571979541e-2,  0.15870308324157e-3, -0.16388568342530e-4,
      0.43613615723811e-1,  0.34994005463765e-1, -0.76788197844621e-1,
      0.22446277332006e-1, -0.62689710414685e-4, -0.55711118565645e-9,
     -0.19905718354408e+0,  0.31777497330738e+0, -0.11841182425981e+0,
     -0.31306260323435e+2,  0.31546140237781e+2, -0.25213154341695e+4,
     -0.14874640856724e+0,  0.31806110878444e+0};

static double [] RES_c = {0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
      1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0,
      2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
      2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 4.0, 6.0, 6.0, 6.0, 6.0};

static double [] RES_d = {0,
      1.0, 1.0, 1.0, 2.0,  2.0,  3.0,  4.0,  1.0,  1.0, 1.0, 2.0,  2.0,  3.0,  4.0,
      4.0, 5.0, 7.0, 9.0, 10.0, 11.0, 13.0, 15.0,  1.0, 2.0, 2.0,  2.0,  3.0,  4.0,
      4.0, 4.0, 5.0, 6.0,  6.0,  7.0,  9.0,  9.0,  9.0, 9.0, 9.0, 10.0, 10.0, 12.0,
      3.0, 4.0, 4.0, 5.0, 14.0,  3.0,  6.0,  6.0,  6.0, 3.0, 3.0,  3.0};

static double [] RES_t = {0,
      -0.5, 0.875,  1.0,  0.5,  0.75, 0.375,  1.0,  4.0,  6.0, 12.0,  1.0,
       5.0, 4.0  ,  2.0, 13.0,  9.0 , 3.0  ,  4.0, 11.0,  4.0, 13.0,  1.0,
       7.0, 1.0  ,  9.0, 10.0, 10.0 , 3.0  ,  7.0, 10.0, 10.0,  6.0, 10.0,
      10.0, 1.0  ,  2.0,  3.0,  4.0 , 8.0  ,  6.0,  9.0,  8.0, 16.0, 22.0,
      23.0,23.0  , 10.0, 50.0, 44.0, 46.0  , 50.0,  0.0,  1.0,  4.0};


static double [] off_RES_alpha = {20.0,20.0,20.0};
static double [] off_RES_beta  = {150.0,150.0,250.0,0.3,0.3};
static double [] off_RES_gamma = {1.21,1.21,1.25};
static double [] off_RES_eps   = {1.0, 1.0 ,1.0 };

static double RES_alpha(int x) { return off_RES_alpha[x-52]; }
static double RES_beta(int x) { return off_RES_beta[x-52]; }
static double RES_gamma(int x) { return off_RES_gamma[x-52]; }
static double RES_eps(int x) { return off_RES_eps[x-52]; }

static double [] off_RES_a = {3.5 ,3.5};
static double [] off_RES_b = {0.85,0.95};
static double [] off_RES_B = {0.2,0.2};
static double [] off_RES_C = {28.0,32.0};
static double [] off_RES_D = {700.0,800.0};
static double [] off_RES_A = {0.32,0.32};

static double RES_a(int x) { return off_RES_a[x-55]; }
static double RES_b(int x) { return off_RES_b[x-55]; }
static double RES_B(int x) { return off_RES_B[x-55]; }
static double RES_C(int x) { return off_RES_C[x-55]; }
static double RES_D(int x) { return off_RES_D[x-55]; }
static double RES_A(int x) { return off_RES_A[x-55]; }

public static double IAPWS95_IG_phi(double delta, double tau)
{
  double t1 = 0.0;
  for(int i=4; i<9; i++)
  {
    t1 += IG_n[i]*Math.Log(1.0 - Math.Exp(-IG_gamma[i]*tau));
  }
  return Math.Log(delta) + IG_n[1] + IG_n[2]*tau + IG_n[3]*Math.Log(tau) + t1;
}

public static double IAPWS95_IG_phi_delta(double delta, double tau)
{
  return 1.0/delta;
}

public static double IAPWS95_IG_phi_delta_delta(double delta,double tau)
{
  return -1.0/(delta*delta);
}

public static double IAPWS95_IG_phi_tau(double delta, double tau)
{
  double t1 = 0.0;
  for(int i=4; i<9; i++)
  {
    t1 += IG_n[i]*IG_gamma[i]*(Math.Pow(1.0 - Math.Exp(-IG_gamma[i]*tau),-1.0) - 1.0);
  }
  return IG_n[2] + IG_n[3]/tau + t1;
}

public static double IAPWS95_IG_phi_tau_tau(double delta,double tau)
{
  double t1 = 0.0;
  for(int i=4; i<9; i++)
  {
    t1 += IG_n[i]*(Math.Pow(IG_gamma[i],2))*Math.Exp(-IG_gamma[i]*tau)*Math.Pow(1.0 - Math.Exp(-IG_gamma[i]*tau),-2.0);
  }
  return -IG_n[3]/(tau*tau) - t1;
}

public static double IAPWS95_IG_phi_delta_tau(double delta, double tau)
{
  return 0.0;
}


public static double IAPWS95_RES_phi(double delta, double tau)
{
  double t1;
  double t2;
  double t3;
  double t4;

  double THETA;
  double DELTA;
  double PSI;

  double p1;
  double p2;

  t1 = 0.0;
  for(int i=1; i<8; i++)
  {
    t1 += RES_n[i]*(Math.Pow(delta,RES_d[i]))*(Math.Pow(tau,RES_t[i]));
  }

  t2 = 0.0;
  for(int i=8; i<52; i++)
  {
    t2 += RES_n[i]*(Math.Pow(delta,RES_d[i]))*(Math.Pow(tau,RES_t[i]))*Math.Exp(-Math.Pow(delta,RES_c[i]));
  }

  t3 = 0.0;
  for(int i=52; i<55; i++)
  {
    p1  = (delta-RES_eps(i))*(delta-RES_eps(i));
    p2  = (tau-RES_gamma(i))*(tau-RES_gamma(i));
    t3 += RES_n[i]*(Math.Pow(delta,RES_d[i]))*(Math.Pow(tau,RES_t[i]))*Math.Exp(-RES_alpha(i)*p1 - RES_beta(i)*p2);
  }

  t4 = 0.0;
  for(int i=55; i<57; i++)
  {
    THETA = (1.0 - tau) + RES_A(i)*Math.Pow((delta - 1.0)*(delta - 1.0),0.5/RES_beta(i));
    DELTA = THETA*THETA + RES_B(i)*Math.Pow((delta - 1.0)*(delta - 1.0),RES_a(i));
    PSI   = Math.Exp(-RES_C(i)*((delta - 1.0)*(delta - 1.0)) - RES_D(i)*((tau - 1.0)*(tau - 1.0)));
    t4 += RES_n[i]*Math.Pow(DELTA,RES_b(i))*delta*PSI;
  }
  return t1 + t2 + t3 + t4;
}


public static double IAPWS95_RES_phi_delta(double delta, double tau)
{
  double t1;
  double t2;
  double t3;
  double t4;

  double THETA;
  double DELTA;
  double D_DELTA_delta;
  double D_DELTA_bi_delta;
  double PSI;
  double D_PSI_delta;

  double p1;
  double p2;

  int i;

  t1 = 0.0;
  for(i=1; i<8; i++)
  {
    t1 += RES_n[i]*RES_d[i]*Math.Pow(delta,RES_d[i]-1.0)*Math.Pow(tau,RES_t[i]);
  }

  t2 = 0.0;
  for(i=8; i<52; i++)
  {
    t2 += RES_n[i]*Math.Exp(-Math.Pow(delta,RES_c[i]))*(Math.Pow(delta,RES_d[i]-1.0)*Math.Pow(tau,RES_t[i])*(RES_d[i] - RES_c[i]*Math.Pow(delta,RES_c[i])));
  }

  t3 = 0.0;
  for(i=52; i<55; i++)
  {
    p1 = Math.Exp(-RES_alpha(i)*((delta - RES_eps(i))*(delta - RES_eps(i))) - RES_beta(i)*((tau - RES_gamma(i))*(tau - RES_gamma(i))));
    p2 = (RES_d[i]/delta) - 2.0*RES_alpha(i)*(delta - RES_eps(i));
    t3 += RES_n[i]*Math.Pow(delta,RES_d[i])*Math.Pow(tau,RES_t[i])*p1*p2;
  }

   t4 = 0.0;
   for (i=55; i<57; i++)
   {
    THETA = (1.0 - tau) + RES_A(i)*Math.Pow((delta - 1.0)*(delta - 1.0),0.5/RES_beta(i));
    DELTA = THETA*THETA + RES_B(i)*Math.Pow((delta - 1.0)*(delta - 1.0),RES_a(i));
    PSI   = Math.Exp(-RES_C(i)*((delta - 1.0)*(delta - 1.0)) - RES_D(i)*((tau - 1.0)*(tau - 1.0)));

    D_DELTA_delta = (delta - 1.0)*(RES_A(i)*THETA*(2.0/RES_beta(i))*Math.Pow((delta - 1.0)*(delta - 1.0),0.5/RES_beta(i) - 1.0) + 2.0*RES_B(i)*RES_a(i)*Math.Pow((delta - 1.0)*(delta - 1.0),RES_a(i)-1.0));
    D_DELTA_bi_delta = RES_b(i)*Math.Pow(DELTA,RES_b(i)-1.0)*D_DELTA_delta;
    D_PSI_delta = -2.0*RES_C(i)*(delta - 1.0)*PSI;

    t4 += RES_n[i]*(Math.Pow(DELTA,RES_b(i))*(PSI + delta*D_PSI_delta) + D_DELTA_bi_delta*delta*PSI);
   }

   return t1 + t2 + t3 + t4;
}



public static double IAPWS95_RES_phi_delta_delta(double delta,double tau)
{
  double t1;
  double t2;
  double t3;
  double t4;

  double THETA;
  double DELTA;
  double D_DELTA_delta;
  double D_DELTA_bi_delta;
  double D2_DELTA_delta;
  double D2_DELTA_bi_delta;
  double PSI;
  double D_PSI_delta;
  double D2_PSI_delta;

  double p1;
  double p2;
  double p3;
  double p4;
  double p5;

  int i;

  t1 = 0.0;
  for(i=1; i<8; i++)
  {
    t1 += RES_n[i]*RES_d[i]*(RES_d[i]-1.0)*Math.Pow(delta,RES_d[i]-2.0)*Math.Pow(tau,RES_t[i]);
  }
  t2 = 0.0;
  for(i=8; i<52; i++)
  {
    p1 = Math.Pow(delta,RES_d[i] - 2.0)*Math.Pow(tau,RES_t[i]);
    p2 = (RES_d[i] - RES_c[i]*Math.Pow(delta,RES_c[i]))*(RES_d[i] - 1.0 - (RES_c[i])*Math.Pow(delta,RES_c[i]));
    t2 += RES_n[i]*Math.Exp(-Math.Pow(delta,RES_c[i]))*(p1*(p2 - (RES_c[i]*RES_c[i])*Math.Pow(delta,RES_c[i])));
  }
  t3 = 0.0;
  for(i=52; i<55; i++)
  {
    p1 = Math.Exp(-RES_alpha(i)*((delta - RES_eps(i))*(delta - RES_eps(i))) - RES_beta(i)*((tau-RES_gamma(i))*(tau-RES_gamma(i))));
    p2 = 2.0*RES_alpha(i)*Math.Pow(delta,RES_d[i]);
    p3 = 4.0*(RES_alpha(i)*RES_alpha(i))*Math.Pow(delta,RES_d[i])*((delta - RES_eps(i))*(delta - RES_eps(i)));
    p4 = 4.0*RES_d[i]*RES_alpha(i)*Math.Pow(delta,RES_d[i] - 1.0)*(delta - RES_eps(i));
    p5 = RES_d[i]*(RES_d[i] - 1.0)*Math.Pow(delta,RES_d[i] - 2.0);
    t3 += RES_n[i]*Math.Pow(tau,RES_t[i])*p1*(-p2 + p3 - p4 + p5);
  }
  t4 = 0.0;
  for (i=55; i<57; i++)
  {
    THETA = (1.0 - tau) + RES_A(i)*Math.Pow((delta - 1.0)*(delta - 1.0),0.5/RES_beta(i));
    DELTA = THETA*THETA + RES_B(i)*Math.Pow((delta - 1.0)*(delta - 1.0),RES_a(i));
    PSI   = Math.Exp(-RES_C(i)*((delta - 1.0)*(delta - 1.0)) - RES_D(i)*((tau - 1.0)*(tau - 1.0)));

    D_DELTA_delta = (delta - 1.0)*(RES_A(i)*THETA*(2.0/RES_beta(i))*Math.Pow((delta - 1.0)*(delta - 1.0),0.5/RES_beta(i) - 1.0) + 2.0*RES_B(i)*RES_a(i)*Math.Pow((delta - 1.0)*(delta - 1.0),RES_a(i)-1.0));
    D_DELTA_bi_delta = RES_b(i)*Math.Pow(DELTA,RES_b(i)-1.0)*D_DELTA_delta;
    D_PSI_delta = -2.0*RES_C(i)*(delta - 1.0)*PSI;

    D2_PSI_delta = 2.0*RES_C(i)*PSI*(2.0*RES_C(i)*(delta - 1.0)*(delta - 1.0) - 1.0);

    p1 = 4.0*RES_B(i)*RES_a(i)*(RES_a(i) - 1.0)*Math.Pow((delta - 1.0)*(delta - 1.0),RES_a(i) - 2.0);
    p2 = 2.0*(RES_A(i)*RES_A(i))*Math.Pow(RES_beta(i),-2.0)*Math.Pow(Math.Pow((delta - 1.0)*(delta - 1.0),0.5/RES_beta(i) - 1.0),2);
    p3 = RES_A(i)*THETA*(4.0/RES_beta(i))*(0.5/RES_beta(i) - 1.0)*Math.Pow((delta - 1.0)*(delta - 1.0),0.5/RES_beta(i) - 2.0);

    D2_DELTA_delta = (D_DELTA_delta/(delta - 1.0)) + ((delta - 1.0)*(delta - 1.0))*(p1 + p2 + p3);
    D2_DELTA_bi_delta = RES_b(i)*(Math.Pow(DELTA,RES_b(i) - 1.0)*D2_DELTA_delta + (RES_b(i) - 1.0)*Math.Pow(DELTA,RES_b(i) - 2.0)*Math.Pow(D_DELTA_delta,2.0));

    p1 = Math.Pow(DELTA,RES_b(i))*(2.0*D_PSI_delta + delta*D2_PSI_delta);
    p2 = 2.0*(D_DELTA_bi_delta)*(PSI + delta*D_PSI_delta);
    p3 = D2_DELTA_bi_delta*delta*PSI;

    t4 += RES_n[i]*(p1 + p2 + p3);
  }

  return t1 + t2 + t3 + t4;
}


public static double IAPWS95_RES_phi_tau(double delta,double tau)
{
  double t1;
  double t2;
  double t3;
  double t4;

  double THETA;
  double DELTA;
  double D_DELTA_bi_tau;
  double PSI;
  double D_PSI_tau;

  double p1;
  double p2;

  int i;

  t1 = 0.0;
  for(i=1; i<8; i++)
  {
    t1 += RES_n[i]*RES_t[i]*Math.Pow(delta,RES_d[i])*Math.Pow(tau,RES_t[i] - 1.0);
  }
  t2 = 0.0;
  for(i=8; i<52; i++)
  {
    t2 += RES_n[i]*RES_t[i]*Math.Pow(delta,RES_d[i])*Math.Pow(tau,RES_t[i] - 1.0)*Math.Exp(-Math.Pow(delta,RES_c[i]));
  }
  t3 = 0.0;
  for(i=52; i<55; i++)
  {
    p1 = -RES_alpha(i)*((delta - RES_eps(i))*(delta - RES_eps(i))) - RES_beta(i)*((tau - RES_gamma(i))*(tau - RES_gamma(i)));
    p2 = RES_t[i]/tau - 2.0*RES_beta(i)*(tau - RES_gamma(i));
    t3 += RES_n[i]*Math.Pow(delta,RES_d[i])*Math.Pow(tau,RES_t[i])*Math.Exp(p1)*p2;
  }
  t4 = 0.0;
  for (i=55; i<57; i++)
  {
    THETA = (1.0 - tau) + RES_A(i)*Math.Pow((delta - 1.0)*(delta - 1.0),0.5/RES_beta(i));
    DELTA = THETA*THETA + RES_B(i)*Math.Pow((delta - 1.0)*(delta - 1.0),RES_a(i));
    PSI   = Math.Exp(-RES_C(i)*((delta - 1.0)*(delta - 1.0)) - RES_D(i)*((tau - 1.0)*(tau - 1.0)));

    D_DELTA_bi_tau = -2.0*THETA*RES_b(i)*Math.Pow(DELTA,RES_b(i) - 1.0);
    D_PSI_tau = -2.0*RES_D(i)*(tau - 1.0)*PSI;
    t4 += RES_n[i]*delta*(D_DELTA_bi_tau*PSI + Math.Pow(DELTA,RES_b(i))*D_PSI_tau);

  }
  return t1 + t2 + t3 + t4;
}

public static double IAPWS95_RES_phi_tau_tau(double delta,double tau)
{
  double t1;
  double t2;
  double t3;
  double t4;

  double THETA;
  double DELTA;
  double D_DELTA_bi_tau;
  double D2_DELTA_bi_tau;
  double PSI;
  double D_PSI_tau;
  double D2_PSI_tau;

  double p1;
  double p2;

  int i;

  t1 = 0.0;
  for(i=1; i<8; i++)
  {
    t1 += RES_n[i]*RES_t[i]*(RES_t[i] - 1.0)*Math.Pow(delta,RES_d[i])*Math.Pow(tau,RES_t[i] - 2.0);
  }
  t2 = 0.0;
  for(i=8; i<52; i++)
  {
    t2 += RES_n[i]*RES_t[i]*(RES_t[i] - 1.0)*Math.Pow(delta,RES_d[i])*Math.Pow(tau,RES_t[i] - 2.0)*Math.Exp(-Math.Pow(delta,RES_c[i]));
  }
  t3 = 0.0;
  for(i=52; i<55; i++)
  {
    p1 = -RES_alpha(i)*((delta - RES_eps(i))*(delta - RES_eps(i))) - RES_beta(i)*((tau - RES_gamma(i))*(tau - RES_gamma(i)));
    p2 = Math.Pow((RES_t[i]/tau - 2.0*RES_beta(i)*(tau - RES_gamma(i))),2.0) - RES_t[i]/(tau*tau) - 2.0*RES_beta(i);
    t3 += RES_n[i]*Math.Pow(delta,RES_d[i])*Math.Pow(tau,RES_t[i])*Math.Exp(p1)*p2;
  }
  t4 = 0.0;
  for (i=55; i<57; i++)
  {
    THETA = (1.0 - tau) + RES_A(i)*Math.Pow((delta - 1.0)*(delta - 1.0),0.5/RES_beta(i));
    DELTA = THETA*THETA + RES_B(i)*Math.Pow((delta - 1.0)*(delta - 1.0),RES_a(i));
    PSI   = Math.Exp(-RES_C(i)*((delta - 1.0)*(delta - 1.0)) - RES_D(i)*((tau - 1.0)*(tau - 1.0)));

    D_DELTA_bi_tau = -2.0*THETA*RES_b(i)*Math.Pow(DELTA,RES_b(i) - 1.0);
    D_PSI_tau = -2.0*RES_D(i)*(tau - 1.0)*PSI;

    D2_DELTA_bi_tau = 2.0*RES_b(i)*Math.Pow(DELTA,RES_b(i) - 1.0) + 4.0*(THETA*THETA)*RES_b(i)*(RES_b(i) - 1.0)*Math.Pow(DELTA,RES_b(i) - 2.0);
    D2_PSI_tau = (2*RES_D(i)*((tau - 1.0)*(tau - 1.0)) - 1.0)*2.0*RES_D(i)*PSI;

    t4 += RES_n[i]*delta*(D2_DELTA_bi_tau*PSI + 2.0*D_DELTA_bi_tau*D_PSI_tau + Math.Pow(DELTA,RES_b(i))*D2_PSI_tau);
  }
  return t1 + t2 + t3 + t4;
}

public static double IAPWS95_RES_phi_delta_tau(double delta,double tau)
{
  double t1;
  double t2;
  double t3;
  double t4;

  double THETA;
  double DELTA;
  double D_DELTA_delta;
  double D_DELTA_bi_delta;
  double D_DELTA_bi_tau;
  double PSI;
  double D_PSI_delta;
  double D_PSI_tau;
  double D2_DELTA_bi_tau;
  double D2_PSI_tau;
  double D2_DELTA_bi_delta_tau;
  double D2_PSI_delta_tau;

  double p1;
  double p2;
  double p3;
  double p4;

  int i;
  t1 = 0.0;
  for(i=1; i<8; i++)
  {
    t1 += RES_n[i]*RES_d[i]*RES_t[i]*Math.Pow(delta,RES_d[i] - 1.0)*Math.Pow(tau,RES_t[i] - 1.0);
  }
  t2 = 0.0;
  for(i=8; i<52; i++)
  {
    p1  = (RES_d[i] - RES_c[i]*Math.Pow(delta,RES_c[i]))*Math.Exp(-Math.Pow(delta,RES_c[i]));
    t2 += RES_n[i]*RES_t[i]*Math.Pow(delta,RES_d[i] - 1.0)*Math.Pow(tau,RES_t[i] - 1.0)*p1;
  }
  t3 = 0.0;
  for(i=52; i<55; i++)
  {
    p1 = -RES_alpha(i)*((delta - RES_eps(i))*(delta - RES_eps(i))) - RES_beta(i)*((tau - RES_gamma(i))*(tau - RES_gamma(i)));
    p2 = (-RES_d[i]/delta - 2.0*RES_alpha(i)*(delta - RES_eps(i)))*(RES_t[i]/tau - 2.0*RES_beta(i)*(tau - RES_gamma(i)));
    t3 += RES_n[i]*Math.Pow(delta,RES_d[i])*Math.Pow(tau,RES_t[i])*Math.Exp(p1)*p2;
  }
  t4 = 0.0;
  for (i=55; i<57; i++)
  {
    THETA = (1.0 - tau) + RES_A(i)*Math.Pow((delta - 1.0)*(delta - 1.0),0.5/RES_beta(i));
    DELTA = THETA*THETA + RES_B(i)*Math.Pow((delta - 1.0)*(delta - 1.0),RES_a(i));
    PSI   = Math.Exp(-RES_C(i)*((delta - 1.0)*(delta - 1.0)) - RES_D(i)*((tau - 1.0)*(tau - 1.0)));

    D_DELTA_bi_tau = -2.0*THETA*RES_b(i)*Math.Pow(DELTA,RES_b(i) - 1.0);
    D_PSI_tau = -2.0*RES_D(i)*(tau - 1.0)*PSI;

    D2_DELTA_bi_tau = 2.0*RES_b(i)*Math.Pow(DELTA,RES_b(i) - 1.0) + 4.0*(THETA*THETA)*RES_b(i)*(RES_b(i) - 1.0)*Math.Pow(DELTA,RES_b(i) - 2.0);
    D2_PSI_tau = (2*RES_D(i)*((tau - 1.0)*(tau - 1.0)) - 1.0)*2.0*RES_D(i)*PSI;

    D_DELTA_delta = (delta - 1.0)*(RES_A(i)*THETA*(2.0/RES_beta(i))*Math.Pow((delta - 1.0)*(delta - 1.0),0.5/RES_beta(i) - 1.0) + 2.0*RES_B(i)*RES_a(i)*Math.Pow((delta - 1.0)*(delta - 1.0),RES_a(i)-1.0));
    D_DELTA_bi_delta = RES_b(i)*Math.Pow(DELTA,RES_b(i)-1.0)*D_DELTA_delta;
    D_PSI_delta = -2.0*RES_C(i)*(delta - 1.0)*PSI;

    D2_PSI_delta_tau = 4.0*RES_C(i)*RES_D(i)*(delta - 1.0)*(tau - 1.0)*PSI;

    p1 = RES_A(i)*RES_b(i)*(2.0/RES_beta(i))*Math.Pow(DELTA,RES_b(i) - 1.0)*(delta - 1.0)*Math.Pow((delta - 1.0)*(delta - 1.0),0.5/RES_beta(i) - 1.0);
    p2 = 2.0*THETA*RES_b(i)*(RES_b(i) - 1.0)*Math.Pow(DELTA,RES_b(i) - 1.0)*D_DELTA_delta;
    D2_DELTA_bi_delta_tau = -p1 - p2;

    p1 = Math.Pow(DELTA,RES_b(i))*(D_PSI_tau + delta*D2_PSI_delta_tau);
    p2 = delta*D_DELTA_bi_tau*D_PSI_tau;
    p3 = D_DELTA_bi_tau*(PSI + delta*D_PSI_delta);
    p4 = D2_DELTA_bi_delta_tau*delta*PSI;

    t4 += RES_n[i]*(p1 + p2 + p3 + p4);
  }
  return t1 + t2 + t3 + t4;
}


public static double IAPWS95_pressure(double rho, double T)
{
  double delta;
  double tau;

  delta = rho/IAPWS95_CRITICAL_DENSITY;
  tau = IAPWS95_CRITICAL_TEMPERATURE/T;
  return rho*IAPWS95_SPECIFIC_GAS_CONSTANT*T*(1.0 + delta*IAPWS95_RES_phi_delta(delta,tau));
}

public static double IAPWS95_internal_energy(double rho, double T)
{
  double delta;
  double tau;

  delta = rho/IAPWS95_CRITICAL_DENSITY;
  tau = IAPWS95_CRITICAL_TEMPERATURE/T;
  return tau*IAPWS95_SPECIFIC_GAS_CONSTANT*T*(IAPWS95_IG_phi_tau(delta,tau) + IAPWS95_RES_phi_tau(delta,tau));
}

public static double IAPWS95_entropy(double rho, double T)
{
  double delta;
  double tau;

  delta = rho/IAPWS95_CRITICAL_DENSITY;
  tau = IAPWS95_CRITICAL_TEMPERATURE/T;
  return IAPWS95_SPECIFIC_GAS_CONSTANT*(tau*(IAPWS95_IG_phi_tau(delta,tau) + IAPWS95_RES_phi_tau(delta,tau)) - IAPWS95_IG_phi(delta,tau) - IAPWS95_RES_phi(delta,tau));
}

public static double IAPWS95_enthalpy(double rho, double T)
{
  double delta;
  double tau;

  delta = rho/IAPWS95_CRITICAL_DENSITY;
  tau = IAPWS95_CRITICAL_TEMPERATURE/T;
  return (IAPWS95_SPECIFIC_GAS_CONSTANT*T)*(1.0 + tau*(IAPWS95_IG_phi_tau(delta,tau) + IAPWS95_RES_phi_tau(delta,tau)) + delta*IAPWS95_RES_phi_delta(delta,tau));
}

public static double IAPWS95_isochoric_heat_capacity(double rho, double T)
{
  double delta;
  double tau;

  delta = rho/IAPWS95_CRITICAL_DENSITY;
  tau = IAPWS95_CRITICAL_TEMPERATURE/T;
  return (-(tau*tau)*(IAPWS95_IG_phi_tau_tau(delta,tau) + IAPWS95_RES_phi_tau_tau(delta,tau)))*IAPWS95_SPECIFIC_GAS_CONSTANT;
}

public static double IAPWS95_isobaric_heat_capacity(double rho, double T)
{
  double delta;
  double tau;

  double p1;
  double p2;
  double p3;

  delta = rho/IAPWS95_CRITICAL_DENSITY;
  tau = IAPWS95_CRITICAL_TEMPERATURE/T;
  p1 = -(tau*tau)*(IAPWS95_IG_phi_tau_tau(delta,tau) + IAPWS95_RES_phi_tau_tau(delta,tau));
  p2 = Math.Pow((1.0 + delta*IAPWS95_RES_phi_delta(delta,tau) - delta*tau*IAPWS95_RES_phi_delta_tau(delta,tau)),2.0);
  p3 = 1.0 + 2.0*delta*IAPWS95_RES_phi_delta(delta,tau) + (delta*delta)*IAPWS95_RES_phi_delta_delta(delta,tau);
  return (p1 + p2/p3)*IAPWS95_SPECIFIC_GAS_CONSTANT;
}

public static double IAPWS95_speed_of_sound(double rho, double T)
{
  double delta;
  double tau;

  double p1;
  double p2;
  double p3;

  delta = rho/IAPWS95_CRITICAL_DENSITY;
  tau = IAPWS95_CRITICAL_TEMPERATURE/T;
  p1 = 1.0 + 2.0*delta*IAPWS95_RES_phi_delta(delta,tau) + (delta*delta)*IAPWS95_RES_phi_delta_delta(delta,tau);
  p2 = Math.Pow((1.0 + delta*IAPWS95_RES_phi_delta(delta,tau) - delta*tau*IAPWS95_RES_phi_delta_tau(delta,tau)),2.0);
  p3 = (tau*tau)*(IAPWS95_IG_phi_tau_tau(delta,tau) + IAPWS95_RES_phi_tau_tau(delta,tau));
  return Math.Sqrt((IAPWS95_SPECIFIC_GAS_CONSTANT*T)*(p1 - p2/p3)*1000.0);
}

public static double IAPWS95_joule_thompson_coefficient(double rho, double T)
{
  double delta;
  double tau;

  double p1;
  double p2;
  double p3;
  double p4;

  delta = rho/IAPWS95_CRITICAL_DENSITY;
  tau = IAPWS95_CRITICAL_TEMPERATURE/T;
  p1 = delta*IAPWS95_RES_phi_delta(delta,tau) + (delta*delta)*IAPWS95_RES_phi_delta_delta(delta,tau) + (delta*tau)*IAPWS95_RES_phi_delta_tau(delta,tau);
  p2 = Math.Pow(1.0 + delta*IAPWS95_RES_phi_delta(delta,tau) - delta*tau*IAPWS95_RES_phi_delta_tau(delta,tau),2.0);
  p3 = IAPWS95_IG_phi_tau_tau(delta,tau) + IAPWS95_RES_phi_tau_tau(delta,tau);
  p4 = 1.0 + 2.0*delta*IAPWS95_RES_phi_delta(delta,tau) + (delta*delta)*IAPWS95_RES_phi_delta_delta(delta,tau);
  return (-p1/(p2 - (tau*tau)*p3*p4))/(IAPWS95_SPECIFIC_GAS_CONSTANT*rho);
}

public static double IAPWS95_isothermal_throttling_coefficient(double rho, double T)
{
  double delta;
  double tau;

  double p1;
  double p2;

  delta = rho/IAPWS95_CRITICAL_DENSITY;
  tau = IAPWS95_CRITICAL_TEMPERATURE/T;
  p1 = 1.0 + delta*IAPWS95_RES_phi_delta(delta,tau)-delta*tau*IAPWS95_RES_phi_delta_tau(delta,tau);
  p2 = 1.0 + 2.0*delta*IAPWS95_RES_phi_delta(delta,tau) + (delta*delta)*IAPWS95_RES_phi_delta_delta(delta,tau);
  return (1.0 - p1/p2)/rho;
}

public static double IAPWS95_isentropic_temperature_pressure_coefficent(double rho, double T)
{
  double delta;
  double tau;

  double p1;
  double p2;
  double p3;
  double p4;

  delta = rho/IAPWS95_CRITICAL_DENSITY;
  tau = IAPWS95_CRITICAL_TEMPERATURE/T;
  p1 = 1.0 + delta*IAPWS95_RES_phi_delta(delta,tau) - delta*tau*IAPWS95_RES_phi_delta_tau(delta,tau);
  p2 = Math.Pow(1.0 + delta*IAPWS95_RES_phi_delta(delta,tau) - delta*tau*IAPWS95_RES_phi_delta_tau(delta,tau),2.0);
  p3 = IAPWS95_IG_phi_tau_tau(delta,tau) + IAPWS95_RES_phi_tau_tau(delta,tau);
  p4 = 1.0 + 2.0*delta*IAPWS95_RES_phi_delta(delta,tau) + (delta*delta)*IAPWS95_RES_phi_delta_delta(delta,tau);

  return (p1/(p2 - (tau*tau)*p3*p4))/(IAPWS95_SPECIFIC_GAS_CONSTANT*rho);
}

public static void print_tables()
{
  double delta = 838.025/IAPWS95_CRITICAL_DENSITY;
  double tau = IAPWS95_CRITICAL_TEMPERATURE/500.0;
  double T;
  double rho;

  double [] Ta = {275.0,450.0,625.0};
  double [] rhoL = {0.999887406e3,0.890341250e3,0.567090385e3};
  double [] rhoV = {0.550664919e-2,0.481200360e1,0.118290280e3};

  Debug.Log("IAPWS95: TABLE 6 (IDEAL AND RESIDUAL)\n");
  Debug.Log("-------------------------------\n");
  Debug.LogFormat("PHI              {0,9:E} \t {1,9:E}\n",IAPWS95_IG_phi(delta,tau),IAPWS95_RES_phi(delta,tau));
  Debug.LogFormat("PHI_delta        {0,9:E} \t {1,9:E}\n",IAPWS95_IG_phi_delta(delta,tau),IAPWS95_RES_phi_delta(delta,tau));
  Debug.LogFormat("PHI_delta_delta  {0,9:E} \t {1,9:E}\n",IAPWS95_IG_phi_delta_delta(delta,tau),IAPWS95_RES_phi_delta_delta(delta,tau));
  Debug.LogFormat("PHI_tau          {0,9:E} \t {1,9:E}\n",IAPWS95_IG_phi_tau(delta,tau),IAPWS95_RES_phi_tau(delta,tau));
  Debug.LogFormat("PHI_tau_tau      {0,9:E} \t {1,9:E}\n",IAPWS95_IG_phi_tau_tau(delta,tau),IAPWS95_RES_phi_tau_tau(delta,tau));
  Debug.LogFormat("PHI_delta_tau    {0,9:E} \t {1,9:E}\n",IAPWS95_IG_phi_delta_tau(delta,tau),IAPWS95_RES_phi_delta_tau(delta,tau));


  Debug.Log("IAPWS95: TABLE 7 (VERIFICATION)\n");
  Debug.Log("-------------------------------\n");
  T = 300.0;
  rho = 0.996556e3;
  Debug.LogFormat("{0} {1,9:E} {2,9:E} {3,9:E} {4,9:E} {5,9:E}\n", T,rho,IAPWS95_pressure(rho,T)*1e-3,IAPWS95_isochoric_heat_capacity(rho,T),IAPWS95_speed_of_sound(rho,T),IAPWS95_entropy(rho,T));
  rho = 0.1005308e4;
  Debug.LogFormat("{0} {1,9:E} {2,9:E} {3,9:E} {4,9:E} {5,9:E}\n", T,rho,IAPWS95_pressure(rho,T)*1e-3,IAPWS95_isochoric_heat_capacity(rho,T),IAPWS95_speed_of_sound(rho,T),IAPWS95_entropy(rho,T));
  rho = 0.1188202e4;
  Debug.LogFormat("{0} {1,9:E} {2,9:E} {3,9:E} {4,9:E} {5,9:E}\n", T,rho,IAPWS95_pressure(rho,T)*1e-3,IAPWS95_isochoric_heat_capacity(rho,T),IAPWS95_speed_of_sound(rho,T),IAPWS95_entropy(rho,T));

  Debug.Log("\n\n");

  T = 500.0;
  rho = 0.435;
  Debug.LogFormat("{0} {1,9:E} {2,9:E} {3,9:E} {4,9:E} {5,9:E}\n", T,rho,IAPWS95_pressure(rho,T)*1e-3,IAPWS95_isochoric_heat_capacity(rho,T),IAPWS95_speed_of_sound(rho,T),IAPWS95_entropy(rho,T));
  rho = 0.4352e1;
  Debug.LogFormat("{0} {1,9:E} {2,9:E} {3,9:E} {4,9:E} {5,9:E}\n", T,rho,IAPWS95_pressure(rho,T)*1e-3,IAPWS95_isochoric_heat_capacity(rho,T),IAPWS95_speed_of_sound(rho,T),IAPWS95_entropy(rho,T));
  rho = 0.838025e3;
  Debug.LogFormat("{0} {1,9:E} {2,9:E} {3,9:E} {4,9:E} {5,9:E}\n", T,rho,IAPWS95_pressure(rho,T)*1e-3,IAPWS95_isochoric_heat_capacity(rho,T),IAPWS95_speed_of_sound(rho,T),IAPWS95_entropy(rho,T));
  rho = 0.1084564e4;
  Debug.LogFormat("{0} {1,9:E} {2,9:E} {3,9:E} {4,9:E} {5,9:E}\n", T,rho,IAPWS95_pressure(rho,T)*1e-3,IAPWS95_isochoric_heat_capacity(rho,T),IAPWS95_speed_of_sound(rho,T),IAPWS95_entropy(rho,T));

  Debug.Log("\n\n");

  T = 647.0;
  rho = 0.358e3;
  Debug.LogFormat("{0} {1,9:E} {2,9:E} {3,9:E} {4,9:E} {5,9:E}\n", T,rho,IAPWS95_pressure(rho,T)*1e-3,IAPWS95_isochoric_heat_capacity(rho,T),IAPWS95_speed_of_sound(rho,T),IAPWS95_entropy(rho,T));

  Debug.Log("\n\n");

  T = 900.0;
  rho = 0.241;
  Debug.LogFormat("{0} {1,9:E} {2,9:E} {3,9:E} {4,9:E} {5,9:E}\n", T,rho,IAPWS95_pressure(rho,T)*1e-3,IAPWS95_isochoric_heat_capacity(rho,T),IAPWS95_speed_of_sound(rho,T),IAPWS95_entropy(rho,T));
  rho = 0.526150e2;
  Debug.LogFormat("{0} {1,9:E} {2,9:E} {3,9:E} {4,9:E} {5,9:E}\n", T,rho,IAPWS95_pressure(rho,T)*1e-3,IAPWS95_isochoric_heat_capacity(rho,T),IAPWS95_speed_of_sound(rho,T),IAPWS95_entropy(rho,T));
  rho = 0.870769e3;
  Debug.LogFormat("{0} {1,9:E} {2,9:E} {3,9:E} {4,9:E} {5,9:E}\n", T,rho,IAPWS95_pressure(rho,T)*1e-3,IAPWS95_isochoric_heat_capacity(rho,T),IAPWS95_speed_of_sound(rho,T),IAPWS95_entropy(rho,T));
  Debug.LogFormat("\n\n");

  Debug.Log("IAPWS95: TABLE 8 (2-PHASE)\n");
  Debug.Log("--------------------------\n");
  Debug.LogFormat("{0,9:E} {1,9:E} {2,9:E}\n", IAPWS95_pressure(rhoL[0],Ta[0]),IAPWS95_pressure(rhoL[1],Ta[1]),IAPWS95_pressure(rhoL[2],Ta[2]));
  Debug.LogFormat("{0,9:E} {1,9:E} {2,9:E}\n", IAPWS95_pressure(rhoV[0],Ta[0]),IAPWS95_pressure(rhoV[1],Ta[1]),IAPWS95_pressure(rhoV[2],Ta[2]));
  Debug.LogFormat("{0,9:E} {1,9:E} {2,9:E}\n", rhoL[0], rhoL[1], rhoL[2]); 
  Debug.LogFormat("{0,9:E} {1,9:E} {2,9:E}\n", rhoV[0], rhoV[1], rhoV[2]);
  Debug.LogFormat("{0,9:E} {1,9:E} {2,9:E}\n", IAPWS95_enthalpy(rhoL[0],Ta[0]),IAPWS95_enthalpy(rhoL[1],Ta[1]),IAPWS95_enthalpy(rhoL[2],Ta[2]));
  Debug.LogFormat("{0,9:E} {1,9:E} {2,9:E}\n", IAPWS95_enthalpy(rhoV[0],Ta[0]),IAPWS95_enthalpy(rhoV[1],Ta[1]),IAPWS95_enthalpy(rhoV[2],Ta[2]));
  Debug.LogFormat("{0,9:E} {1,9:E} {2,9:E}\n", IAPWS95_entropy(rhoL[0],Ta[0]),IAPWS95_entropy(rhoL[1],Ta[1]),IAPWS95_entropy(rhoL[2],Ta[2]));
  Debug.LogFormat("{0,9:E} {1,9:E} {2,9:E}\n", IAPWS95_entropy(rhoV[0],Ta[0]),IAPWS95_entropy(rhoV[1],Ta[1]),IAPWS95_entropy(rhoV[2],Ta[2]));
}

}

