
public class IAPWS
{
  /*
  Miscelaneous IAPWS standards. This module include:

    * :func:`_Ice`: Ice Ih state equation
    * :func:`_Liquid`: Properties of liquid water at 0.1 MPa
    * :func:`_Supercooled`: Thermodynamic properties of supercooled water
    * :func:`_Sublimation_Pressure`: Sublimation pressure correlation
    * :func:`_Melting_Pressure`: Melting pressure correlation
    * :func:`_Viscosity`: Viscosity correlation
    * :func:`_ThCond`: Themal conductivity correlation
    * :func:`_Tension`: Surface tension correlation
    * :func:`_Dielectric`: Dielectric constant correlation
    * :func:`_Refractive`: Refractive index correlation
    * :func:`_Kw`: Ionization constant correlation for ordinary water
    * :func:`_Conductivity`: Electrolytic conductivity correlation
    * :func:`_D2O_Viscosity`: Viscosity correlation for heavy water
    * :func:`_D2O_ThCond`: Thermal conductivity correlation for heavy water
    * :func:`_D2O_Tension`: Surface tension correlation for heavy water
    * :func:`_D2O_Sublimation_Pressure`: Sublimation Pressure correlation for
    heavy water
    * :func:`_D2O_Melting_Pressure`: Melting Pressure correlation for heavy
    water
    * :func:`_Henry`: Henry constant for liquid-gas equilibrium
    * :func:`_Kvalue`: Vapor-liquid distribution constant
  */

  // Constants
  private double M = 18.015268;   // g/mol
  private double R = 0.461526;    // kJ/kg·K

  // Table 1 from Release on the Values of Temperature, Pressure and Density of
  // Ordinary and Heavy Water Substances at their Respective Critical Points
  private double Tc = 647.096;    // K
  private double Pc = 22.064;     // MPa
  private double rhoc = 322.;     // kg/m³
  private double Tc_D2O = 643.847;  // K
  private double Pc_D2O = 21.671;   // MPa
  private double rhoc_D2O = 355.9999698294;  // kg/m³

  private double Tt = 273.16;     // K
  private double Pt = 611.657e-6;   // MPa
  private double Tb = 373.1243;   // K
  private double f_acent = 0.3443;

  // IAPWS, Guideline on the Use of Fundamental Physical Constants and Basic
  // Constants of Water, http://www.iapws.org/relguide/fundam.pdf
  private double Dipole = 1.85498;  // Debye

  // IAPWS-06 for Ice
  public double _Ice(double T, double P)
  {
    /*
    """Basic state equation for Ice Ih

    Parameters
    ----------
    T : float
      Temperature, [K]
    P : float
      Pressure, [MPa]

    Returns
    -------
    prop : dict
      Dict with calculated properties of ice. The available properties are:

        * rho: Density, [kg/m³]
        * h: Specific enthalpy, [kJ/kg]
        * u: Specific internal energy, [kJ/kg]
        * a: Specific Helmholtz energy, [kJ/kg]
        * g: Specific Gibbs energy, [kJ/kg]
        * s: Specific entropy, [kJ/kgK]
        * cp: Specific isobaric heat capacity, [kJ/kgK]
        * alfav: Cubic expansion coefficient, [1/K]
        * beta: Pressure coefficient, [MPa/K]
        * xkappa: Isothermal compressibility, [1/MPa]
        * ks: Isentropic compressibility, [1/MPa]
        * gt: [∂g/∂T]P
        * gtt: [∂²g/∂T²]P
        * gp: [∂g/∂P]T
        * gpp: [∂²g/∂P²]T
        * gtp: [∂²g/∂T∂P]

    Notes
    ------
    Raise :class:`NotImplementedError` if input isn't in limit:

      * T ≤ 273.16
      * P ≤ 208.566
      * State below the melting and sublimation lines

    Examples
    --------
    >>> st1 = _Ice(100, 100)
    >>> st1["rho"], st1["h"], st1["s"]
    941.678203297 -483.491635676 -2.61195122589

    >>> st2 = _Ice(273.152519,0.101325)
    >>> st2["a"], st2["u"], st2["cp"]
    -0.00918701567 -333.465403393 2.09671391024

    >>> st3 = _Ice(273.16,611.657e-6)
    >>> st3["alfav"], st3["beta"], st3["xkappa"], st3["ks"]
    0.000159863102566 1.35714764659 1.17793449348e-04 1.14161597779e-04

    References
    ----------
    IAPWS, Revised Release on the Equation of State 2006 for H2O Ice Ih
    September 2009, http://iapws.org/relguide/Ice-2009.html
    """
    */

    // Check input in range of validity
    if(T > 273.16)
    {
      // No Ice Ih stable
      warnings.warn("Metastable ice");
    }
    else if(P > 208.566)
    {
      // Ice Ih limit upper pressure
      raise NotImplementedError("Incoming out of bound");
    }
    else if(P < Pt)
    {
      Psub = _Sublimation_Pressure(T);
      if(Psub > P)
      {
        // Zone Gas
        warnings.warn("Metastable ice in vapor region");
      }
    }
    else if(251.165 < T)
    {
      Pmel = _Melting_Pressure(T);
      if(Pmel < P)
      {
        // Zone Liquid
        warnings.warn("Metastable ice in liquid region");
      }
    }

    double Tr = T/Tt;
    double Pr = P/Pt;
    double P0 = 101325e-6/Pt;
    double s0 = -0.332733756492168e4*1e-3;  // Express in kJ/kgK

    gok = [-0.632020233335886e6, 0.655022213658955, -0.189369929326131e-7, 0.339746123271053e-14, -0.556464869058991e-21];
    r2k = [complex(-0.725974574329220e2, -0.781008427112870e2)*1e-3, complex(-0.557107698030123e-4, 0.464578634580806e-4)*1e-3, complex(0.234801409215913e-10, -0.285651142904972e-10)*1e-3];
    t1 = complex(0.368017112855051e-1, 0.510878114959572e-1);
    t2 = complex(0.337315741065416, 0.335449415919309);
    r1 = complex(0.447050716285388e2, 0.656876847463481e2)*1e-3;

    double go, gop, gopp = 0;
    for(int k = 0; k < 5; k++)
    {
      go += gok[k]*1e-3*(Pr-P0)**k;
    }
    for(int k = 1; k < 5; k++)
    {
      gop += gok[k]*1e-3*k/Pt*(Pr-P0)**(k-1);
    }
    for(int k = 2; k < 5; k++)
    {
      gopp += gok[k]*1e-3*k*(k-1)/Pt**2*(Pr-P0)**(k-2);
    }
    r2 = r2p = 0;
    for(int k = 0; k < 3; k++)
    {
      r2 += r2k[k]*(Pr-P0)**k;
    }
    for(int k = 1; k < 3; k++)
    {
      r2p += r2k[k]*k/Pt*(Pr-P0)**(k-1);
    }
    r2pp = r2k[2]*2/Pt**2;

    double c = r1*((t1-Tr)*log_c(t1-Tr)+(t1+Tr)*log_c(t1+Tr)-2*t1*log_c(t1)-Tr**2/t1)+r2*((t2-Tr)*log_c(t2-Tr)+(t2+Tr)*log_c(t2+Tr)-2*t2*log_c(t2)-Tr**2/t2);
    double ct = r1*(-log_c(t1-Tr)+log_c(t1+Tr)-2*Tr/t1)+r2*(-log_c(t2-Tr)+log_c(t2+Tr)-2*Tr/t2);
    double ctt = r1*(1/(t1-Tr)+1/(t1+Tr)-2/t1) + r2*(1/(t2-Tr)+1/(t2+Tr)-2/t2);
    double cp = r2p*((t2-Tr)*log_c(t2-Tr)+(t2+Tr)*log_c(t2+Tr)-2*t2*log_c(t2)-Tr**2/t2);
    double ctp = r2p*(-log_c(t2-Tr)+log_c(t2+Tr)-2*Tr/t2);
    double cpp = r2pp*((t2-Tr)*log_c(t2-Tr)+(t2+Tr)*log_c(t2+Tr)-2*t2*log_c(t2)-Tr**2/t2);

    double g = go-s0*Tt*Tr+Tt*c.real;
    double gt = -s0+ct.real;
    double gp = gop+Tt*cp.real;
    double gtt = ctt.real/Tt;
    double gtp = ctp.real;
    double gpp = gopp+Tt*cpp.real;

    propiedades = {}
    propiedades["gt"] = gt;
    propiedades["gp"] = gp;
    propiedades["gtt"] = gtt;
    propiedades["gpp"] = gpp;
    propiedades["gtp"] = gtp;
    propiedades["T"] = T;
    propiedades["P"] = P;
    propiedades["v"] = gp/1000;
    propiedades["rho"] = 1000./gp;
    propiedades["h"] = g-T*gt;
    propiedades["s"] = -gt;
    propiedades["cp"] = -T*gtt;
    propiedades["u"] = g-T*gt-P*gp;
    propiedades["g"] = g;
    propiedades["a"] = g-P*gp;
    propiedades["alfav"] = gtp/gp;
    propiedades["beta"] = -gtp/gpp;
    propiedades["xkappa"] = -gpp/gp;
    propiedades["ks"] = (gtp**2-gtt*gpp)/gp/gtt;
    return propiedades;
  }

  // IAPWS-08 for Liquid water at 0.1 MPa
  public double _Liquid(double T, double P=0.1)
  {
    /*
    """Supplementary release on properties of liquid water at 0.1 MPa

    Parameters
    ----------
    T : float
      Temperature, [K]
    P : float
      Pressure, [MPa]
      Although this relation is for P=0.1MPa, can be extrapoled at pressure
      0.3 MPa

    Returns
    -------
    prop : dict
      Dict with calculated properties of water. The available properties are:

        * h: Specific enthalpy, [kJ/kg]
        * u: Specific internal energy, [kJ/kg]
        * a: Specific Helmholtz energy, [kJ/kg]
        * g: Specific Gibbs energy, [kJ/kg]
        * s: Specific entropy, [kJ/kgK]
        * cp: Specific isobaric heat capacity, [kJ/kgK]
        * cv: Specific isochoric heat capacity, [kJ/kgK]
        * w: Speed of sound, [m/s²]
        * rho: Density, [kg/m³]
        * v: Specific volume, [m³/kg]
        * vt: [∂v/∂T]P, [m³/kgK]
        * vtt: [∂²v/∂T²]P, [m³/kgK²]
        * vp: [∂v/∂P]T, [m³/kg/MPa]
        * vtp: [∂²v/∂T∂P], [m³/kg/MPa]
        * alfav: Cubic expansion coefficient, [1/K]
        * xkappa : Isothermal compressibility, [1/MPa]
        * ks: Isentropic compressibility, [1/MPa]
        * mu: Viscosity, [mPas]
        * k: Thermal conductivity, [W/mK]
        * epsilon: Dielectric constant, [-]

    Notes
    ------
    Raise :class:`NotImplementedError` if input isn't in limit:

      * 253.15 ≤ T ≤ 383.15
      * 0.1 ≤ P ≤ 0.3

    Examples
    --------
    >>> st1 = _Liquid(260)
    >>> st1["rho"], st1["h"], st1["s"]
    997.0683602710492 -55.86223174460868 -0.20998554842619535

    References
    ----------
    IAPWS, Revised Supplementary Release on Properties of Liquid Water at 0.1
    MPa, http://www.iapws.org/relguide/LiquidWater.html
    """
    */

    // Check input in range of validity
    if(T <= 253.15 || T >= 383.15 || P < 0.1 || P > 0.3)
    {
      raise NotImplementedError("Incoming out of bound");
    }
    else if(P != 0.1)
    {
      // Raise a warning if the P value is extrapolated
      warnings.warn("Using extrapolated values");
    }

    R = 0.46151805;   // kJ/kgK
    Po = 0.1;
    Tr = 10;
    tau = T/Tr;
    alfa = Tr/(593-T);
    beta = Tr/(T-232);

    a = [None, -1.661470539e5, 2.708781640e6, -1.557191544e8, None, 1.93763157e-2, 6.74458446e3, -2.22521604e5, 1.00231247e8, -1.63552118e9, 8.32299658e9, -7.5245878e-6, -1.3767418e-2, 1.0627293e1, -2.0457795e2, 1.2037414e3];
    b = [None, -8.237426256e-1, 1.908956353, -2.017597384, 8.546361348e-1, 5.78545292e-3, -1.53195665E-2, 3.11337859e-2, -4.23546241e-2, 3.38713507e-2, -1.19946761e-2, -3.1091470e-6, 2.8964919e-5, -1.3112763e-4, 3.0410453e-4, -3.9034594e-4, 2.3403117e-4, -4.8510101e-5];
    c = [None, -2.452093414e2, 3.869269598e1, -8.983025854];
    n = [None, 4, 5, 7, None, None, 4, 5, 7, 8, 9, 1, 3, 5, 6, 7];
    m = [None, 2, 3, 4, 5, 1, 2, 3, 4, 5, 6, 1, 3, 4, 5, 6, 7, 9];

    sum = 0;
    for(int i = 1; i < 4; i++)
    {
      sum += a[i]*alfa**n[i];
    }
    suma1 = sum;
    sum = 0;
    for(int i = 1; i < 5; i++)
    {
      sum += b[i]*beta**m[i];
    }
    suma2 = sum;
    go = R*Tr*(c[1]+c[2]*tau+c[3]*tau*log(tau)+suma1+suma2);

    sum = 0;
    for(int i = 6; i < 11; i++)
    {
      sum += a[i]*alfa**n[i];
    }
    suma1 = sum;
    sum = 0;;
    for(int i = 5; i < 11; i++)
    {
      sum += b[i]*beta**m[i];
    }
    suma2 = sum;
    vo = R*Tr/Po/1000*(a[5]+suma1+suma2);

    sum = 0;
    for(int i = 11; i < 16; i++)
    {
      sum += a[i]*alfa**n[i];
    }
    suma1 = sum;
    sum = 0;;
    for(int i = 11; i < 18; i++)
    {
      sum += b[i]*beta**m[i];
    }
    suma2 = sum;
    vpo = R*Tr/Po**2/1000*(suma1+suma2)

    sum = 0;
    for(int i = 1; i < 4; i++)
    {
      sum += n[i]*a[i]*alfa**(n[i]+1);
    }
    suma1 = sum;
    sum = 0;
    for(int i = 1; i < 5; i++)
    {
      sum += m[i]*b[i]*beta**(m[i]+1);
    }
    suma2 = sum;
    so = -R*(c[2]+c[3]*(1+log(tau))+suma1-suma2);

    sum = 0;
    for(int i = 1; i < 4; i++)
    {
      sum += n[i]*(n[i]+1)*a[i]*alfa**(n[i]+2);
    }
    suma1 = sum;
    sum = 0;
    for(int i = 1; i < 5; i++)
    {
      sum += m[i]*(m[i]+1)*b[i]*beta**(m[i]+2);
    }
    suma2 = sum;
    cpo = -R*(c[3]+tau*suma1+tau*suma2);

    sum = 0;
    for(int i = 6; i < 11; i++);
    {
      sum += n[i]*a[i]*alfa**(n[i]+1);
    }
    suma1 = sum;
    sum = 0;
    for(int i = 5; i < 11; i++)
    {
      sum += m[i]*b[i]*beta**(m[i]+1);
    }
    suma2 = sum;
    vto = R/Po/1000*(suma1-suma2);

    // This properties are only neccessary for computing thermodynamic
    // properties at pressures different from 0.1 MPa
    sum = 0;
    for(int i = 6; i < 11; i++);
    {
      sum += n[i]*(n[i]+1)*a[i]*alfa**(n[i]+2);
    }
    suma1 = sum;
    sum = 0;
    for(int i = 5; i < 11; i++)
    {
      sum += m[i]*(m[i]+1)*b[i]*beta**(m[i]+2);
    }
    suma2 = sum;
    vtto = R/Tr/Po/1000*(suma1+suma2);

    sum = 0;
    for(int i = 11; i < 16; i++)
    {
      sum += n[i]*a[i]*alfa**(n[i]+1);
    }
    suma1 = sum;
    sum = 0;
    for(int i = 11; i < 18; i++)
    {
      sum += m[i]*b[i]*beta**(m[i]+1);
    }
    suma2 = sum;
    vpto = R/Po**2/1000*(suma1-suma2);

    if(P != 0.1)
    {
      go += vo*(P-0.1);
      so -= vto*(P-0.1);
      cpo -= T*vtto*(P-0.1);
      vo -= vpo*(P-0.1);
      vto += vpto*(P-0.1);
      vppo = 3.24e-10*R*Tr/0.1**3;
      vpo += vppo*(P-0.1);
    }

    h = go+T*so;
    u = h-P*vo;
    a = go-P*vo;
    cv = cpo+T*vto**2/vpo;
    xkappa = -vpo/vo;
    alfa = vto/vo;
    ks = -(T*vto**2/cpo+vpo)/vo;
    w = (-vo**2*1e9/(vpo*1e3+T*vto**2*1e6/cpo))**0.5;

    propiedades = {};
    propiedades["g"] = go;
    propiedades["T"] = T;
    propiedades["P"] = P;
    propiedades["v"] = vo;
    propiedades["vt"] = vto;
    propiedades["vp"] = vpo;
    propiedades["vpt"] = vpto;
    propiedades["vtt"] = vtto;
    propiedades["rho"] = 1/vo;
    propiedades["h"] = h;
    propiedades["s"] = so;
    propiedades["cp"] = cpo;
    propiedades["cv"] = cv;
    propiedades["u"] = u;
    propiedades["a"] = a;
    propiedades["xkappa"] = xkappa;
    propiedades["alfav"] = vto/vo;
    propiedades["ks"] = ks;
    propiedades["w"] = w;

    // Viscosity correlation, Eq 7
    a = [None, 280.68, 511.45, 61.131, 0.45903];
    b = [None, -1.9, -7.7, -19.6, -40];
    T_ = T/300;
    sum = 0;;
    for(int i = 1; i < 5; i++)
    {
      sum += a[i]*T_**b[i];
    }
    mu = sum/1e6;
    propiedades["mu"] = mu;

    // Thermal conductivity correlation, Eq 8
    c = [None, 1.6630, -1.7781, 1.1567, -0.432115];
    d = [None, -1.15, -3.4, -6.0, -7.6];
    sum = 0;;
    for(int i = 1; i < 5; i++)
    {
      sum += c[i]*T_**d[i];
    }
    k = sum;
    propiedades["k"] = k;

    // Dielectric constant correlation, Eq 9
    e = [None, -43.7527, 299.504, -399.364, 221.327];
    f = [None, -0.05, -1.47, -2.11, -2.31];
    sum = 0;
    for(int i = 1; i < 5; i++)
    {
      sum += e[i]*T_**f[i];
    }
    epsilon = sum;
    propiedades["epsilon"] = epsilon;

    return propiedades;
  }
}

//!/usr/bin/python
// -*- coding: utf-8 -*-

from __future__ import division

from cmath import log as log_c
from math import log, exp, tan, atan, acos, sin, pi, log10
import warnings

from scipy.optimize import minimize

// IAPWS-15 for supercooled liquid water
void _Supercooled(T, P)
{
  /*
  """Guideline on thermodynamic properties of supercooled water

  Parameters
  ----------
  T : float
    Temperature, [K]
  P : float
    Pressure, [MPa]

  Returns
  -------
  prop : dict
    Dict with calculated properties of water. The available properties are:

      * L: Ordering field, [-]
      * x: Mole fraction of low-density structure, [-]
      * rho: Density, [kg/m³]
      * s: Specific entropy, [kJ/kgK]
      * h: Specific enthalpy, [kJ/kg]
      * u: Specific internal energy, [kJ/kg]
      * a: Specific Helmholtz energy, [kJ/kg]
      * g: Specific Gibbs energy, [kJ/kg]
      * alfap: Thermal expansion coefficient, [1/K]
      * xkappa : Isothermal compressibility, [1/MPa]
      * cp: Specific isobaric heat capacity, [kJ/kgK]
      * cv: Specific isochoric heat capacity, [kJ/kgK]
      * w: Speed of sound, [m/s²]

  Notes
  ------
  Raise :class:`NotImplementedError` if input isn't in limit:

    * Tm ≤ T ≤ 300
    * 0 < P ≤ 1000

  The minimum temperature in range of validity is the melting temperature, it
  depend of pressure

  Examples
  --------
  >>> liq = _Supercooled(235.15, 0.101325)
  >>> liq["rho"], liq["cp"], liq["w"]
  968.09999 5.997563 1134.5855

  References
  ----------
  IAPWS, Guideline on Thermodynamic Properties of Supercooled Water,
  http://iapws.org/relguide/Supercooled.html
  """
  */

  // Check input in range of validity
  if(P < 198.9)
  {
    Tita = T/235.15
    Ph = 0.1+228.27*(1-Tita**6.243)+15.724*(1-Tita**79.81)
    if(P < Ph || T > 300)
    {
      raise NotImplementedError("Incoming out of bound")
    }
  }
  else
  {
    Th = 172.82+0.03718*P+3.403e-5*P**2-1.573e-8*P**3
    if(T < Th || T > 300 || P > 1000)
    {
      raise NotImplementedError("Incoming out of bound")
    }
  }
  // Parameters, Table 1
  Tll = 228.2
  rho0 = 1081.6482
  R = 0.461523087
  pi0 = 300e3/rho0/R/Tll
  omega0 = 0.5212269
  L0 = 0.76317954
  k0 = 0.072158686
  k1 = -0.31569232
  k2 = 5.2992608

  // Reducing parameters, Eq 2
  tau = T/Tll-1
  p = P*1000/rho0/R/Tll
  tau_ = tau+1
  p_ = p+pi0

  // Eq 3
  ci = [-8.1570681381655, 1.2875032, 7.0901673598012, -3.2779161e-2, 7.3703949e-1, -2.1628622e-1, -5.1782479, 4.2293517e-4, 2.3592109e-2, 4.3773754, -2.9967770e-3, -9.6558018e-1, 3.7595286, 1.2632441, 2.8542697e-1, -8.5994947e-1, -3.2916153e-1, 9.0019616e-2, 8.1149726e-2, -3.2788213]
  ai = [0, 0, 1, -0.2555, 1.5762, 1.6400, 3.6385, -0.3828, 1.6219, 4.3287, 3.4763, 5.1556, -0.3593, 5.0361, 2.9786, 6.2373, 4.0460, 5.3558, 9.0157, 1.2194]
  bi = [0, 1, 0, 2.1051, 1.1422, 0.9510, 0, 3.6402, 2.0760, -0.0016, 2.2769, 0.0008, 0.3706, -0.3975, 2.9730, -0.3180, 2.9805, 2.9265, 0.4456, 0.1298]
  di = [0, 0, 0, -0.0016, 0.6894, 0.0130, 0.0002, 0.0435, 0.0500, 0.0004, 0.0528, 0.0147, 0.8584, 0.9924, 1.0041, 1.0961, 1.0228, 1.0303, 1.6180, 0.5213]
  phir = phirt = phirp = phirtt = phirtp = phirpp = 0
  for(int i = 0; i < ci.length; i++)
  {
    double c = ci[i];
    double a = ai[i];
    double b = bi[i];
    double d = di[i];
    phir += c*tau_**a*p_**b*exp(-d*p_)
    phirt += c*a*tau_**(a-1)*p_**b*exp(-d*p_)
    phirp += c*tau_**a*p_**(b-1)*(b-d*p_)*exp(-d*p_)
    phirtt += c*a*(a-1)*tau_**(a-2)*p_**b*exp(-d*p_)
    phirtp += c*a*tau_**(a-1)*p_**(b-1)*(b-d*p_)*exp(-d*p_)
    phirpp += c*tau_**a*p_**(b-2)*((d*p_-b)**2-b)*exp(-d*p_)
  }

  // Eq 5
  K1 = ((1+k0*k2+k1*(p-k2*tau))**2-4*k0*k1*k2*(p-k2*tau))**0.5
  K2 = (1+k2**2)**0.5

  // Eq 6
  omega = 2+omega0*p

  // Eq 4
  L = L0*K2/2/k1/k2*(1+k0*k2+k1*(p+k2*tau)-K1)

  // Define interval of solution, Table 4
  if(omega < 10/9*(log(19)-L))
  {
    xmin = 0.049
    xmax = 0.5
  }
  else if(10/9*(log(19)-L) <= omega < 50/49*(log(99)-L))
  {
    xmin = 0.0099
    xmax = 0.051
  }
  else
  {
    xmin = 0.99*exp(-50/49*L-omega)
    xmax = min(1.1*exp(-L-omega), 0.0101)

  }
  void f(x)
  {
    return abs(L+log(x/(1-x))+omega*(1-2*x))
  }
  x = minimize(f, ((xmin+xmax)/2,), bounds=((xmin, xmax),))["x"][0]

  // Eq 12
  fi = 2*x-1
  Xi = 1/(2/(1-fi**2)-omega)

  // Derivatives, Table 3
  Lt = L0*K2/2*(1+(1-k0*k2+k1*(p-k2*tau))/K1)
  Lp = L0*K2*(K1+k0*k2-k1*p+k1*k2*tau-1)/2/k2/K1
  Ltt = -2*L0*K2*k0*k1*k2**2/K1**3
  Ltp = 2*L0*K2*k0*k1*k2/K1**3
  Lpp = -2*L0*K2*k0*k1/K1**3

  prop = {}
  prop["L"] = L
  prop["x"] = x

  // Eq 13
  prop["rho"] = rho0/((tau+1)/2*(omega0/2*(1-fi**2)+Lp*(fi+1))+phirp)

  // Eq 1
  prop["g"] = phir+(tau+1)*(x*L+x*log(x)+(1-x)*log(1-x)+omega*x*(1-x))

  // Eq 14
  prop["s"] = -R*((tau+1)/2*Lt*(fi+1) +
          (x*L+x*log(x)+(1-x)*log(1-x)+omega*x*(1-x))+phirt)

  // Basic derived state properties
  prop["h"] = prop["g"]+T*prop["s"]
  prop["u"] = prop["h"]+P/prop["rho"]
  prop["a"] = prop["u"]-T*prop["s"]

  // Eq 15
  prop["xkappa"] = prop["rho"]/rho0**2/R*1000/Tll*(
    (tau+1)/2*(Xi*(Lp-omega0*fi)**2-(fi+1)*Lpp)-phirpp)
  prop["alfap"] = prop["rho"]/rho0/Tll*(
    Ltp/2*(tau+1)*(fi+1) + (omega0*(1-fi**2)/2+Lp*(fi+1))/2 -
    (tau+1)*Lt/2*Xi*(Lp-omega0*fi) + phirtp)
  prop["cp"] = -R*(tau+1)*(Lt*(fi+1)+(tau+1)/2*(Ltt*(fi+1)-Lt**2*Xi)+phirtt)

  // Eq 16
  prop["cv"] = prop["cp"]-T*prop["alfap"]**2/prop["rho"]/prop["xkappa"]*1e3

  // Eq 17
  prop["w"] = (prop["rho"]*prop["xkappa"]*1e-6*prop["cv"]/prop["cp"])**-0.5
  return prop
}

void _Sublimation_Pressure(T)
{
  /*
  """Sublimation Pressure correlation

  Parameters
  ----------
  T : float
    Temperature, [K]

  Returns
  -------
  P : float
    Pressure at sublimation line, [MPa]

  Notes
  ------
  Raise :class:`NotImplementedError` if input isn't in limit:

    * 50 ≤ T ≤ 273.16

  Examples
  --------
  >>> _Sublimation_Pressure(230)
  8.947352740189152e-06

  References
  ----------
  IAPWS, Revised Release on the Pressure along the Melting and Sublimation
  Curves of Ordinary Water Substance, http://iapws.org/relguide/MeltSub.html.
  """
  */
  if(50 <= T <= 273.16)
  {
    Tita = T/Tt
    suma = 0
    a = [-0.212144006e2, 0.273203819e2, -0.61059813e1]
    expo = [0.333333333e-2, 1.20666667, 1.70333333]
    for(int i = 0; i < a.length; i++)
    {
      double ai = a[i];
      double expi = expo[i];
      suma += ai*Tita**expi
    }
    return exp(suma/Tita)*Pt
  }
  else
  {
    raise NotImplementedError("Incoming out of bound")
  }


}
void _Melting_Pressure(T, ice="Ih")
{
  /*
  """Melting Pressure correlation

  Parameters
  ----------
  T : float
    Temperature, [K]
  ice: string
    Type of ice: Ih, III, V, VI, VII.
    Below 273.15 is a mandatory input, the ice Ih is the default value.
    Above 273.15, the ice type is unnecesary.

  Returns
  -------
  P : float
    Pressure at sublimation line, [MPa]

  Notes
  ------
  Raise :class:`NotImplementedError` if input isn't in limit:

    * 251.165 ≤ T ≤ 715

  Examples
  --------
  >>> _Melting_Pressure(260)
  8.947352740189152e-06
  >>> _Melting_Pressure(254, "III")
  268.6846466336108

  References
  ----------
  IAPWS, Revised Release on the Pressure along the Melting and Sublimation
  Curves of Ordinary Water Substance, http://iapws.org/relguide/MeltSub.html.
  """
  */
  if(ice == "Ih" && 251.165 <= T <= 273.16)
  {
    // Ice Ih
    Tref = Tt
    Pref = Pt
    Tita = T/Tref
    a = [0.119539337e7, 0.808183159e5, 0.33382686e4]
    expo = [3., 0.2575e2, 0.10375e3]
    suma = 1
    for(int i = 0; i < a.length; i++)
    {
      ai = a[i];
      expi= expo[i];
      suma += ai*(1-Tita**expi)
    }
    P = suma*Pref
  }
  else if(ice == "III" && 251.165 < T <= 256.164)
  {
    // Ice III
    Tref = 251.165
    Pref = 208.566
    Tita = T/Tref
    P = Pref*(1-0.299948*(1-Tita**60.))
  }
  else if((ice == "V" && 256.164 < T <= 273.15) || 273.15 < T <= 273.31)
  {
    // Ice V
    Tref = 256.164
    Pref = 350.100
    Tita = T/Tref
    P = Pref*(1-1.18721*(1-Tita**8.))
  }
  else if(273.31 < T <= 355)
  {
    // Ice VI
    Tref = 273.31
    Pref = 632.400
    Tita = T/Tref
    P = Pref*(1-1.07476*(1-Tita**4.6))
  }
  else if(355. < T <= 715)
  {
    // Ice VII
    Tref = 355
    Pref = 2216.000
    Tita = T/Tref
    P = Pref*exp(1.73683*(1-1./Tita)-0.544606e-1*(1-Tita**5) +
           0.806106e-7*(1-Tita**22))
  }
  else
  {
    raise NotImplementedError("Incoming out of bound")
  }
  return P
}

// Transport properties
void _Viscosity(rho, T, fase=None, drho=None)
{
  /*
  """Equation for the Viscosity

  Parameters
  ----------
  rho : float
    Density, [kg/m³]
  T : float
    Temperature, [K]
  fase: dict, optional for calculate critical enhancement
    phase properties
  drho: float, optional for calculate critical enhancement
    [∂ρ/∂P]T at reference state,

  Returns
  -------
  μ : float
    Viscosity, [Pa·s]

  Examples
  --------
  >>> _Viscosity(998, 298.15)
  0.0008897351001498108
  >>> _Viscosity(600, 873.15)
  7.743019522728247e-05

  References
  ----------
  IAPWS, Release on the IAPWS Formulation 2008 for the Viscosity of Ordinary
  Water Substance, http://www.iapws.org/relguide/viscosity.html
  """
  */
  Tr = T/Tc
  Dr = rho/rhoc

  // Eq 11
  H = [1.67752, 2.20462, 0.6366564, -0.241605]
  sum = 0;
  for(int i = 0; i < H.length; i++)
  {
    Hi = H[i];
    sum += Hi/Tr**i;
  }
  mu0 = 100*Tr**0.5/sum

  // Eq 12
  I = [0, 1, 2, 3, 0, 1, 2, 3, 5, 0, 1, 2, 3, 4, 0, 1, 0, 3, 4, 3, 5]
  J = [0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 4, 4, 5, 6, 6]
  Hij = [0.520094, 0.850895e-1, -0.108374e1, -0.289555, 0.222531, 0.999115, 0.188797e1, 0.126613e1, 0.120573, -0.281378, -0.906851, -0.772479, -0.489837, -0.257040, 0.161913, 0.257399, -0.325372e-1, 0.698452e-1, 0.872102e-2, -0.435673e-2, -0.593264e-3]
  sum = 0;
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    h = Hij[iter];
    sum += (1/Tr-1)**i*h*(Dr-1)**j;
  }
  mu1 = exp(Dr*sum)

  // Critical enhancement
  if(fase && drho)
  {
    qc = 1/1.9
    qd = 1/1.1

    // Eq 21
    DeltaX = Pc*Dr**2*(fase.drhodP_T/rho-drho/rho*1.5/Tr)
    if(DeltaX < 0)
    {
      DeltaX = 0
    }

    // Eq 20
    X = 0.13*(DeltaX/0.06)**(0.63/1.239)

    if(X <= 0.3817016416)
    {
      // Eq 15
      Y = qc/5*X*(qd*X)**5*(1-qc*X+(qc*X)**2-765./504*(qd*X)**2)
    }
    else
    {
      Fid = acos((1+qd**2*X**2)**-0.5)              // Eq 17
      w = abs((qc*X-1)/(qc*X+1))**0.5*tan(Fid/2)          // Eq 19

      // Eq 18
      if(qc*X > 1)
      {
        Lw = log((1+w)/(1-w))
      }
      else
      {
        Lw = 2*atan(abs(w))
      }

      // Eq 16
      Y = sin(3*Fid)/12-sin(2*Fid)/4/qc/X+(1-5/4*(qc*X)**2)/(
        qc*X)**2*sin(Fid)-((1-3/2*(qc*X)**2)*Fid-abs((
          qc*X)**2-1)**1.5*Lw)/(qc*X)**3

    }
    // Eq 14
    mu2 = exp(0.068*Y)
  }
  else
  {
    mu2 = 1

  }
  // Eq 10
  mu = mu0*mu1*mu2
  return mu*1e-6
}

void _ThCond(rho, T, fase=None, drho=None)
{
  /*
  """Equation for the thermal conductivity

  Parameters
  ----------
  rho : float
    Density, [kg/m³]
  T : float
    Temperature, [K]
  fase: dict, optional for calculate critical enhancement
    phase properties
  drho: float, optional for calculate critical enhancement
    [∂ρ/∂P]T at reference state,

  Returns
  -------
  k : float
    Thermal conductivity, [W/mK]

  Examples
  --------
  >>> _ThCond(998, 298.15)
  0.6077128675880629
  >>> _ThCond(0, 873.15)
  0.07910346589648833

  References
  ----------
  IAPWS, Release on the IAPWS Formulation 2011 for the Thermal Conductivity
  of Ordinary Water Substance, http://www.iapws.org/relguide/ThCond.html
  """
  */
  d = rho/rhoc
  Tr = T/Tc

  // Eq 16
  no = [2.443221e-3, 1.323095e-2, 6.770357e-3, -3.454586e-3, 4.096266e-4]
  sum = 0;
  for(int i = 0; i < no.length; i++)
  {
    n = no[i];
    sum += n/Tr**i;
  }
  k0 = Tr**0.5/sum

  // Eq 17
  I = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4]
  J = [0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 0, 1, 2, 3, 4, 5]
  nij = [1.60397357, -0.646013523, 0.111443906, 0.102997357, -0.0504123634, 0.00609859258, 2.33771842, -2.78843778, 1.53616167, -0.463045512, 0.0832827019, -0.00719201245, 2.19650529, -4.54580785, 3.55777244, -1.40944978, 0.275418278, -0.0205938816, -1.21051378, 1.60812989, -0.621178141, 0.0716373224, -2.7203370, 4.57586331, -3.18369245, 1.1168348, -0.19268305, 0.012913842]
  sum = 0;
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    n = n[iter];
    sum += (1/Tr-1)**i*n*(d-1)**j;
  }
  k1 = exp(d*sum)

  // Critical enhancement
  if(fase)
  {
    R = 0.46151805

    if(!drho)
    {
      // Industrial formulation
      // Eq 25
      if(d <= 0.310559006)
      {
        ai = [6.53786807199516, -5.61149954923348, 3.39624167361325, -2.27492629730878, 10.2631854662709, 1.97815050331519]
      }
      else if(d <= 0.776397516)
      {
        ai = [6.52717759281799, -6.30816983387575, 8.08379285492595, -9.82240510197603, 12.1358413791395, -5.54349664571295]
      }
      else if(d <= 1.242236025)
      {
        ai = [5.35500529896124, -3.96415689925446, 8.91990208918795, -12.0338729505790, 9.19494865194302, -2.16866274479712]
      }
      else if(d <= 1.863354037)
      {
        ai = [1.55225959906681, 0.464621290821181, 8.93237374861479, -11.0321960061126, 6.16780999933360, -0.965458722086812]
      }
      else
      {
        ai = [1.11999926419994, 0.595748562571649, 9.88952565078920, -10.3255051147040, 4.66861294457414, -0.503243546373828]
      }
      sum = 0;
      for(int i = 0; i < ai.length; i++)
      {
        a = ai[i];
        sum += a*d**i;
      }
      drho = 1/sum*rhoc/Pc
    }

    DeltaX = d*(Pc/rhoc*fase.drhodP_T-Pc/rhoc*drho*1.5/Tr)
    if(DeltaX < 0)
    {
      DeltaX = 0
    }
    X = 0.13*(DeltaX/0.06)**(0.63/1.239)              // Eq 22
    y = X/0.4                             // Eq 20

    // Eq 19
    if(y < 1.2e-7)
    {
      Z = 0
    }
    else
    {
      Z = 2/pi/y*(((1-1/fase.cp_cv)*atan(y)+y/fase.cp_cv)-(
        1-exp(-1/(1/y+y**2/3/d**2))))
    }

    // Eq 18
    k2 = 177.8514*d*fase.cp/R*Tr/fase.mu*1e-6*Z
  }
  else
  {
    // No critical enhancement
    k2 = 0
  }
  // Eq 10
  k = k0*k1+k2
  return 1e-3*k
}

void _Tension(T)
{
  /*
  """Equation for the surface tension

  Parameters
  ----------
  T : float
    Temperature, [K]

  Returns
  -------
  σ : float
    Surface tension, [N/m]

  Notes
  ------
  Raise :class:`NotImplementedError` if input isn't in limit:

    * 248.15 ≤ T ≤ 647
    * Estrapolate to -25ºC in supercooled liquid metastable state

  Examples
  --------
  >>> _Tension(300)
  0.0716859625
  >>> _Tension(450)
  0.0428914992

  References
  ----------
  IAPWS, Revised Release on Surface Tension of Ordinary Water Substance
  June 2014, http://www.iapws.org/relguide/Surf-H2O.html
  """
  */
  if(248.15 <= T <= Tc)
  {
    Tr = T/Tc
    return 1e-3*(235.8*(1-Tr)**1.256*(1-0.625*(1-Tr)))
  }
  else
  {
    raise NotImplementedError("Incoming out of bound")
  }

}
void _Dielectric(rho, T)
{
  /*
  """Equation for the Dielectric constant

  Parameters
  ----------
  rho : float
    Density, [kg/m³]
  T : float
    Temperature, [K]

  Returns
  -------
  epsilon : float
    Dielectric constant, [-]

  Notes
  ------
  Raise :class:`NotImplementedError` if input isn't in limit:

    * 238 ≤ T ≤ 1200

  Examples
  --------
  >>> _Dielectric(999.242866, 298.15)
  78.5907250
  >>> _Dielectric(26.0569558, 873.15)
  1.12620970

  References
  ----------
  IAPWS, Release on the Static Dielectric Constant of Ordinary Water
  Substance for Temperatures from 238 K to 873 K and Pressures up to 1000
  MPa, http://www.iapws.org/relguide/Dielec.html
  """
  */
  // Check input parameters
  if(T < 238 || T > 1200)
  {
    raise NotImplementedError("Incoming out of bound")
  }

  k = 1.380658e-23
  Na = 6.0221367e23
  alfa = 1.636e-40
  epsilon0 = 8.854187817e-12
  mu = 6.138e-30

  d = rho/rhoc
  Tr = Tc/T
  I = [1, 1, 1, 2, 3, 3, 4, 5, 6, 7, 10, None]
  J = [0.25, 1, 2.5, 1.5, 1.5, 2.5, 2, 2, 5, 0.5, 10, None]
  n = [0.978224486826, -0.957771379375, 0.237511794148, 0.714692244396, -0.298217036956, -0.108863472196, .949327488264e-1, -.980469816509e-2, .165167634970e-4, .937359795772e-4, -.12317921872e-9, .196096504426e-2]

  g = 1+n[11]*d/(Tc/228/Tr-1)**1.2
  for(int i = 0; i < 11; i++)
  {
    g += n[i]*d**I[i]*Tr**J[i]
  }
  A = Na*mu**2*rho*g/M*1000/epsilon0/k/T
  B = Na*alfa*rho/3/M*1000/epsilon0
  e = (1+A+5*B+(9+2*A+18*B+A**2+10*A*B+9*B**2)**0.5)/4/(1-B)
  return e
}

void _Refractive(rho, T, l=0.5893)
{
  /*
  """Equation for the refractive index

  Parameters
  ----------
  rho : float
    Density, [kg/m³]
  T : float
    Temperature, [K]
  l : float, optional
    Light Wavelength, [μm]

  Returns
  -------
  n : float
    Refractive index, [-]

  Notes
  ------
  Raise :class:`NotImplementedError` if input isn't in limit:

    * 0 ≤ ρ ≤ 1060
    * 261.15 ≤ T ≤ 773.15
    * 0.2 ≤ λ ≤ 1.1

  Examples
  --------
  >>> _Refractive(997.047435, 298.15, 0.2265)
  1.39277824
  >>> _Refractive(30.4758534, 773.15, 0.5893)
  1.00949307

  References
  ----------
  IAPWS, Release on the Refractive Index of Ordinary Water Substance as a
  Function of Wavelength, Temperature and Pressure,
  http://www.iapws.org/relguide/rindex.pdf
  """
  */
  // Check input parameters
  if(rho < 0 || rho > 1060 || T < 261.15 || T > 773.15 || l < 0.2 || l > 1.1)
  {
    raise NotImplementedError("Incoming out of bound")
  }

  Lir = 5.432937
  Luv = 0.229202
  d = rho/1000.
  Tr = T/273.15
  L = l/0.589
  a = [0.244257733, 0.974634476e-2, -0.373234996e-2, 0.268678472e-3, 0.158920570e-2, 0.245934259e-2, 0.900704920, -0.166626219e-1]
  A = d*(a[0]+a[1]*d+a[2]*Tr+a[3]*L**2*Tr+a[4]/L**2+a[5]/(L**2-Luv**2)+a[6]/(
    L**2-Lir**2)+a[7]*d**2)
  return ((2*A+1)/(1-A))**0.5
}

void _Kw(rho, T)
{
  /*
  """Equation for the ionization constant of ordinary water

  Parameters
  ----------
  rho : float
    Density, [kg/m³]
  T : float
    Temperature, [K]

  Returns
  -------
  pKw : float
    Ionization constant in -log10(kw), [-]

  Notes
  ------
  Raise :class:`NotImplementedError` if input isn't in limit:

    * 0 ≤ ρ ≤ 1250
    * 273.15 ≤ T ≤ 1073.15

  Examples
  --------
  >>> _Kw(1000, 300)
  13.906565

  References
  ----------
  IAPWS, Release on the Ionization Constant of H2O,
  http://www.iapws.org/relguide/Ionization.pdf
  """
  */
  // Check input parameters
  if(rho < 0 || rho > 1250 || T < 273.15 || T > 1073.15)
  {
    raise NotImplementedError("Incoming out of bound")
  }

  // The internal method of calculation use rho in g/cm³
  d = rho/1000.

  // Water molecular weight different
  Mw = 18.015268

  gamma = [6.1415e-1, 4.825133e4, -6.770793e4, 1.01021e7]
  pKg = 0
  for(int i = 0; i < gamma.length; i++)
  {
    g = gamma[i];
    pKg += g/T**i;
  }

  Q = d*exp(-0.864671+8659.19/T-22786.2/T**2*d**(2./3))
  pKw = -12*(log10(1+Q)-Q/(Q+1)*d*(0.642044-56.8534/T-0.375754*d)) + \
  {
    pKg+2*log10(Mw/1000)
  }
  return pKw
}

void _Conductivity(rho, T)
{
  /*
  """Equation for the electrolytic conductivity of liquid and dense
  supercrítical water

  Parameters
  ----------
  rho : float
    Density, [kg/m³]
  T : float
    Temperature, [K]

  Returns
  -------
  K : float
    Electrolytic conductivity, [S/m]

  Notes
  ------
  Raise :class:`NotImplementedError` if input isn't in limit:

    * 600 ≤ ρ ≤ 1200
    * 273.15 ≤ T ≤ 1073.15

  Examples
  --------
  >>> _Conductivity(1000, 373.15)
  1.13

  References
  ----------
  IAPWS, Electrolytic Conductivity (Specific Conductance) of Liquid and Dense
  Supercritical Water from 0°C to 800°C and Pressures up to 1000 MPa,
  http://www.iapws.org/relguide/conduct.pdf
  """
  */
  // FIXME: Dont work
  rho_ = rho/1000
  kw = 10**-_Kw(rho, T)

  A = [1850., 1410., 2.16417e-6, 1.81609e-7, -1.75297e-9, 7.20708e-12]
  B = [16., 11.6, 3.26e-4, -2.3e-6, 1.1e-8]
  t = T-273.15

  sum = 0;
  for(int i = 0; i < 4; i++)
  {
    sum += A[i+2]*t**(i+1) ;
  }
  Loo = A[0]-1/(1/A[1]+sum)    // Eq 5
  sum = 0;
  for(int i = 0; i < 3; i++)
  {
    sum += B[i+2]*t**(i+1) ;
  }
  rho_h = B[0]-1/(1/B[1]+sum)  // Eq 6

  // Eq 4
  L_o = (rho_h-rho_)*Loo/rho_h

  // Eq 1
  k = 100*1e-3*L_o*kw**0.5*rho_
  return k
}

// Heavy water transport properties
void _D2O_Viscosity(rho, T)
{
  /*
  """Equation for the Viscosity of heavy water

  Parameters
  ----------
  rho : float
    Density, [kg/m³]
  T : float
    Temperature, [K]

  Returns
  -------
  μ : float
    Viscosity, [Pa·s]

  Examples
  --------
  >>> _D2O_Viscosity(998, 298.15)
  0.0008897351001498108
  >>> _D2O_Viscosity(600, 873.15)
  7.743019522728247e-05

  References
  ----------
  IAPWS, Revised Release on Viscosity and Thermal Conductivity of Heavy
  Water Substance, http://www.iapws.org/relguide/TransD2O-2007.pdf
  """
  */
  Tr = T/643.847
  rhor = rho/358.0

  no = [1.0, 0.940695, 0.578377, -0.202044]
  sum = 0;
  for(int i = 0; i < no.length; i++)
  {
    n = no[i];
    sum += n/Tr**i;
  }
  fi0 = Tr**0.5/sum

  Li = [0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 0, 1, 2, 5, 0, 1, 2, 3, 0, 1, 3, 5, 0, 1, 5, 3]
  Lj = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6]
  Lij = [0.4864192, -0.2448372, -0.8702035, 0.8716056, -1.051126, 0.3458395, 0.3509007, 1.315436, 1.297752, 1.353448, -0.2847572, -1.037026, -1.287846, -0.02148229, 0.07013759, 0.4660127, 0.2292075, -0.4857462, 0.01641220, -0.02884911, 0.1607171, -.009603846, -.01163815, -.008239587, 0.004559914, -0.003886659]

  sum = 0;
  for(int iter = 0; iter < Li.length; iter++)
  {
    i = Li[iter];
    j = Lj[iter];
    lij = Lij[iter];
    sum += lij*(1./Tr-1)**i*(rhor-1)**j;
  }
  fi1 = exp(rhor*sum)

  return 55.2651e-6*fi0*fi1
}

void _D2O_ThCond(rho, T)
{
  /*
  """Equation for the thermal conductivity of heavy water

  Parameters
  ----------
  rho : float
    Density, [kg/m³]
  T : float
    Temperature, [K]

  Returns
  -------
  k : float
    Thermal conductivity, [W/mK]

  Examples
  --------
  >>> _D2O_ThCond(998, 298.15)
  0.6077128675880629
  >>> _D2O_ThCond(0, 873.15)
  0.07910346589648833

  References
  ----------
  IAPWS, Revised Release on Viscosity and Thermal Conductivity of Heavy
  Water Substance, http://www.iapws.org/relguide/TransD2O-2007.pdf
  """
  */
  rhor = rho/358
  Tr = T/643.847
  tau = Tr/(abs(Tr-1.1)+1.1)

  no = [1.0, 37.3223, 22.5485, 13.0465, 0.0, -2.60735]
  sum = 0;
  for(int i = 0; i < no.length; i++)
  {
    Li = no[i];
    sum += Li*Tr**i;
  }
  Lo = sum

  nr = [483.656, -191.039, 73.0358, -7.57467]
  sum = 0;
  for(int i = 0; i < nr.length; i++)
  {
    Li = nr[i];
    sum += Li*rhor**(i+1)
  }
  Lr = -167.31*(1-exp(-2.506*rhor))+sum

  f1 = exp(0.144847*Tr-5.64493*Tr**2)
  f2 = exp(-2.8*(rhor-1)**2)-0.080738543*exp(-17.943*(rhor-0.125698)**2)
  f3 = 1+exp(60*(tau-1)+20)
  f4 = 1+exp(100*(tau-1)+15)
  Lc = 35429.6*f1*f2*(1+f2**2*(5e9*f1**4/f3+3.5*f2/f4))

  Ll = -741.112*f1**1.2*(1-exp(-(rhor/2.5)**10))

  return 0.742128e-3*(Lo+Lr+Lc+Ll)
}

void _D2O_Tension(T)
{
  /*
  """Equation for the surface tension of heavy water

  Parameters
  ----------
  T : float
    Temperature, [K]

  Returns
  -------
  σ : float
    Surface tension, [N/m]

  Notes
  ------
  Raise :class:`NotImplementedError` if input isn't in limit:

    * 269.65 ≤ T ≤ 643.847

  Examples
  --------
  >>> _D2O_Tension(298.15)
  0.07186
  >>> _D2O_Tension(573.15)
  0.01399

  References
  ----------
  IAPWS, Release on Surface Tension of Heavy Water Substance,
  http://www.iapws.org/relguide/surfd2o.pdf
  """
  */
  Tr = T/643.847
  if(269.65 <= T < 643.847)
  {
    return 1e-3*(238*(1-Tr)**1.25*(1-0.639*(1-Tr)))
  }
  else
  {
    raise NotImplementedError("Incoming out of bound")
  }
}

void _D2O_Sublimation_Pressure(T)
{
  /*
  """Sublimation Pressure correlation for heavy water

  Parameters
  ----------
  T : float
    Temperature, [K]

  Returns
  -------
  P : float
    Pressure at sublimation line, [MPa]

  Notes
  ------
  Raise :class:`NotImplementedError` if(input isn't in limit)

    * 210 ≤ T ≤ 276.969

  Examples
  --------
  >>> _Sublimation_Pressure(245)
  3.27390934e-5

  References
  ----------
  IAPWS, Revised Release on the IAPWS Formulation 2017 for the Thermodynamic
  Properties of Heavy Water, http://www.iapws.org/relguide/Heavy.html.
  """
  */
  if(210 <= T <= 276.969)
  {
    Tita = T/276.969
    suma = 0
    ai = [-0.1314226e2, 0.3212969e2]
    ti = [-1.73, -1.42]
    for(int i = 0; i < ai.length; i++)
    {
      a = ai[i];
      t = ti[i];
      suma += a*(1-Tita**t)
    }
    return exp(suma)*0.00066159
  }
  else
  {
    raise NotImplementedError("Incoming out of bound")
  }
}

void _D2O_Melting_Pressure(T, ice="Ih")
{
  /*
  """Melting Pressure correlation for heavy water

  Parameters
  ----------
  T : float
    Temperature, [K]
  ice: string
    Type of ice: Ih, III, V, VI, VII.
    Below 276.969 is a mandatory input, the ice Ih is the default value.
    Above 276.969, the ice type is unnecesary.

  Returns
  -------
  P : float
    Pressure at melting line, [MPa]

  Notes
  ------
  Raise :class:`NotImplementedError` if input isn't in limit:

    * 254.415 ≤ T ≤ 315

  Examples
  --------
  >>> _D2O__Melting_Pressure(260)
  8.947352740189152e-06
  >>> _D2O__Melting_Pressure(254, "III")
  268.6846466336108

  References
  ----------
  IAPWS, Revised Release on the Pressure along the Melting and Sublimation
  Curves of Ordinary Water Substance, http://iapws.org/relguide/MeltSub.html.
  """
  */
  if(ice == "Ih" && 254.415 <= T <= 276.969)
  {
    // Ice Ih, Eq 9
    Tita = T/276.969
    ai = [-0.30153e5, 0.692503e6]
    ti = [5.5, 8.2]
    suma = 1
    for(int i = 0; i < ai.length; i++)
    {
      a = ai[i];
      t = ti[i];
      suma += a*(1-Tita**t)
    }
    P = suma*0.00066159
  }
  else if(ice == "III" && 254.415 < T <= 258.661)
  {
    // Ice III, Eq 10
    Tita = T/254.415
    P = 222.41*(1-0.802871*(1-Tita**33))
  }
  else if(ice == "V" && 258.661 < T <= 275.748)
  {
    // Ice V, Eq 11
    Tita = T/258.661
    P = 352.19*(1-1.280388*(1-Tita**7.6))
  }
  else if((ice == "VI" && 275.748 < T <= 276.969) || 276.969 < T <= 315)
  {
    // Ice VI
    Tita = T/275.748
    P = 634.53*(1-1.276026*(1-Tita**4))
  }
  else
  {
    raise NotImplementedError("Incoming out of bound")
  }
  return P
}

void _Henry(T, gas, liquid="H2O")
{
  /*
  """Equation for the calculation of Henry's constant

  Parameters
  ----------
  T : float
    Temperature, [K]
  gas : string
    Name of gas to calculate solubility
  liquid : string
    Name of liquid solvent, can be H20 (default) or D2O

  Returns
  -------
  kw : float
    Henry's constant, [MPa]

  Notes
  -----
  The gas availables for H2O solvent are He, Ne, Ar, Kr, Xe, H2, N2, O2, CO,
  CO2, H2S, CH4, C2H6, SF6
  For D2O as solvent He, Ne, Ar, Kr, Xe, D2, CH4

  Raise :class:`NotImplementedError` if input gas or liquid are unsupported

  Examples
  --------
  >>> _Henry(500, "He")
  1.1973
  >>> _Henry(300, "D2", "D2O")
  1.6594

  References
  ----------
  IAPWS, Guideline on the Henry's Constant and Vapor-Liquid Distribution
  Constant for Gases in H2O and D2O at High Temperatures,
  http://www.iapws.org/relguide/HenGuide.html
  """
  */
  if(liquid == "D2O")
  {
    gas += "(D2O)"
  }

  limit =
  {
    "He": (273.21, 553.18),
    "Ne": (273.20, 543.36),
    "Ar": (273.19, 568.36),
    "Kr": (273.19, 525.56),
    "Xe": (273.22, 574.85),
    "H2": (273.15, 636.09),
    "N2": (278.12, 636.46),
    "O2": (274.15, 616.52),
    "CO": (278.15, 588.67),
    "CO2": (274.19, 642.66),
    "H2S": (273.15, 533.09),
    "CH4": (275.46, 633.11),
    "C2H6": (275.44, 473.46),
    "SF6": (283.14, 505.55),
    "He(D2O)": (288.15, 553.18),
    "Ne(D2O)": (288.18, 549.96),
    "Ar(D2O)": (288.30, 583.76),
    "Kr(D2O)": (288.19, 523.06),
    "Xe(D2O)": (295.39, 574.85),
    "D2(D2O)": (288.17, 581.00),
    "CH4(D2O)": (288.16, 517.46)
  }

  // Check input parameters
  if(liquid != "D2O" && liquid != "H2O")
  {
    raise NotImplementedError("Solvent liquid unsupported")
  }
  if(gas !in limit)
  {
    raise NotImplementedError("Gas unsupported")
  }

  Tmin, Tmax = limit[gas]
  if(T < Tmin || T > Tmax)
  {
    warnings.warn("Temperature out of data of correlation")
  }

  if(liquid == "D2O")
  {
    Tc = 643.847
    Pc = 21.671
  }
  else
  {
    Tc = 647.096
    Pc = 22.064
  }

  Tr = T/Tc
  tau = 1-Tr

  // Eq 4
  if(liquid == "H2O")
  {
    ai = [-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502]
    bi = [1, 1.5, 3, 3.5, 4, 7.5]
  }
  else
  {
    ai = [-7.896657, 24.73308, -27.81128, 9.355913, -9.220083]
    bi = [1, 1.89, 2, 3, 3.6]
  }
  sum = 0;
  for(int iter = 0; iter < a.length; iter++)
  {
    a = ai[iter];
    b = bi[iter];
    sum += a*tau**b;
  }
  ps = Pc*exp(1/Tr*sum)

  // Select values from Table 2
  par =
  {
    "He": (-3.52839, 7.12983, 4.47770),
    "Ne": (-3.18301, 5.31448, 5.43774),
    "Ar": (-8.40954, 4.29587, 10.52779),
    "Kr": (-8.97358, 3.61508, 11.29963),
    "Xe": (-14.21635, 4.00041, 15.60999),
    "H2": (-4.73284, 6.08954, 6.06066),
    "N2": (-9.67578, 4.72162, 11.70585),
    "O2": (-9.44833, 4.43822, 11.42005),
    "CO": (-10.52862, 5.13259, 12.01421),
    "CO2": (-8.55445, 4.01195, 9.52345),
    "H2S": (-4.51499, 5.23538, 4.42126),
    "CH4": (-10.44708, 4.66491, 12.12986),
    "C2H6": (-19.67563, 4.51222, 20.62567),
    "SF6": (-16.56118, 2.15289, 20.35440),
    "He(D2O)": (-0.72643, 7.02134, 2.04433),
    "Ne(D2O)": (-0.91999, 5.65327, 3.17247),
    "Ar(D2O)": (-7.17725, 4.48177, 9.31509),
    "Kr(D2O)": (-8.47059, 3.91580, 10.69433),
    "Xe(D2O)": (-14.46485, 4.42330, 15.60919),
    "D2(D2O)": (-5.33843, 6.15723, 6.53046),
    "CH4(D2O)": (-10.01915, 4.73368, 11.75711)
  }
  A, B, C = par[gas]

  // Eq 3
  kh = ps*exp(A/Tr+B*tau**0.355/Tr+C*Tr**-0.41*exp(tau))
  return kh
}

void _Kvalue(T, gas, liquid="H2O")
{
  /*
  """Equation for the vapor-liquid distribution constant

  Parameters
  ----------
  T : float
    Temperature, [K]
  gas : string
    Name of gas to calculate solubility
  liquid : string
    Name of liquid solvent, can be H20 (default) or D2O

  Returns
  -------
  kd : float
    Vapor-liquid distribution constant, [-]

  Notes
  -----
  The gas availables for H2O solvent are He, Ne, Ar, Kr, Xe, H2, N2, O2, CO,
  CO2, H2S, CH4, C2H6, SF6

  For D2O as solvent He, Ne, Ar, Kr, Xe, D2, CH4

  Raise :class:`NotImplementedError` if input gas or liquid are unsupported

  Examples
  --------
  >>> _Kvalue(600, "He")
  3.8019
  >>> _Kvalue(300, "D2", "D2O")
  14.3520

  References
  ----------
  IAPWS, Guideline on the Henry's Constant and Vapor-Liquid Distribution
  Constant for Gases in H2O and D2O at High Temperatures,
  http://www.iapws.org/relguide/HenGuide.html
  """
  */
  if(liquid == "D2O")
  {
    gas += "(D2O)"
  }

  limit =
  {
    "He": (273.21, 553.18),
    "Ne": (273.20, 543.36),
    "Ar": (273.19, 568.36),
    "Kr": (273.19, 525.56),
    "Xe": (273.22, 574.85),
    "H2": (273.15, 636.09),
    "N2": (278.12, 636.46),
    "O2": (274.15, 616.52),
    "CO": (278.15, 588.67),
    "CO2": (274.19, 642.66),
    "H2S": (273.15, 533.09),
    "CH4": (275.46, 633.11),
    "C2H6": (275.44, 473.46),
    "SF6": (283.14, 505.55),
    "He(D2O)": (288.15, 553.18),
    "Ne(D2O)": (288.18, 549.96),
    "Ar(D2O)": (288.30, 583.76),
    "Kr(D2O)": (288.19, 523.06),
    "Xe(D2O)": (295.39, 574.85),
    "D2(D2O)": (288.17, 581.00),
    "CH4(D2O)": (288.16, 517.46)
  }

  // Check input parameters
  if(liquid != "D2O" && liquid != "H2O")
  {
    raise NotImplementedError("Solvent liquid unsupported")
  }
  if(gas !in limit)
  {
    raise NotImplementedError("Gas unsupported")
  }

  Tmin, Tmax = limit[gas]
  if(T < Tmin || T > Tmax)
  {
    warnings.warn("Temperature out of data of correlation")
  }
  if(liquid == "D2O")
  {
    Tc = 643.847
  }
  else
  {
    Tc = 647.096
  }

  Tr = T/Tc
  tau = 1-Tr

  // Eq 6
  if(liquid == "H2O")
  {
    ci = [1.99274064, 1.09965342, -0.510839303, -1.75493479, -45.5170352, -6.7469445e5]
    di = [1/3, 2/3, 5/3, 16/3, 43/3, 110/3]
    q = -0.023767
  }
  else
  {
    ci = [2.7072, 0.58662, -1.3069, -45.663]
    di = [0.374, 1.45, 2.6, 12.3]
    q = -0.024552
  }
  sum = 0;
  for(int iter = 0; iter < c.length; iter++)
  {
    c = ci[iter];
    d = di[iter];
    sum += c*tau**d;
  }
  f = sum

  // Select values from Table 2
  par =
  {
     "He": (2267.4082, -2.9616, -3.2604, 7.8819),
     "Ne": (2507.3022, -38.6955, 110.3992, -71.9096),
     "Ar": (2310.5463, -46.7034, 160.4066, -118.3043),
     "Kr": (2276.9722, -61.1494, 214.0117, -159.0407),
     "Xe": (2022.8375, 16.7913, -61.2401, 41.9236),
     "H2": (2286.4159, 11.3397, -70.7279, 63.0631),
     "N2": (2388.8777, -14.9593, 42.0179, -29.4396),
     "O2": (2305.0674, -11.3240, 25.3224, -15.6449),
     "CO": (2346.2291, -57.6317, 204.5324, -152.6377),
     "CO2": (1672.9376, 28.1751, -112.4619, 85.3807),
     "H2S": (1319.1205, 14.1571, -46.8361, 33.2266),
     "CH4": (2215.6977, -0.1089, -6.6240, 4.6789),
     "C2H6": (2143.8121, 6.8859, -12.6084, 0),
     "SF6": (2871.7265, -66.7556, 229.7191, -172.7400),
     "He(D2O)": (2293.2474, -54.7707, 194.2924, -142.1257),
     "Ne(D2O)": (2439.6677, -93.4934, 330.7783, -243.0100),
     "Ar(D2O)": (2269.2352, -53.6321, 191.8421, -143.7659),
     "Kr(D2O)": (2250.3857, -42.0835, 140.7656, -102.7592),
     "Xe(D2O)": (2038.3656, 68.1228, -271.3390, 207.7984),
     "D2(D2O)": (2141.3214, -1.9696, 1.6136, 0),
     "CH4(D2O)": (2216.0181, -40.7666, 152.5778, -117.7430),
  }
  E, F, G, H = par[gas]

  // Eq 5
  kd = exp(q*F+E/T*f+(F+G*tau**(2./3)+H*tau)*exp((273.15-T)/100))
  return kd
}

//!/usr/bin/python
// -*- coding: utf-8 -*-
/*
"""IAPWS-IF97 standard implementation

.. image:: https://raw.githubusercontent.com/jjgomera/iapws/master/images/iapws97.png
  :alt: iapws97

The module implement the fundamental equation for the five regions (rectangular
boxes) and the backward equation (marked in grey).

:class:`IAPWS97`: Global module class with all the functionality integrated

Fundamental equations:
   * :func:`_Region1`
   * :func:`_Region2`
   * :func:`_Region3`
   * :func:`_Region4`
   * :func:`_TSat_P`
   * :func:`_PSat_T`
   * :func:`_Region5`

Backward equations:
   * :func:`_Backward1_T_Ph`
   * :func:`_Backward1_T_Ps`
   * :func:`_Backward1_P_hs`
   * :func:`_Backward2_T_Ph`
   * :func:`_Backward2_T_Ps`
   * :func:`_Backward2_P_hs`
   * :func:`_Backward3_T_Ph`
   * :func:`_Backward3_T_Ps`
   * :func:`_Backward3_P_hs`
   * :func:`_Backward3_v_Ph`
   * :func:`_Backward3_v_Ps`
   * :func:`_Backward3_v_PT`
   * :func:`_Backward4_T_hs`

Boundary equations:
   * :func:`_h13_s`
   * :func:`_h3a_s`
   * :func:`_h1_s`
   * :func:`_t_hs`
   * :func:`_PSat_h`
   * :func:`_h2ab_s`
   * :func:`_h_3ab`
   * :func:`_h2c3b_s`
   * :func:`_hab_s`
   * :func:`_hbc_P`


References:

IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
Thermodynamic Properties of Water and Steam August 2007,
http://www.iapws.org/relguide/IF97-Rev.html

IAPWS, Revised Supplementary Release on Backward Equations for Pressure
as a Function of Enthalpy and Entropy p(h,s) for Regions 1 and 2 of the IAPWS
Industrial Formulation 1997 for the Thermodynamic Properties of Water and
Steam, http://www.iapws.org/relguide/Supp-PHS12-2014.pdf

IAPWS, Revised Supplementary Release on Backward Equations for the
Functions T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS
Industrial Formulation 1997 for the Thermodynamic Properties of Water and
Steam, http://www.iapws.org/relguide/Supp-Tv%28ph,ps%293-2014.pdf

IAPWS, Revised Supplementary Release on Backward Equations p(h,s) for
Region 3, Equations as a Function of h and s for the Region Boundaries, and an
Equation Tsat(h,s) for Region 4 of the IAPWS Industrial Formulation 1997 for
the Thermodynamic Properties of Water and Steam,
http://www.iapws.org/relguide/Supp-phs3-2014.pdf

IAPWS, Revised Supplementary Release on Backward Equations for Specific
Volume as a Function of Pressure and Temperature v(p,T) for Region 3 of the
IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and
Steam, http://www.iapws.org/relguide/Supp-VPT3-2016.pdf

IAPWS, Revised Advisory Note No. 3: Thermodynamic Derivatives from IAPWS
Formulations, http://www.iapws.org/relguide/Advise3.pdf

Wagner, W; Kretzschmar, H-J: International Steam Tables: Properties of
Water and Steam Based on the Industrial Formulation IAPWS-IF97; Springer, 2008;
doi: 10.1007/978-3-540-74234-0
"""
*/

from __future__ import division
from math import sqrt, log, exp

from scipy.optimize import fsolve, newton

from ._iapws import R, Tc, Pc, rhoc, Tt, Pt, Tb, Dipole, f_acent
from ._iapws import _Viscosity, _ThCond, _Tension, _Dielectric, _Refractive
from ._utils import getphase, deriv_G, _fase


// Critic properties
sc = 4.41202148223476
hc = 2087.5468451171537

// Pmin = _PSat_T(273.15)   // Minimum pressure
Pmin = 0.000611212677444
// Ps_623 = _PSat_T(623.15)  // P Saturation at 623.15 K, boundary region 1-3
Ps_623 = 16.5291642526


// Boundary Region1-Region3
void _h13_s(s)
{
  /*
  """Define the boundary between Region 1 and 3, h=f(s)

  Parameters
  ----------
  s : float
    Specific entropy, [kJ/kgK]

  Returns
  -------
  h : float
    Specific enthalpy, [kJ/kg]

  Notes
  ------
  Raise :class:`NotImplementedError` if input isn't in limit:

    * s(100MPa,623.15K) ≤ s ≤ s'(623.15K)

  References
  ----------
  IAPWS, Revised Supplementary Release on Backward Equations p(h,s) for
  Region 3, Equations as a Function of h and s for the Region Boundaries, and
  an Equation Tsat(h,s) for Region 4 of the IAPWS Industrial Formulation 1997
  for the Thermodynamic Properties of Water and Steam,
  http://www.iapws.org/relguide/Supp-phs3-2014.pdf. Eq 7

  Examples
  --------
  >>> _h13_s(3.7)
  1632.525047
  >>> _h13_s(3.5)
  1566.104611
  """
  */
  // Check input parameters
  if(s < 3.397782955 || s > 3.77828134)
  {
    raise NotImplementedError("Incoming out of bound")
  }

  sigma = s/3.8
  I = [0, 1, 1, 3, 5, 6]
  J = [0, -2, 2, -12, -4, -3]
  n = [0.913965547600543, -0.430944856041991e-4, 0.603235694765419e2, 0.117518273082168e-17, 0.220000904781292, -0.690815545851641e2]

  suma = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    ni = n[iter];
    suma += ni * (sigma-0.884)**i * (sigma-0.864)**j
  }
  return 1700 * suma
}

// Boundary Region2-Region3
void _P23_T(T)
{
  /*
  """Define the boundary between Region 2 and 3, P=f(T)

  Parameters
  ----------
  T : float
    Temperature, [K]

  Returns
  -------
  P : float
    Pressure, [MPa]

  References
  ----------
  IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  Thermodynamic Properties of Water and Steam August 2007,
  http://www.iapws.org/relguide/IF97-Rev.html, Eq 5

  Examples
  --------
  >>> _P23_T(623.15)
  16.52916425
  """
  */
  n = [0.34805185628969e3, -0.11671859879975e1, 0.10192970039326e-2]
  return n[0]+n[1]*T+n[2]*T**2
}

void _t_P(P)
{
  /*
  """Define the boundary between Region 2 and 3, T=f(P)

  Parameters
  ----------
  P : float
    Pressure, [MPa]

  Returns
  -------
  T : float
    Temperature, [K]

  References
  ----------
  IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  Thermodynamic Properties of Water and Steam August 2007,
  http://www.iapws.org/relguide/IF97-Rev.html, Eq 5

  Examples
  --------
  >>> _t_P(16.52916425)
  623.15
  """
  */
  n = [0.10192970039326e-2, 0.57254459862746e3, 0.1391883977870e2]
  return n[1]+((P-n[2])/n[0])**0.5
}

void _t_hs(h, s)
{
  /*
  """Define the boundary between Region 2 and 3, T=f(h,s)

  Parameters
  ----------
  h : float
    Specific enthalpy, [kJ/kg]
  s : float
    Specific entropy, [kJ/kgK]

  Returns
  -------
  T : float
    Temperature, [K]

  Notes
  ------
  Raise :class:`NotImplementedError` if input isn't in limit:

    * 5.048096828 ≤ s ≤ 5.260578707
    * 2.563592004e3 ≤ h ≤ 2.812942061e3

  References
  ----------
  IAPWS, Revised Supplementary Release on Backward Equations p(h,s) for
  Region 3, Equations as a Function of h and s for the Region Boundaries, and
  an Equation Tsat(h,s) for Region 4 of the IAPWS Industrial Formulation 1997
  for the Thermodynamic Properties of Water and Steam,
  http://www.iapws.org/relguide/Supp-phs3-2014.pdf. Eq 8

  Examples
  --------
  >>> _t_hs(2600, 5.1)
  713.5259364
  >>> _t_hs(2800, 5.2)
  817.6202120
  """
  */
  // Check input parameters
  if(s < 5.048096828 || s > 5.260578707 || h < 2.563592004e3 || h > 2.812942061e3)
  {
    raise NotImplementedError("Incoming out of bound")
  }

  nu = h/3000
  sigma = s/5.3
  I = [-12, -10, -8, -4, -3, -2, -2, -2, -2, 0, 1, 1, 1, 3, 3, 5, 6, 6, 8, 8, 8, 12, 12, 14, 14]
  J = [10, 8, 3, 4, 3, -6, 2, 3, 4, 0, -3, -2, 10, -2, -1, -5, -6, -3, -8, -2, -1, -12, -1, -12, 1]
  n = [0.629096260829810e-3, -0.823453502583165e-3, 0.515446951519474e-7, -0.117565945784945e1, 0.348519684726192e1, -0.507837382408313e-11, -0.284637670005479e1, -0.236092263939673e1, 0.601492324973779e1, 0.148039650824546e1, 0.360075182221907e-3, -0.126700045009952e-1, -0.122184332521413e7, 0.149276502463272, 0.698733471798484, -0.252207040114321e-1, 0.147151930985213e-1, -0.108618917681849e1, -0.936875039816322e-3, 0.819877897570217e2, -0.182041861521835e3, 0.261907376402688e-5, -0.291626417025961e5, 0.140660774926165e-4, 0.783237062349385e7]

  suma = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    ni = n[iter];
    suma += ni * (nu-0.727)**i * (sigma-0.864)**j
  }
  return 900*suma
}

// Saturated line
void _PSat_T(T)
{
  /*
  """Define the saturated line, P=f(T)

  Parameters
  ----------
  T : float
    Temperature, [K]

  Returns
  -------
  P : float
    Pressure, [MPa]

  Notes
  ------
  Raise :class:`NotImplementedError` if input isn't in limit:

    * 273.15 ≤ T ≤ 647.096

  References
  ----------
  IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  Thermodynamic Properties of Water and Steam August 2007,
  http://www.iapws.org/relguide/IF97-Rev.html, Eq 30

  Examples
  --------
  >>> _PSat_T(500)
  2.63889776
  """
  */
  // Check input parameters
  if(T < 273.15 || T > Tc)
  {
    raise NotImplementedError("Incoming out of bound")
  }

  n = [0, 0.11670521452767E+04, -0.72421316703206E+06, -0.17073846940092E+02, 0.12020824702470E+05, -0.32325550322333E+07, 0.14915108613530E+02, -0.48232657361591E+04, 0.40511340542057E+06, -0.23855557567849E+00, 0.65017534844798E+03]
  tita = T+n[9]/(T-n[10])
  A = tita**2+n[1]*tita+n[2]
  B = n[3]*tita**2+n[4]*tita+n[5]
  C = n[6]*tita**2+n[7]*tita+n[8]
  return (2*C/(-B+(B**2-4*A*C)**0.5))**4
}

void _TSat_P(P)
{
  /*
  """Define the saturated line, T=f(P)

  Parameters
  ----------
  P : float
    Pressure, [MPa]

  Returns
  -------
  T : float
    Temperature, [K]

  Notes
  ------
  Raise :class:`NotImplementedError` if input isn't in limit:

    * 0.00061121 ≤ P ≤ 22.064

  References
  ----------
  IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  Thermodynamic Properties of Water and Steam August 2007,
  http://www.iapws.org/relguide/IF97-Rev.html, Eq 31

  Examples
  --------
  >>> _TSat_P(10)
  584.149488
  """
  */
  // Check input parameters
  if(P < 611.212677/1e6 || P > 22.064)
  {
    raise NotImplementedError("Incoming out of bound")
  }
  n = [0, 0.11670521452767E+04, -0.72421316703206E+06, -0.17073846940092E+02, 0.12020824702470E+05, -0.32325550322333E+07, 0.14915108613530E+02, -0.48232657361591E+04, 0.40511340542057E+06, -0.23855557567849E+00, 0.65017534844798E+03]
  beta = P**0.25
  E = beta**2+n[3]*beta+n[6]
  F = n[1]*beta**2+n[4]*beta+n[7]
  G = n[2]*beta**2+n[5]*beta+n[8]
  D = 2*G/(-F-(F**2-4*E*G)**0.5)
  return (n[10]+D-((n[10]+D)**2-4*(n[9]+n[10]*D))**0.5)/2
}

void _PSat_h(h)
{
  /*
  """Define the saturated line, P=f(h) for region 3

  Parameters
  ----------
  h : float
    Specific enthalpy, [kJ/kg]

  Returns
  -------
  P : float
    Pressure, [MPa]

  Notes
  ------
  Raise :class:`NotImplementedError` if input isn't in limit:

    * h'(623.15K) ≤ h ≤ h''(623.15K)

  References
  ----------
  IAPWS, Revised Supplementary Release on Backward Equations for the
  Functions T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS
  Industrial Formulation 1997 for the Thermodynamic Properties of Water and
  Steam, http://www.iapws.org/relguide/Supp-Tv%28ph,ps%293-2014.pdf, Eq 10

  Examples
  --------
  >>> _PSat_h(1700)
  17.24175718
  >>> _PSat_h(2400)
  20.18090839
  """
  */
  // Check input parameters
  hmin_Ps3 = _Region1(623.15, Ps_623)["h"]
  hmax_Ps3 = _Region2(623.15, Ps_623)["h"]
  if(h < hmin_Ps3 || h > hmax_Ps3)
  {
    raise NotImplementedError("Incoming out of bound")
  }

  nu = h/2600
  I = [0, 1, 1, 1, 1, 5, 7, 8, 14, 20, 22, 24, 28, 36]
  J = [0, 1, 3, 4, 36, 3, 0, 24, 16, 16, 3, 18, 8, 24]
  n = [0.600073641753024, -0.936203654849857e1, 0.246590798594147e2, -0.107014222858224e3, -0.915821315805768e14, -0.862332011700662e4, -0.235837344740032e2, 0.252304969384128e18, -0.389718771997719e19, -0.333775713645296e23, 0.356499469636328e11, -0.148547544720641e27, 0.330611514838798e19, 0.813641294467829e38]

  suma = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    ni = n[iter];
    suma += ni * (nu-1.02)**i * (nu-0.608)**j
  }
  return 22*suma
}

void _PSat_s(s)
{
  /*
  """Define the saturated line, P=f(s) for region 3

  Parameters
  ----------
  s : float
    Specific entropy, [kJ/kgK]

  Returns
  -------
  P : float
    Pressure, [MPa]

  Notes
  ------
  Raise :class:`NotImplementedError` if input isn't in limit:

    * s'(623.15K) ≤ s ≤ s''(623.15K)

  References
  ----------
  IAPWS, Revised Supplementary Release on Backward Equations for the
  Functions T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS
  Industrial Formulation 1997 for the Thermodynamic Properties of Water and
  Steam, http://www.iapws.org/relguide/Supp-Tv%28ph,ps%293-2014.pdf, Eq 11

  Examples
  --------
  >>> _PSat_s(3.8)
  16.87755057
  >>> _PSat_s(5.2)
  16.68968482
  """
  */
  // Check input parameters
  smin_Ps3 = _Region1(623.15, Ps_623)["s"]
  smax_Ps3 = _Region2(623.15, Ps_623)["s"]
  if(s < smin_Ps3 || s > smax_Ps3)
  {
    raise NotImplementedError("Incoming out of bound")
  }

  sigma = s/5.2
  I = [0, 1, 1, 4, 12, 12, 16, 24, 28, 32]
  J = [0, 1, 32, 7, 4, 14, 36, 10, 0, 18]
  n = [0.639767553612785, -0.129727445396014e2, -0.224595125848403e16, 0.177466741801846e7, 0.717079349571538e10, -0.378829107169011e18, -0.955586736431328e35, 0.187269814676188e24, 0.119254746466473e12, 0.110649277244882e37]

  suma = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    ni = n[iter];
    suma += ni * (sigma-1.03)**i * (sigma-0.699)**j
  }
  return 22*suma
}

void _h1_s(s)
{
  /*
  """Define the saturated line boundary between Region 1 and 4, h=f(s)

  Parameters
  ----------
  s : float
    Specific entropy, [kJ/kgK]

  Returns
  -------
  h : float
    Specific enthalpy, [kJ/kg]

  Notes
  ------
  Raise :class:`NotImplementedError` if input isn't in limit:

    * s'(273.15K) ≤ s ≤ s'(623.15K)

  References
  ----------
  IAPWS, Revised Supplementary Release on Backward Equations p(h,s) for
  Region 3, Equations as a Function of h and s for the Region Boundaries, and
  an Equation Tsat(h,s) for Region 4 of the IAPWS Industrial Formulation 1997
  for the Thermodynamic Properties of Water and Steam,
  http://www.iapws.org/relguide/Supp-phs3-2014.pdf. Eq 3

  Examples
  --------
  >>> _h1_s(1)
  308.5509647
  >>> _h1_s(3)
  1198.359754
  """
  */
  // Check input parameters
  if(s < -1.545495919e-4 || s > 3.77828134)
  {
    raise NotImplementedError("Incoming out of bound")
  }

  sigma = s/3.8
  I = [0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 7, 8, 12, 12, 14, 14, 16, 20, 20, 22, 24, 28, 32, 32]
  J = [14, 36, 3, 16, 0, 5, 4, 36, 4, 16, 24, 18, 24, 1, 4, 2, 4, 1, 22, 10, 12, 28, 8, 3, 0, 6, 8]
  n = [0.332171191705237, 0.611217706323496e-3, -0.882092478906822e1, -0.455628192543250, -0.263483840850452e-4, -0.223949661148062e2, -0.428398660164013e1, -0.616679338856916, -0.146823031104040e2, 0.284523138727299e3, -0.113398503195444e3, 0.115671380760859e4, 0.395551267359325e3, -0.154891257229285e1, 0.194486637751291e2, -0.357915139457043e1, -0.335369414148819e1, -0.664426796332460, 0.323321885383934e5, 0.331766744667084e4, -0.223501257931087e5, 0.573953875852936e7, 0.173226193407919e3, -0.363968822121321e-1, 0.834596332878346e-6, 0.503611916682674e1, 0.655444787064505e2]

  suma = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    ni = n[iter];
    suma += ni * (sigma-1.09)**i * (sigma+0.366e-4)**j
  }
  return 1700*suma
}

void _h3a_s(s)
{
  /*
  """Define the saturated line boundary between Region 4 and 3a, h=f(s)

  Parameters
  ----------
  s : float
    Specific entropy, [kJ/kgK]

  Returns
  -------
  h : float
    Specific enthalpy, [kJ/kg]

  Notes
  ------
  Raise :class:`NotImplementedError` if input isn't in limit:

    * s'(623.15K) ≤ s ≤ sc

  References
  ----------
  IAPWS, Revised Supplementary Release on Backward Equations p(h,s) for
  Region 3, Equations as a Function of h and s for the Region Boundaries, and
  an Equation Tsat(h,s) for Region 4 of the IAPWS Industrial Formulation 1997
  for the Thermodynamic Properties of Water and Steam,
  http://www.iapws.org/relguide/Supp-phs3-2014.pdf. Eq 4

  Examples
  --------
  >>> _h3a_s(3.8)
  1685.025565
  >>> _h3a_s(4.2)
  1949.352563
  """
  */
  // Check input parameters
  if(s < 3.77828134 || s > 4.41202148223476)
  {
    raise NotImplementedError("Incoming out of bound")
  }

  sigma = s/3.8
  I = [0, 0, 0, 0, 2, 3, 4, 4, 5, 5, 6, 7, 7, 7, 10, 10, 10, 32, 32]
  J = [1, 4, 10, 16, 1, 36, 3, 16, 20, 36, 4, 2, 28, 32, 14, 32, 36, 0, 6]
  n = [0.822673364673336, 0.181977213534479, -0.112000260313624e-1, -0.746778287048033e-3, -0.179046263257381, 0.424220110836657e-1, -0.341355823438768, -0.209881740853565e1, -0.822477343323596e1, -0.499684082076008e1, 0.191413958471069, 0.581062241093136e-1, -0.165505498701029e4, 0.158870443421201e4, -0.850623535172818e2, -0.317714386511207e5, -0.945890406632871e5, -0.139273847088690e-5, 0.631052532240980]

  suma = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    ni = n[iter];
    suma += ni * (sigma-1.09)**i * (sigma+0.366e-4)**j
  }
  return 1700*suma
}

void _h2ab_s(s)
{
  /*
  """Define the saturated line boundary between Region 4 and 2a-2b, h=f(s)

  Parameters
  ----------
  s : float
    Specific entropy, [kJ/kgK]

  Returns
  -------
  h : float
    Specific enthalpy, [kJ/kg]

  Notes
  ------
  Raise :class:`NotImplementedError` if input isn't in limit:

    * 5.85 ≤ s ≤ s"(273.15K)

  References
  ----------
  IAPWS, Revised Supplementary Release on Backward Equations p(h,s) for
  Region 3, Equations as a Function of h and s for the Region Boundaries, and
  an Equation Tsat(h,s) for Region 4 of the IAPWS Industrial Formulation 1997
  for the Thermodynamic Properties of Water and Steam,
  http://www.iapws.org/relguide/Supp-phs3-2014.pdf. Eq 5

  Examples
  --------
  >>> _h2ab_s(7)
  2723.729985
  >>> _h2ab_s(9)
  2511.861477
  """
  */
  // Check input parameters
  if(s < 5.85 || s > 9.155759395)
  {
    raise NotImplementedError("Incoming out of bound")
  }

  sigma1 = s/5.21
  sigma2 = s/9.2
  I = [1, 1, 2, 2, 4, 4, 7, 8, 8, 10, 12, 12, 18, 20, 24, 28, 28, 28, 28, 28, 32, 32, 32, 32, 32, 36, 36, 36, 36, 36]
  J = [8, 24, 4, 32, 1, 2, 7, 5, 12, 1, 0, 7, 10, 12, 32, 8, 12, 20, 22, 24, 2, 7, 12, 14, 24, 10, 12, 20, 22, 28]
  n = [-0.524581170928788e3, -0.926947218142218e7, -0.237385107491666e3, 0.210770155812776e11, -0.239494562010986e2, 0.221802480294197e3, -0.510472533393438e7, 0.124981396109147e7, 0.200008436996201e10, -0.815158509791035e3, -0.157612685637523e3, -0.114200422332791e11, 0.662364680776872e16, -0.227622818296144e19, -0.171048081348406e32, 0.660788766938091e16, 0.166320055886021e23, -0.218003784381501e30, -0.787276140295618e30, 0.151062329700346e32, 0.795732170300541e7, 0.131957647355347e16, -0.325097068299140e24, -0.418600611419248e26, 0.297478906557467e35, -0.953588761745473e20, 0.166957699620939e25, -0.175407764869978e33, 0.347581490626396e35, -0.710971318427851e39]

  suma = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    ni = n[iter];
    suma += ni * (1/sigma1-0.513)**i * (sigma2-0.524)**j
  }
  return 2800*exp(suma)
}

void _h2c3b_s(s)
{
  /*
  """Define the saturated line boundary between Region 4 and 2c-3b, h=f(s)

  Parameters
  ----------
  s : float
    Specific entropy, [kJ/kgK]

  Returns
  -------
  h : float
    Specific enthalpy, [kJ/kg]

  Notes
  ------
  Raise :class:`NotImplementedError` if input isn't in limit:

    * sc ≤ s ≤ 5.85

  References
  ----------
  IAPWS, Revised Supplementary Release on Backward Equations p(h,s) for
  Region 3, Equations as a Function of h and s for the Region Boundaries, and
  an Equation Tsat(h,s) for Region 4 of the IAPWS Industrial Formulation 1997
  for the Thermodynamic Properties of Water and Steam,
  http://www.iapws.org/relguide/Supp-phs3-2014.pdf. Eq 6

  Examples
  --------
  >>> _h2c3b_s(5.5)
  2687.693850
  >>> _h2c3b_s(4.5)
  2144.360448
  """
  */
  // Check input parameters
  if(s < 4.41202148223476 || s > 5.85)
  {
    raise NotImplementedError("Incoming out of bound")
  }

  sigma = s/5.9
  I = [0, 0, 0, 1, 1, 5, 6, 7, 8, 8, 12, 16, 22, 22, 24, 36]
  J = [0, 3, 4, 0, 12, 36, 12, 16, 2, 20, 32, 36, 2, 32, 7, 20]
  n = [0.104351280732769e1, -0.227807912708513e1, 0.180535256723202e1, 0.420440834792042, -0.105721244834660e6, 0.436911607493884e25, -0.328032702839753e12, -0.678686760804270e16, 0.743957464645363e4, -0.356896445355761e20, 0.167590585186801e32, -0.355028625419105e38, 0.396611982166538e12, -0.414716268484468e41, 0.359080103867382e19, -0.116994334851995e41]

  suma = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    ni = n[iter];
    suma += ni * (sigma-1.02)**i * (sigma-0.726)**j
  }
  return 2800*suma**4
}

// Region 1
void _Region1(T, P)
{
  /*
  """Basic equation for region 1

  Parameters
  ----------
  T : float
    Temperature, [K]
  P : float
    Pressure, [MPa]

  Returns
  -------
  prop : dict
    Dict with calculated properties. The available properties are:

      * v: Specific volume, [m³/kg]
      * h: Specific enthalpy, [kJ/kg]
      * s: Specific entropy, [kJ/kgK]
      * cp: Specific isobaric heat capacity, [kJ/kgK]
      * cv: Specific isocoric heat capacity, [kJ/kgK]
      * w: Speed of sound, [m/s]
      * alfav: Cubic expansion coefficient, [1/K]
      * kt: Isothermal compressibility, [1/MPa]

  References
  ----------
  IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  Thermodynamic Properties of Water and Steam August 2007,
  http://www.iapws.org/relguide/IF97-Rev.html, Eq 7

  Examples
  --------
  >>> _Region1(300,3)["v"]
  0.00100215168
  >>> _Region1(300,3)["h"]
  115.331273
  >>> _Region1(300,3)["h"]-3000*_Region1(300,3)["v"]
  112.324818
  >>> _Region1(300,80)["s"]
  0.368563852
  >>> _Region1(300,80)["cp"]
  4.01008987
  >>> _Region1(300,80)["cv"]
  3.91736606
  >>> _Region1(500,3)["w"]
  1240.71337
  >>> _Region1(500,3)["alfav"]
  0.00164118128
  >>> _Region1(500,3)["kt"]
  0.00112892188
  """
  */

  if(P < 0)
  {
    P = Pmin
  }

  I = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29, 30, 31, 32]
  J = [-2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1, 3, 17, -4, 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38, -39, -40, -41]
  n = [0.14632971213167, -0.84548187169114, -0.37563603672040e1, 0.33855169168385e1, -0.95791963387872, 0.15772038513228, -0.16616417199501e-1, 0.81214629983568e-3, 0.28319080123804e-3, -0.60706301565874e-3, -0.18990068218419e-1, -0.32529748770505e-1, -0.21841717175414e-1, -0.52838357969930e-4, -0.47184321073267e-3, -0.30001780793026e-3, 0.47661393906987e-4, -0.44141845330846e-5, -0.72694996297594e-15, -0.31679644845054e-4, -0.28270797985312e-5, -0.85205128120103e-9, -0.22425281908000e-5, -0.65171222895601e-6, -0.14341729937924e-12, -0.40516996860117e-6, -0.12734301741641e-8, -0.17424871230634e-9, -0.68762131295531e-18, 0.14478307828521e-19, 0.26335781662795e-22, -0.11947622640071e-22, 0.18228094581404e-23, -0.93537087292458e-25]
  Tr = 1386/T
  Pr = P/16.53
  g = gp = gpp = gt = gtt = gpt = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    ni = n[iter];
    g += ni * (7.1-Pr)**i * (Tr-1.222)**j
    gp -= ni*i * (7.1-Pr)**(i-1) * (Tr-1.222)**j
    gpp += ni*i*(i-1) * (7.1-Pr)**(i-2) * (Tr-1.222)**j
    gt += ni*j * (7.1-Pr)**i * (Tr-1.222)**(j-1)
    gtt += ni*j*(j-1) * (7.1-Pr)**i * (Tr-1.222)**(j-2)
    gpt -= ni*i*j * (7.1-Pr)**(i-1) * (Tr-1.222)**(j-1)
  }

  propiedades = {}
  propiedades["T"] = T
  propiedades["P"] = P
  propiedades["v"] = Pr*gp*R*T/P/1000
  propiedades["h"] = Tr*gt*R*T
  propiedades["s"] = R*(Tr*gt-g)
  propiedades["cp"] = -R*Tr**2*gtt
  propiedades["cv"] = R*(-Tr**2*gtt+(gp-Tr*gpt)**2/gpp)
  propiedades["w"] = sqrt(R*T*1000*gp**2/((gp-Tr*gpt)**2/(Tr**2*gtt)-gpp))
  propiedades["alfav"] = (1-Tr*gpt/gp)/T
  propiedades["kt"] = -Pr*gpp/gp/P
  propiedades["region"] = 1
  propiedades["x"] = 0
  return propiedades
}

void _Backward1_T_Ph(P, h)
{
  /*
  """
  Backward equation for region 1, T=f(P,h)

  Parameters
  ----------
  P : float
    Pressure, [MPa]
  h : float
    Specific enthalpy, [kJ/kg]

  Returns
  -------
  T : float
    Temperature, [K]

  References
  ----------
  IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  Thermodynamic Properties of Water and Steam August 2007,
  http://www.iapws.org/relguide/IF97-Rev.html, Eq 11

  Examples
  --------
  >>> _Backward1_T_Ph(3,500)
  391.798509
  >>> _Backward1_T_Ph(80,1500)
  611.041229
  """
  */
  I = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 4, 5, 6]
  J = [0, 1, 2, 6, 22, 32, 0, 1, 2, 3, 4, 10, 32, 10, 32, 10, 32, 32, 32, 32]
  n = [-0.23872489924521e3, 0.40421188637945e3, 0.11349746881718e3, -0.58457616048039e1, -0.15285482413140e-3, -0.10866707695377e-5, -0.13391744872602e2, 0.43211039183559e2, -0.54010067170506e2, 0.30535892203916e2, -0.65964749423638e1, 0.93965400878363e-2, 0.11573647505340e-6, -0.25858641282073e-4, -0.40644363084799e-8, 0.66456186191635e-7, 0.80670734103027e-10, -0.93477771213947e-12, 0.58265442020601e-14, -0.15020185953503e-16]

  Pr = P/1
  nu = h/2500
  T = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    ni = n[iter];
    T += ni * Pr**i * (nu+1)**j
  }
  return T
}

void _Backward1_T_Ps(P, s)
{
  /*
  """Backward equation for region 1, T=f(P,s)

  Parameters
  ----------
  P : float
    Pressure, [MPa]
  s : float
    Specific entropy, [kJ/kgK]

  Returns
  -------
  T : float
    Temperature, [K]

  References
  ----------
  IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  Thermodynamic Properties of Water and Steam August 2007,
  http://www.iapws.org/relguide/IF97-Rev.html, Eq 13

  Examples
  --------
  >>> _Backward1_T_Ps(3,0.5)
  307.842258
  >>> _Backward1_T_Ps(80,3)
  565.899909
  """
  */
  I = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 4]
  J = [0, 1, 2, 3, 11, 31, 0, 1, 2, 3, 12, 31, 0, 1, 2, 9, 31, 10, 32, 32]
  n = [0.17478268058307e3, 0.34806930892873e2, 0.65292584978455e1, 0.33039981775489, -0.19281382923196e-6, -0.24909197244573e-22, -0.26107636489332, 0.22592965981586, -0.64256463395226e-1, 0.78876289270526e-2, 0.35672110607366e-9, 0.17332496994895e-23, 0.56608900654837e-3, -0.32635483139717e-3, 0.44778286690632e-4, -0.51322156908507e-9, -0.42522657042207e-25, 0.26400441360689e-12, 0.78124600459723e-28, -0.30732199903668e-30]

  Pr = P/1
  sigma = s/1
  T = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    ni = n[iter];
    T += ni * Pr**i * (sigma+2)**j
  }
  return T
}

void _Backward1_P_hs(h, s)
{
  /*
  """Backward equation for region 1, P=f(h,s)

  Parameters
  ----------
  h : float
    Specific enthalpy, [kJ/kg]
  s : float
    Specific entropy, [kJ/kgK]

  Returns
  -------
  P : float
    Pressure, [MPa]

  References
  ----------
  IAPWS, Revised Supplementary Release on Backward Equations for Pressure
  as a Function of Enthalpy and Entropy p(h,s) for Regions 1 and 2 of the
  IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of
  Water and Steam, http://www.iapws.org/relguide/Supp-PHS12-2014.pdf, Eq 1

  Examples
  --------
  >>> _Backward1_P_hs(0.001,0)
  0.0009800980612
  >>> _Backward1_P_hs(90,0)
  91.92954727
  >>> _Backward1_P_hs(1500,3.4)
  58.68294423
  """
  */
  I = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 4, 4, 5]
  J = [0, 1, 2, 4, 5, 6, 8, 14, 0, 1, 4, 6, 0, 1, 10, 4, 1, 4, 0]
  n = [-0.691997014660582, -0.183612548787560e2, -0.928332409297335e1, 0.659639569909906e2, -0.162060388912024e2, 0.450620017338667e3, 0.854680678224170e3, 0.607523214001162e4, 0.326487682621856e2, -0.269408844582931e2, -0.319947848334300e3, -0.928354307043320e3, 0.303634537455249e2, -0.650540422444146e2, -0.430991316516130e4, -0.747512324096068e3, 0.730000345529245e3, 0.114284032569021e4, -0.436407041874559e3]

  nu = h/3400
  sigma = s/7.6
  P = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    ni = n[iter];
    P += ni * (nu+0.05)**i * (sigma+0.05)**j
  }
  return 100*P
}

// Region 2
void _Region2(T, P)
{
  /*
  """Basic equation for region 2

  Parameters
  ----------
  T : float
    Temperature, [K]
  P : float
    Pressure, [MPa]

  Returns
  -------
  prop : dict
    Dict with calculated properties. The available properties are:

      * v: Specific volume, [m³/kg]
      * h: Specific enthalpy, [kJ/kg]
      * s: Specific entropy, [kJ/kgK]
      * cp: Specific isobaric heat capacity, [kJ/kgK]
      * cv: Specific isocoric heat capacity, [kJ/kgK]
      * w: Speed of sound, [m/s]
      * alfav: Cubic expansion coefficient, [1/K]
      * kt: Isothermal compressibility, [1/MPa]

  References
  ----------
  IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  Thermodynamic Properties of Water and Steam August 2007,
  http://www.iapws.org/relguide/IF97-Rev.html, Eq 15-17

  Examples
  --------
  >>> _Region2(700,30)["v"]
  0.00542946619
  >>> _Region2(700,30)["h"]
  2631.49474
  >>> _Region2(700,30)["h"]-30000*_Region2(700,30)["v"]
  2468.61076
  >>> _Region2(700,0.0035)["s"]
  10.1749996
  >>> _Region2(700,0.0035)["cp"]
  2.08141274
  >>> _Region2(700,0.0035)["cv"]
  1.61978333
  >>> _Region2(300,0.0035)["w"]
  427.920172
  >>> _Region2(300,0.0035)["alfav"]
  0.00337578289
  >>> _Region2(300,0.0035)["kt"]
  286.239651
  """
  */
  if(P < 0)
  {
    P = Pmin
  }

  Tr = 540/T
  Pr = P/1

  go, gop, gopp, got, gott, gopt = Region2_cp0(Tr, Pr)

  Ir = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 5, 6, 6, 6, 7, 7, 7, 8, 8, 9, 10, 10, 10, 16, 16, 18, 20, 20, 20, 21, 22, 23, 24, 24, 24]
  Jr = [0, 1, 2, 3, 6, 1, 2, 4, 7, 36, 0, 1, 3, 6, 35, 1, 2, 3, 7, 3, 16, 35, 0, 11, 25, 8, 36, 13, 4, 10, 14, 29, 50, 57, 20, 35, 48, 21, 53, 39, 26, 40, 58]
  nr = [-0.0017731742473212999, -0.017834862292357999, -0.045996013696365003, -0.057581259083432, -0.050325278727930002, -3.3032641670203e-05, -0.00018948987516315, -0.0039392777243355001, -0.043797295650572998, -2.6674547914087001e-05, 2.0481737692308999e-08, 4.3870667284435001e-07, -3.2277677238570002e-05, -0.0015033924542148, -0.040668253562648998, -7.8847309559367001e-10, 1.2790717852285001e-08, 4.8225372718507002e-07, 2.2922076337661001e-06, -1.6714766451061001e-11, -0.0021171472321354998, -23.895741934103999, -5.9059564324270004e-18, -1.2621808899101e-06, -0.038946842435739003, 1.1256211360459e-11, -8.2311340897998004, 1.9809712802088e-08, 1.0406965210174e-19, -1.0234747095929e-13, -1.0018179379511e-09, -8.0882908646984998e-11, 0.10693031879409, -0.33662250574170999, 8.9185845355420999e-25, 3.0629316876231997e-13, -4.2002467698208001e-06, -5.9056029685639003e-26, 3.7826947613457002e-06, -1.2768608934681e-15, 7.3087610595061e-29, 5.5414715350778001e-17, -9.4369707241209998e-07]

  gr = grp = grpp = grt = grtt = grpt = 0
  for(int iter = 0; iter < Ir.length; iter++)
  {
    i = Ir[iter];
    j = Jr[iter];
    ni = nr[iter];
    gr += ni * Pr**i * (Tr-0.5)**j
    grp += ni*i * Pr**(i-1) * (Tr-0.5)**j
    grpp += ni*i*(i-1) * Pr**(i-2) * (Tr-0.5)**j
    grt += ni*j * Pr**i * (Tr-0.5)**(j-1)
    grtt += ni*j*(j-1) * Pr**i * (Tr-0.5)**(j-2)
    grpt += ni*i*j * Pr**(i-1) * (Tr-0.5)**(j-1)
  }

  propiedades = {}
  propiedades["T"] = T
  propiedades["P"] = P
  propiedades["v"] = Pr*(gop+grp)*R*T/P/1000
  propiedades["h"] = Tr*(got+grt)*R*T
  propiedades["s"] = R*(Tr*(got+grt)-(go+gr))
  propiedades["cp"] = -R*Tr**2*(gott+grtt)
  propiedades["cv"] = R*(-Tr**2*(gott+grtt)-(1+Pr*grp-Tr*Pr*grpt)**2 /
               (1-Pr**2*grpp))
  propiedades["w"] = (R*T*1000*(1+2*Pr*grp+Pr**2*grp**2)/(1-Pr**2*grpp+(
    1+Pr*grp-Tr*Pr*grpt)**2/Tr**2/(gott+grtt)))**0.5
  propiedades["alfav"] = (1+Pr*grp-Tr*Pr*grpt)/(1+Pr*grp)/T
  propiedades["kt"] = (1-Pr**2*grpp)/(1+Pr*grp)/P
  propiedades["region"] = 2
  propiedades["x"] = 1
  return propiedades
}

void Region2_cp0(Tr, Pr)
{
  /*
  """Ideal properties for Region 2

  Parameters
  ----------
  Tr : float
    Reduced temperature, [-]
  Pr : float
    Reduced pressure, [-]

  Returns
  -------
  prop : array
    Array with ideal Gibbs energy partial derivatives:

      * g: Ideal Specific Gibbs energy [kJ/kg]
      * gp: ∂g/∂P|T
      * gpp: ∂²g/∂P²|T
      * gt: ∂g/∂T|P
      * gtt: ∂²g/∂T²|P
      * gpt: ∂²g/∂T∂P

  References
  ----------
  IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  Thermodynamic Properties of Water and Steam August 2007,
  http://www.iapws.org/relguide/IF97-Rev.html, Eq 16

  """
  */
  Jo = [0, 1, -5, -4, -3, -2, -1, 2, 3]
  no = [-0.96927686500217E+01, 0.10086655968018E+02, -0.56087911283020E-02, 0.71452738081455E-01, -0.40710498223928E+00, 0.14240819171444E+01, -0.43839511319450E+01, -0.28408632460772E+00, 0.21268463753307E-01]
  go = log(Pr)
  gop = Pr**-1
  gopp = -Pr**-2
  got = gott = gopt = 0
  for(int i = 0; i < Jo.length; i++)
  {
    j = Jo[i];
    ni= no[i];
    go += ni * Tr**j
    got += ni*j * Tr**(j-1)
    gott += ni*j*(j-1) * Tr**(j-2)
  }
  return go, gop, gopp, got, gott, gopt
}

void _P_2bc(h)
{
  /*
  """Define the boundary between Region 2b and 2c, P=f(h)

  Parameters
  ----------
  h : float
    Specific enthalpy, [kJ/kg]

  Returns
  -------
  P : float
    Pressure, [MPa]

  References
  ----------
  IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  Thermodynamic Properties of Water and Steam August 2007,
  http://www.iapws.org/relguide/IF97-Rev.html, Eq 20

  Examples
  --------
  >>> _P_2bc(3516.004323)
  100.0
  """
  */
  return 905.84278514723-0.67955786399241*h+1.2809002730136e-4*h**2
}

void _hbc_P(P)
{
  /*
  """Define the boundary between Region 2b and 2c, h=f(P)

  Parameters
  ----------
  P : float
    Pressure, [MPa]

  Returns
  -------
  h : float
    Specific enthalpy, [kJ/kg]

  References
  ----------
  IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  Thermodynamic Properties of Water and Steam August 2007,
  http://www.iapws.org/relguide/IF97-Rev.html, Eq 21

  Examples
  --------
  >>> _hbc_P(100)
  3516.004323
  """
  */
  return 0.26526571908428e4+((P-0.45257578905948e1)/1.2809002730136e-4)**0.5
}

void _hab_s(s)
{
  /*
  """Define the boundary between Region 2a and 2b, h=f(s)

  Parameters
  ----------
  s : float
    Specific entropy, [kJ/kgK]

  Returns
  -------
  h : float
    Specific enthalpy, [kJ/kg]

  References
  ----------
  IAPWS, Revised Supplementary Release on Backward Equations for Pressure
  as a Function of Enthalpy and Entropy p(h,s) for Regions 1 and 2 of the
  IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of
  Water and Steam, http://www.iapws.org/relguide/Supp-PHS12-2014.pdf, Eq 2

  Examples
  --------
  >>> _hab_s(7)
  3376.437884
  """
  */
  smin = _Region2(_TSat_P(4), 4)["s"]
  smax = _Region2(1073.15, 4)["s"]
  if(s < smin)
  {
    h = 0
  }
  else if(s > smax)
  {
    h = 5000
  }
  else
  {
    h = -0.349898083432139e4 + 0.257560716905876e4*s - \
      0.421073558227969e3*s**2+0.276349063799944e2*s**3
  }
  return h
}

void _Backward2a_T_Ph(P, h)
{
  /*
  """Backward equation for region 2a, T=f(P,h)

  Parameters
  ----------
  P : float
    Pressure, [MPa]
  h : float
    Specific enthalpy, [kJ/kg]

  Returns
  -------
  T : float
    Temperature, [K]

  References
  ----------
  IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  Thermodynamic Properties of Water and Steam August 2007,
  http://www.iapws.org/relguide/IF97-Rev.html, Eq 22

  Examples
  --------
  >>> _Backward2a_T_Ph(0.001,3000)
  534.433241
  >>> _Backward2a_T_Ph(3,4000)
  1010.77577
  """
  */
  I = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 7]
  J = [0, 1, 2, 3, 7, 20, 0, 1, 2, 3, 7, 9, 11, 18, 44, 0, 2, 7, 36, 38, 40, 42, 44, 24, 44, 12, 32, 44, 32, 36, 42, 34, 44, 28]
  n = [0.10898952318288e4, 0.84951654495535e3, -0.10781748091826e3, 0.33153654801263e2, -0.74232016790248e1, 0.11765048724356e2, 0.18445749355790e1, -0.41792700549624e1, 0.62478196935812e1, -0.17344563108114e2, -0.20058176862096e3, 0.27196065473796e3, -0.45511318285818e3, 0.30919688604755e4, 0.25226640357872e6, -0.61707422868339e-2, -0.31078046629583, 0.11670873077107e2, 0.12812798404046e9, -0.98554909623276e9, 0.28224546973002e10, -0.35948971410703e10, 0.17227349913197e10, -0.13551334240775e5, 0.12848734664650e8, 0.13865724283226e1, 0.23598832556514e6, -0.13105236545054e8, 0.73999835474766e4, -0.55196697030060e6, 0.37154085996233e7, 0.19127729239660e5, -0.41535164835634e6, -0.62459855192507e2]

  Pr = P/1
  nu = h/2000
  T = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    ni = n[iter];
    T += ni * Pr**i * (nu-2.1)**j
  }
  return T
}

void _Backward2b_T_Ph(P, h)
{
  /*
  """Backward equation for region 2b, T=f(P,h)

  Parameters
  ----------
  P : float
    Pressure, [MPa]
  h : float
    Specific enthalpy, [kJ/kg]

  Returns
  -------
  T : float
    Temperature, [K]

  References
  ----------
  IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  Thermodynamic Properties of Water and Steam August 2007,
  http://www.iapws.org/relguide/IF97-Rev.html, Eq 23

  Examples
  --------
  >>> _Backward2b_T_Ph(5,4000)
  1015.31583
  >>> _Backward2b_T_Ph(25,3500)
  875.279054
  """
  */
  I = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 6, 7, 7, 9, 9]
  J = [0, 1, 2, 12, 18, 24, 28, 40, 0, 2, 6, 12, 18, 24, 28, 40, 2, 8, 18, 40, 1, 2, 12, 24, 2, 12, 18, 24, 28, 40, 18, 24, 40, 28, 2, 28, 1, 40]
  n = [0.14895041079516e4, 0.74307798314034e3, -0.97708318797837e2, 0.24742464705674e1, -0.63281320016026, 0.11385952129658e1, -0.47811863648625, 0.85208123431544e-2, 0.93747147377932, 0.33593118604916e1, 0.33809355601454e1, 0.16844539671904, 0.73875745236695, -0.47128737436186, 0.15020273139707, -0.21764114219750e-2, -0.21810755324761e-1, -0.10829784403677, -0.46333324635812e-1, 0.71280351959551e-4, 0.11032831789999e-3, 0.18955248387902e-3, 0.30891541160537e-2, 0.13555504554949e-2, 0.28640237477456e-6, -0.10779857357512e-4, -0.76462712454814e-4, 0.14052392818316e-4, -0.31083814331434e-4, -0.10302738212103e-5, 0.28217281635040e-6, 0.12704902271945e-5, 0.73803353468292e-7, -0.11030139238909e-7, -0.81456365207833e-13, -0.25180545682962e-10, -0.17565233969407e-17, 0.86934156344163e-14]

  Pr = P/1
  nu = h/2000
  T = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    ni = n[iter];
    T += ni * (Pr-2)**i * (nu-2.6)**j
  }
  return T
}

void _Backward2c_T_Ph(P, h)
{
  /*
  """Backward equation for region 2c, T=f(P,h)

  Parameters
  ----------
  P : float
    Pressure, [MPa]
  h : float
    Specific enthalpy, [kJ/kg]

  Returns
  -------
  T : float
    Temperature, [K]

  References
  ----------
  IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  Thermodynamic Properties of Water and Steam August 2007,
  http://www.iapws.org/relguide/IF97-Rev.html, Eq 24

  Examples
  --------
  >>> _Backward2c_T_Ph(40,2700)
  743.056411
  >>> _Backward2c_T_Ph(60,3200)
  882.756860
  """
  */
  I = [-7, -7, -6, -6, -5, -5, -2, -2, -1, -1, 0, 0, 1, 1, 2, 6, 6, 6, 6, 6, 6, 6, 6]
  J = [0, 4, 0, 2, 0, 2, 0, 1, 0, 2, 0, 1, 4, 8, 4, 0, 1, 4, 10, 12, 16, 20, 22]
  n = [-0.32368398555242e13, 0.73263350902181e13, 0.35825089945447e12, -0.58340131851590e12, -0.10783068217470e11, 0.20825544563171e11, 0.61074783564516e6, 0.85977722535580e6, -0.25745723604170e5, 0.31081088422714e5, 0.12082315865936e4, 0.48219755109255e3, 0.37966001272486e1, -0.10842984880077e2, -0.45364172676660e-1, 0.14559115658698e-12, 0.11261597407230e-11, -0.17804982240686e-10, 0.12324579690832e-6, -0.11606921130984e-5, 0.27846367088554e-4, -0.59270038474176e-3, 0.12918582991878e-2]

  Pr = P/1
  nu = h/2000
  T = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    ni = n[iter];
    T += ni * (Pr+25)**i * (nu-1.8)**j
  }
  return T
}

void _Backward2_T_Ph(P, h)
{
  /*
  """Backward equation for region 2, T=f(P,h)

  Parameters
  ----------
  P : float
    Pressure, [MPa]
  h : float
    Specific enthalpy, [kJ/kg]

  Returns
  -------
  T : float
    Temperature, [K]
  """
  */
  if(P <= 4)
  {
    T = _Backward2a_T_Ph(P, h)
  }
  else if(4 < P <= 6.546699678)
  {
    T = _Backward2b_T_Ph(P, h)
  }
  else
  {
    hf = _hbc_P(P)
    if(h >= hf)
    {
      T = _Backward2b_T_Ph(P, h)
    }
    else
    {
      T = _Backward2c_T_Ph(P, h)
    }
  }

  if(P <= 22.064)
  {
    Tsat = _TSat_P(P)
    T = max(Tsat, T)
  }
  return T
}

void _Backward2a_T_Ps(P, s)
{
  /*
  """Backward equation for region 2a, T=f(P,s)

  Parameters
  ----------
  P : float
    Pressure, [MPa]
  s : float
    Specific entropy, [kJ/kgK]

  Returns
  -------
  T : float
    Temperature, [K]

  References
  ----------
  IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  Thermodynamic Properties of Water and Steam August 2007,
  http://www.iapws.org/relguide/IF97-Rev.html, Eq 25

  Examples
  --------
  >>> _Backward2a_T_Ps(0.1,7.5)
  399.517097
  >>> _Backward2a_T_Ps(2.5,8)
  1039.84917
  """
  */
  I = [-1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.25, -1.25, -1.25, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -0.75, -0.75, -0.5, -0.5, -0.5, -0.5, -0.25, -0.25, -0.25, -0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.75, 0.75, 0.75, 0.75, 1.0, 1.0, 1.25, 1.25, 1.5, 1.5]
  J = [-24, -23, -19, -13, -11, -10, -19, -15, -6, -26, -21, -17, -16, -9, -8, -15, -14, -26, -13, -9, -7, -27, -25, -11, -6, 1, 4, 8, 11, 0, 1, 5, 6, 10, 14, 16, 0, 4, 9, 17, 7, 18, 3, 15, 5, 18]
  n = [-0.39235983861984e6, 0.51526573827270e6, 0.40482443161048e5, -0.32193790923902e3, 0.96961424218694e2, -0.22867846371773e2, -0.44942914124357e6, -0.50118336020166e4, 0.35684463560015, 0.44235335848190e5, -0.13673388811708e5, 0.42163260207864e6, 0.22516925837475e5, 0.47442144865646e3, -0.14931130797647e3, -0.19781126320452e6, -0.23554399470760e5, -0.19070616302076e5, 0.55375669883164e5, 0.38293691437363e4, -0.60391860580567e3, 0.19363102620331e4, 0.42660643698610e4, -0.59780638872718e4, -0.70401463926862e3, 0.33836784107553e3, 0.20862786635187e2, 0.33834172656196e-1, -0.43124428414893e-4, 0.16653791356412e3, -0.13986292055898e3, -0.78849547999872, 0.72132411753872e-1, -0.59754839398283e-2, -0.12141358953904e-4, 0.23227096733871e-6, -0.10538463566194e2, 0.20718925496502e1, -0.72193155260427e-1, 0.20749887081120e-6, -0.18340657911379e-1, 0.29036272348696e-6, 0.21037527893619, 0.25681239729999e-3, -0.12799002933781e-1, -0.82198102652018e-5]

  Pr = P/1
  sigma = s/2
  T = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    ni = n[iter];
    T += ni * Pr**i * (sigma-2)**j
  }
  return T
}

void _Backward2b_T_Ps(P, s)
{
  /*
  """Backward equation for region 2b, T=f(P,s)

  Parameters
  ----------
  P : float
    Pressure, [MPa]
  s : float
    Specific entropy, [kJ/kgK]

  Returns
  -------
  T : float
    Temperature, [K]

  References
  ----------
  IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  Thermodynamic Properties of Water and Steam August 2007,
  http://www.iapws.org/relguide/IF97-Rev.html, Eq 26

  Examples
  --------
  >>> _Backward2b_T_Ps(8,6)
  600.484040
  >>> _Backward2b_T_Ps(90,6)
  1038.01126
  """
  */
  I = [-6, -6, -5, -5, -4, -4, -4, -3, -3, -3, -3, -2, -2, -2, -2, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5]
  J = [0, 11, 0, 11, 0, 1, 11, 0, 1, 11, 12, 0, 1, 6, 10, 0, 1, 5, 8, 9, 0, 1, 2, 4, 5, 6, 9, 0, 1, 2, 3, 7, 8, 0, 1, 5, 0, 1, 3, 0, 1, 0, 1, 2]
  n = [0.31687665083497e6, 0.20864175881858e2, -0.39859399803599e6, -0.21816058518877e2, 0.22369785194242e6, -0.27841703445817e4, 0.99207436071480e1, -0.75197512299157e5, 0.29708605951158e4, -0.34406878548526e1, 0.38815564249115, 0.17511295085750e5, -0.14237112854449e4, 0.10943803364167e1, 0.89971619308495, -0.33759740098958e4, 0.47162885818355e3, -0.19188241993679e1, 0.41078580492196, -0.33465378172097, 0.13870034777505e4, -0.40663326195838e3, 0.41727347159610e2, 0.21932549434532e1, -0.10320050009077e1, 0.35882943516703, 0.52511453726066e-2, 0.12838916450705e2, -0.28642437219381e1, 0.56912683664855, -0.99962954584931e-1, -0.32632037778459e-2, 0.23320922576723e-3, -0.15334809857450, 0.29072288239902e-1, 0.37534702741167e-3, 0.17296691702411e-2, -0.38556050844504e-3, -0.35017712292608e-4, -0.14566393631492e-4, 0.56420857267269e-5, 0.41286150074605e-7, -0.20684671118824e-7, 0.16409393674725e-8]

  Pr = P/1
  sigma = s/0.7853
  T = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    ni = n[iter];
    T += ni * Pr**i * (10-sigma)**j
  }
  return T
}

void _Backward2c_T_Ps(P, s)
{
  /*
  """Backward equation for region 2c, T=f(P,s)

  Parameters
  ----------
  P : float
    Pressure, [MPa]
  s : float
    Specific entropy, [kJ/kgK]

  Returns
  -------
  T : float
    Temperature, [K]

  References
  ----------
  IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  Thermodynamic Properties of Water and Steam August 2007,
  http://www.iapws.org/relguide/IF97-Rev.html, Eq 27

  Examples
  --------
  >>> _Backward2c_T_Ps(20,5.75)
  697.992849
  >>> _Backward2c_T_Ps(80,5.75)
  949.017998
  """
  */
  I = [-2, -2, -1, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 7, 7, 7, 7, 7]
  J = [0, 1, 0, 0, 1, 2, 3, 0, 1, 3, 4, 0, 1, 2, 0, 1, 5, 0, 1, 4, 0, 1, 2, 0, 1, 0, 1, 3, 4, 5]
  n = [0.90968501005365e3, 0.24045667088420e4, -0.59162326387130e3, 0.54145404128074e3, -0.27098308411192e3, 0.97976525097926e3, -0.46966772959435e3, 0.14399274604723e2, -0.19104204230429e2, 0.53299167111971e1, -0.21252975375934e2, -0.31147334413760, 0.60334840894623, -0.42764839702509e-1, 0.58185597255259e-2, -0.14597008284753e-1, 0.56631175631027e-2, -0.76155864584577e-4, 0.22440342919332e-3, -0.12561095013413e-4, 0.63323132660934e-6, -0.20541989675375e-5, 0.36405370390082e-7, -0.29759897789215e-8, 0.10136618529763e-7, 0.59925719692351e-11, -0.20677870105164e-10, -0.20874278181886e-10, 0.10162166825089e-9, -0.16429828281347e-9]

  Pr = P/1
  sigma = s/2.9251
  T = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    ni = n[iter];
    T += ni * Pr**i * (2-sigma)**j
  }
  return T
}

void _Backward2_T_Ps(P, s)
{
  /*
  """Backward equation for region 2, T=f(P,s)

  Parameters
  ----------
  P : float
    Pressure, [MPa]
  s : float
    Specific entropy, [kJ/kgK]

  Returns
  -------
  T : float
    Temperature, [K]
  """
  */

  if(P <= 4)
  {
    T = _Backward2a_T_Ps(P, s)
  }
  else if(s >= 5.85)
  {
    T = _Backward2b_T_Ps(P, s)
  }
  else
  {
    T = _Backward2c_T_Ps(P, s)
  }

  if(P <= 22.064)
  {
    Tsat = _TSat_P(P)
    T = max(Tsat, T)
  }
  return T
}

void _Backward2a_P_hs(h, s)
{
  /*
  """Backward equation for region 2a, P=f(h,s)

  Parameters
  ----------
  h : float
    Specific enthalpy, [kJ/kg]
  s : float
    Specific entropy, [kJ/kgK]

  Returns
  -------
  P : float
    Pressure, [MPa]

  References
  ----------
  IAPWS, Revised Supplementary Release on Backward Equations for Pressure
  as a Function of Enthalpy and Entropy p(h,s) for Regions 1 and 2 of the
  IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of
  Water and Steam, http://www.iapws.org/relguide/Supp-PHS12-2014.pdf, Eq 3

  Examples
  --------
  >>> _Backward2a_P_hs(2800,6.5)
  1.371012767
  >>> _Backward2a_P_hs(2800,9.5)
  0.001879743844
  >>> _Backward2a_P_hs(4100,9.5)
  0.1024788997
  """
  */

  I = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3, 4, 5, 5, 6, 7]
  J = [1, 3, 6, 16, 20, 22, 0, 1, 2, 3, 5, 6, 10, 16, 20, 22, 3, 16, 20, 0, 2, 3, 6, 16, 16, 3, 16, 3, 1]
  n = [-0.182575361923032e-1, -0.125229548799536, 0.592290437320145, 0.604769706185122e1, 0.238624965444474e3, -0.298639090222922e3, 0.512250813040750e-1, -0.437266515606486, 0.413336902999504, -0.516468254574773e1, -0.557014838445711e1, 0.128555037824478e2, 0.114144108953290e2, -0.119504225652714e3, -0.284777985961560e4, 0.431757846408006e4, 0.112894040802650e1, 0.197409186206319e4, 0.151612444706087e4, 0.141324451421235e-1, 0.585501282219601, -0.297258075863012e1, 0.594567314847319e1, -0.623656565798905e4, 0.965986235133332e4, 0.681500934948134e1, -0.633207286824489e4, -0.558919224465760e1, 0.400645798472063e-1]

  nu = h/4200
  sigma = s/12
  suma = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    ni = n[iter];
    suma += ni * (nu-0.5)**i * (sigma-1.2)**j
  }
  return 4*suma**4
}

void _Backward2b_P_hs(h, s)
{
  /*
  """Backward equation for region 2b, P=f(h,s)

  Parameters
  ----------
  h : float
    Specific enthalpy, [kJ/kg]
  s : float
    Specific entropy, [kJ/kgK]

  Returns
  -------
  P : float
    Pressure, [MPa]

  References
  ----------
  IAPWS, Revised Supplementary Release on Backward Equations for Pressure
  as a Function of Enthalpy and Entropy p(h,s) for Regions 1 and 2 of the
  IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of
  Water and Steam, http://www.iapws.org/relguide/Supp-PHS12-2014.pdf, Eq 4

  Examples
  --------
  >>> _Backward2b_P_hs(2800,6)
  4.793911442
  >>> _Backward2b_P_hs(3600,6)
  83.95519209
  >>> _Backward2b_P_hs(3600,7)
  7.527161441
  """
  */
  I = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 6, 6, 6, 7, 7, 8, 8, 8, 8, 12, 14]
  J = [0, 1, 2, 4, 8, 0, 1, 2, 3, 5, 12, 1, 6, 18, 0, 1, 7, 12, 1, 16, 1, 12, 1, 8, 18, 1, 16, 1, 3, 14, 18, 10, 16]
  n = [0.801496989929495e-1, -0.543862807146111, 0.337455597421283, 0.890555451157450e1, 0.313840736431485e3, 0.797367065977789, -0.121616973556240e1, 0.872803386937477e1, -0.169769781757602e2, -0.186552827328416e3, 0.951159274344237e5, -0.189168510120494e2, -0.433407037194840e4, 0.543212633012715e9, 0.144793408386013, 0.128024559637516e3, -0.672309534071268e5, 0.336972380095287e8, -0.586634196762720e3, -0.221403224769889e11, 0.171606668708389e4, -0.570817595806302e9, -0.312109693178482e4, -0.207841384633010e7, 0.305605946157786e13, 0.322157004314333e4, 0.326810259797295e12, -0.144104158934487e4, 0.410694867802691e3, 0.109077066873024e12, -0.247964654258893e14, 0.188801906865134e10, -0.123651009018773e15]

  nu = h/4100
  sigma = s/7.9
  suma = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    ni = n[iter];
    suma += ni * (nu-0.6)**i * (sigma-1.01)**j
  }
  return 100*suma**4
}

void _Backward2c_P_hs(h, s)
{
  /*
  """Backward equation for region 2c, P=f(h,s)

  Parameters
  ----------
  h : float
    Specific enthalpy, [kJ/kg]
  s : float
    Specific entropy, [kJ/kgK]

  Returns
  -------
  P : float
    Pressure, [MPa]

  References
  ----------
  IAPWS, Revised Supplementary Release on Backward Equations for Pressure
  as a Function of Enthalpy and Entropy p(h,s) for Regions 1 and 2 of the
  IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of
  Water and Steam, http://www.iapws.org/relguide/Supp-PHS12-2014.pdf, Eq 5

  Examples
  --------
  >>> _Backward2c_P_hs(2800,5.1)
  94.39202060
  >>> _Backward2c_P_hs(2800,5.8)
  8.414574124
  >>> _Backward2c_P_hs(3400,5.8)
  83.76903879
  """
  */

  I = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 5, 5, 5, 5, 6, 6, 10, 12, 16]
  J = [0, 1, 2, 3, 4, 8, 0, 2, 5, 8, 14, 2, 3, 7, 10, 18, 0, 5, 8, 16, 18, 18, 1, 4, 6, 14, 8, 18, 7, 7, 10]
  n = [0.112225607199012, -0.339005953606712e1, -0.320503911730094e2, -0.197597305104900e3, -0.407693861553446e3, 0.132943775222331e5, 0.170846839774007e1, 0.373694198142245e2, 0.358144365815434e4, 0.423014446424664e6, -0.751071025760063e9, 0.523446127607898e2, -0.228351290812417e3, -0.960652417056937e6, -0.807059292526074e8, 0.162698017225669e13, 0.772465073604171, 0.463929973837746e5, -0.137317885134128e8, 0.170470392630512e13, -0.251104628187308e14, 0.317748830835520e14, 0.538685623675312e2, -0.553089094625169e5, -0.102861522421405e7, 0.204249418756234e13, 0.273918446626977e9, -0.263963146312685e16, -0.107890854108088e10, -0.296492620980124e11, -0.111754907323424e16]

  nu = h/3500
  sigma = s/5.9
  suma = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    ni = n[iter];
    suma += ni * (nu-0.7)**i * (sigma-1.1)**j
  }
  return 100*suma**4
}

void _Backward2_P_hs(h, s)
{
  /*
  """Backward equation for region 2, P=f(h,s)

  Parameters
  ----------
  h : float
    Specific enthalpy, [kJ/kg]
  s : float
    Specific entropy, [kJ/kgK]

  Returns
  -------
  P : float
    Pressure, [MPa]
  """
  */

  sfbc = 5.85
  hamin = _hab_s(s)
  if(h <= hamin)
  {
    P = _Backward2a_P_hs(h, s)
  }
  else if(s >= sfbc)
  {
    P = _Backward2b_P_hs(h, s)
  }
  else
  {
    P = _Backward2c_P_hs(h, s)
  }
  return P
}

// Region 3
void _Region3(rho, T)
{
  /*
  """Basic equation for region 3

  Parameters
  ----------
  rho : float
    Density, [kg/m³]
  T : float
    Temperature, [K]

  Returns
  -------
  prop : dict
    Dict with calculated properties. The available properties are:

      * v: Specific volume, [m³/kg]
      * h: Specific enthalpy, [kJ/kg]
      * s: Specific entropy, [kJ/kgK]
      * cp: Specific isobaric heat capacity, [kJ/kgK]
      * cv: Specific isocoric heat capacity, [kJ/kgK]
      * w: Speed of sound, [m/s]
      * alfav: Cubic expansion coefficient, [1/K]
      * kt: Isothermal compressibility, [1/MPa]

  References
  ----------
  IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  Thermodynamic Properties of Water and Steam August 2007,
  http://www.iapws.org/relguide/IF97-Rev.html, Eq 28

  Examples
  --------
  >>> _Region3(500,650)["P"]
  25.5837018
  >>> _Region3(500,650)["h"]
  1863.43019
  >>> p = _Region3(500, 650)
  >>> p["h"]-p["P"]*1000*p["v"]
  1812.26279
  >>> _Region3(200,650)["s"]
  4.85438792
  >>> _Region3(200,650)["cp"]
  44.6579342
  >>> _Region3(200,650)["cv"]
  4.04118076
  >>> _Region3(200,650)["w"]
  383.444594
  >>> _Region3(500,750)["alfav"]
  0.00441515098
  >>> _Region3(500,750)["kt"]
  0.00806710817
  """
  */

  I = [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 8, 9, 9, 10, 10, 11]
  J = [0, 1, 2, 7, 10, 12, 23, 2, 6, 15, 17, 0, 2, 6, 7, 22, 26, 0, 2, 4, 16, 26, 0, 2, 4, 26, 1, 3, 26, 0, 2, 26, 2, 26, 2, 26, 0, 1, 26]
  n = [-0.15732845290239e2, 0.20944396974307e2, -0.76867707878716e1, 0.26185947787954e1, -0.28080781148620e1, 0.12053369696517e1, -0.84566812812502e-2, -0.12654315477714e1, -0.11524407806681e1, 0.88521043984318, -0.64207765181607, 0.38493460186671, -0.85214708824206, 0.48972281541877e1, -0.30502617256965e1, 0.39420536879154e-1, 0.12558408424308, -0.27999329698710, 0.13899799569460e1, -0.20189915023570e1, -0.82147637173963e-2, -0.47596035734923, 0.43984074473500e-1, -0.44476435428739, 0.90572070719733, .70522450087967, .10770512626332, -0.32913623258954, -0.50871062041158, -0.22175400873096e-1, 0.94260751665092e-1, 0.16436278447961, -0.13503372241348e-1, -0.14834345352472e-1, 0.57922953628084e-3, 0.32308904703711e-2, 0.80964802996215e-4, -0.16557679795037e-3, -0.44923899061815e-4]

  d = rho/rhoc
  Tr = Tc/T
  g = 1.0658070028513*log(d)
  gd = 1.0658070028513/d
  gdd = -1.0658070028513/d**2
  gt = gtt = gdt = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    ni = n[iter];
    g += ni * d**i * Tr**j
    gd += ni*i * d**(i-1) * Tr**j
    gdd += ni*i*(i-1) * d**(i-2) * Tr**j
    gt += ni*j * d**i * Tr**(j-1)
    gtt += ni*j*(j-1) * d**i * Tr**(j-2)
    gdt += ni*i*j * d**(i-1) * Tr**(j-1)
  }
  propiedades = {}
  propiedades["T"] = T
  propiedades["P"] = d*gd*R*T*rho/1000
  propiedades["v"] = 1/rho
  propiedades["h"] = R*T*(Tr*gt+d*gd)
  propiedades["s"] = R*(Tr*gt-g)
  propiedades["cp"] = R*(-Tr**2*gtt+(d*gd-d*Tr*gdt)**2/(2*d*gd+d**2*gdd))
  propiedades["cv"] = -R*Tr**2*gtt
  propiedades["w"] = sqrt(R*T*1000*(2*d*gd+d**2*gdd-(d*gd-d*Tr*gdt)**2 /
                    Tr**2/gtt))
  propiedades["alfav"] = (gd-Tr*gdt)/(2*gd+d*gdd)/T
  propiedades["kt"] = 1/(2*d*gd+d**2*gdd)/rho/R/T*1000
  propiedades["region"] = 3
  propiedades["x"] = 1
  return propiedades
}

void _h_3ab(P)
{
  /*
  """Define the boundary between Region 3a-3b, h=f(P)

  Parameters
  ----------
  P : float
    Pressure, [MPa]

  Returns
  -------
  h : float
    Specific enthalpy, [kJ/kg]

  Examples
  --------
  >>> _h_3ab(25)
  2095.936454
  """
  */
  return 0.201464004206875e4 + 3.74696550136983*P - \
    0.0219921901054187*P**2+0.875131686009950e-4*P**3
}

void _tab_P(P)
{
  /*
  """Define the boundary between Region 3a-3b, T=f(P)

  Parameters
  ----------
  P : float
    Pressure, [MPa]

  Returns
  -------
  T : float
    Temperature, [K]

  References
  ----------
  IAPWS, Revised Supplementary Release on Backward Equations for Specific
  Volume as a Function of Pressure and Temperature v(p,T) for Region 3 of the
  IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water
  and Steam, http://www.iapws.org/relguide/Supp-VPT3-2016.pdf, Eq. 2

  Examples
  --------
  >>> _tab_P(40)
  693.0341408
  """
  */
  I = [0, 1, 2, -1, -2]
  n = [0.154793642129415e4, -0.187661219490113e3, 0.213144632222113e2, -0.191887498864292e4, 0.918419702359447e3]

  Pr = P/1
  T = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    ni= n[iter];
    T += ni * log(Pr)**i
  }
  return T
}

void _top_P(P)
{
  /*
  """Define the boundary between Region 3o-3p, T=f(P)

  Parameters
  ----------
  P : float
    Pressure, [MPa]

  Returns
  -------
  T : float
    Temperature, [K]

  References
  ----------
  IAPWS, Revised Supplementary Release on Backward Equations for Specific
  Volume as a Function of Pressure and Temperature v(p,T) for Region 3 of the
  IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water
  and Steam, http://www.iapws.org/relguide/Supp-VPT3-2016.pdf, Eq. 2

  Examples
  --------
  >>> _top_P(22.8)
  650.0106943
  """
  */
  I = [0, 1, 2, -1, -2]
  n = [0.969461372400213e3, -0.332500170441278e3, 0.642859598466067e2, 0.773845935768222e3, -0.152313732937084e4]

  Pr = P/1
  T = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    ni= n[iter];
    T += ni * log(Pr)**i
  }
  return T
}

void _twx_P(P)
{
  /*
  """Define the boundary between Region 3w-3x, T=f(P)

  Parameters
  ----------
  P : float
    Pressure, [MPa]

  Returns
  -------
  T : float
    Temperature, [K]

  References
  ----------
  IAPWS, Revised Supplementary Release on Backward Equations for Specific
  Volume as a Function of Pressure and Temperature v(p,T) for Region 3 of the
  IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water
  and Steam, http://www.iapws.org/relguide/Supp-VPT3-2016.pdf, Eq. 2

  Examples
  --------
  >>> _twx_P(22.3)
  648.2049480
  """
  */
  I = [0, 1, 2, -1, -2]
  n = [0.728052609145380e1, 0.973505869861952e2, 0.147370491183191e2, 0.329196213998375e3, 0.873371668682417e3]

  Pr = P/1
  T = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    ni= n[iter];
    T += ni * log(Pr)**i
  }
  return T
}

void _tef_P(P)
{
  /*
  """Define the boundary between Region 3e-3f, T=f(P)

  Parameters
  ----------
  P : float
    Pressure, [MPa]

  Returns
  -------
  T : float
    Temperature, [K]

  References
  ----------
  IAPWS, Revised Supplementary Release on Backward Equations for Specific
  Volume as a Function of Pressure and Temperature v(p,T) for Region 3 of the
  IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water
  and Steam, http://www.iapws.org/relguide/Supp-VPT3-2016.pdf, Eq. 3

  Examples
  --------
  >>> _tef_P(40)
  713.9593992
  """
  */
  return 3.727888004*(P-22.064)+647.096
}

void _txx_P(P, xy)
{
  /*
  """Define the boundary between 3x-3y, T=f(P)

  Parameters
  ----------
  P : float
    Pressure, [MPa]
  xy: string
    Subregions options: cd, gh, ij, jk, mn, qu, rx, uv

  Returns
  -------
  T : float
    Temperature, [K]

  References
  ----------
  IAPWS, Revised Supplementary Release on Backward Equations for Specific
  Volume as a Function of Pressure and Temperature v(p,T) for Region 3 of the
  IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water
  and Steam, http://www.iapws.org/relguide/Supp-VPT3-2016.pdf, Eq. 1

  Examples
  --------
  >>> _txx_P(25,"cd")
  649.3659208
  >>> _txx_P(23,"gh")
  649.8873759
  >>> _txx_P(23,"ij")
  651.5778091
  >>> _txx_P(23,"jk")
  655.8338344
  >>> _txx_P(22.8,"mn")
  649.6054133
  >>> _txx_P(22,"qu")
  645.6355027
  >>> _txx_P(22,"rx")
  648.2622754
  >>> _txx_P(22.3,"uv")
  647.7996121
  """
  */

  ng =
  {
    "cd": [0.585276966696349e3, 0.278233532206915e1, -0.127283549295878e-1, 0.159090746562729e-3],
    "gh": [-0.249284240900418e5, 0.428143584791546e4, -0.269029173140130e3, 0.751608051114157e1, -0.787105249910383e-1],
    "ij": [0.584814781649163e3, -0.616179320924617, 0.260763050899562, -0.587071076864459e-2, 0.515308185433082e-4],
    "jk": [0.617229772068439e3, -0.770600270141675e1, 0.697072596851896, -0.157391839848015e-1, 0.137897492684194e-3],
    "mn": [0.535339483742384e3, 0.761978122720128e1, -0.158365725441648, 0.192871054508108e-2],
    "qu": [0.565603648239126e3, 0.529062258221222e1, -0.102020639611016, 0.122240301070145e-2],
    "rx": [0.584561202520006e3, -0.102961025163669e1, 0.243293362700452, -0.294905044740799e-2],
    "uv": [0.528199646263062e3, 0.890579602135307e1, -0.222814134903755, 0.286791682263697e-2]
  }

  n = ng[xy]
  Pr = P/1
  T = 0
  for(int i = 0; i < n.length; i++)
  {
    ni = n[i];
    T += ni * Pr**i;
  }
  return T
}

void _Backward3a_v_Ph(P, h)
{
  /*
  """Backward equation for region 3a, v=f(P,h)

  Parameters
  ----------
  P : float
    Pressure, [MPa]
  h : float
    Specific enthalpy, [kJ/kg]

  References
  ----------
  IAPWS, Revised Supplementary Release on Backward Equations for the
  Functions T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS
  Industrial Formulation 1997 for the Thermodynamic Properties of Water and
  Steam, http://www.iapws.org/relguide/Supp-Tv%28ph,ps%293-2014.pdf, Eq 4

  Returns
  -------
  v : float
    Specific volume, [m³/kg]

  Examples
  --------
  >>> _Backward3a_v_Ph(20,1700)
  0.001749903962
  >>> _Backward3a_v_Ph(100,2100)
  0.001676229776
  """
  */
  I = [-12, -12, -12, -12, -10, -10, -10, -8, -8, -6, -6, -6, -4, -4, -3, -2, -2, -1, -1, -1, -1, 0, 0, 1, 1, 1, 2, 2, 3, 4, 5, 8]
  J = [6, 8, 12, 18, 4, 7, 10, 5, 12, 3, 4, 22, 2, 3, 7, 3, 16, 0, 1, 2, 3, 0, 1, 0, 1, 2, 0, 2, 0, 2, 2, 2]
  n = [0.529944062966028e-2, -0.170099690234461, 0.111323814312927e2, -0.217898123145125e4, -0.506061827980875e-3, 0.556495239685324, -0.943672726094016e1, -0.297856807561527, 0.939353943717186e2, 0.192944939465981e-1, 0.421740664704763, -0.368914126282330e7, -0.737566847600639e-2, -0.354753242424366, -0.199768169338727e1, 0.115456297059049e1, 0.568366875815960e4, 0.808169540124668e-2, 0.172416341519307, 0.104270175292927e1, -0.297691372792847, 0.560394465163593, 0.275234661176914, -0.148347894866012, -0.651142513478515e-1, -0.292468715386302e1, 0.664876096952665e-1, 0.352335014263844e1, -0.146340792313332e-1, -0.224503486668184e1, 0.110533464706142e1, -0.408757344495612e-1]

  Pr = P/100
  nu = h/2100
  suma = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    ni = n[iter];
    suma += ni * (Pr+0.128)**i * (nu-0.727)**j
  }
  return 0.0028*suma
}

void _Backward3b_v_Ph(P, h)
{
  /*
  """Backward equation for region 3b, v=f(P,h)

  Parameters
  ----------
  P : float
    Pressure, [MPa]
  h : float
    Specific enthalpy, [kJ/kg]

  Returns
  -------
  v : float
    Specific volume, [m³/kg]

  References
  ----------
  IAPWS, Revised Supplementary Release on Backward Equations for the
  Functions T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS
  Industrial Formulation 1997 for the Thermodynamic Properties of Water and
  Steam, http://www.iapws.org/relguide/Supp-Tv%28ph,ps%293-2014.pdf, Eq 5

  Examples
  --------
  >>> _Backward3b_v_Ph(20,2500)
  0.006670547043
  >>> _Backward3b_v_Ph(100,2700)
  0.002404234998
  """
  */
  I = [-12, -12, -8, -8, -8, -8, -8, -8, -6, -6, -6, -6, -6, -6, -4, -4, -4, -3, -3, -2, -2, -1, -1, -1, -1, 0, 1, 1, 2, 2]
  J = [0, 1, 0, 1, 3, 6, 7, 8, 0, 1, 2, 5, 6, 10, 3, 6, 10, 0, 2, 1, 2, 0, 1, 4, 5, 0, 0, 1, 2, 6]
  n = [-0.225196934336318e-8, 0.140674363313486e-7, 0.233784085280560e-5, -0.331833715229001e-4, 0.107956778514318e-2, -0.271382067378863, 0.107202262490333e1, -0.853821329075382, -0.215214194340526e-4, 0.769656088222730e-3, -0.431136580433864e-2, 0.453342167309331, -0.507749535873652, -0.100475154528389e3, -0.219201924648793, -0.321087965668917e1, 0.607567815637771e3, 0.557686450685932e-3, 0.187499040029550, 0.905368030448107e-2, 0.285417173048685, 0.329924030996098e-1, 0.239897419685483, 0.482754995951394e1, -0.118035753702231e2, 0.169490044091791, -0.179967222507787e-1, 0.371810116332674e-1, -0.536288335065096e-1, 0.160697101092520e1]

  Pr = P/100
  nu = h/2800
  suma = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    ni = n[iter];
    suma += ni * (Pr+0.0661)**i * (nu-0.72)**j
  }
  return 0.0088*suma
}

void _Backward3_v_Ph(P, h)
{
  /*
  """Backward equation for region 3, v=f(P,h)

  Parameters
  ----------
  P : float
    Pressure, [MPa]
  h : float
    Specific enthalpy, [kJ/kg]

  Returns
  -------
  v : float
    Specific volume, [m³/kg]
  """
  */

  hf = _h_3ab(P)
  if(h <= hf)
  {
    return _Backward3a_v_Ph(P, h)
  }
  else
  {
    return _Backward3b_v_Ph(P, h)
  }
}

void _Backward3a_T_Ph(P, h)
{
  /*
  """Backward equation for region 3a, T=f(P,h)

  Parameters
  ----------
  P : float
    Pressure, [MPa]
  h : float
    Specific enthalpy, [kJ/kg]

  Returns
  -------
  T : float
    Temperature, [K]

  References
  ----------
  IAPWS, Revised Supplementary Release on Backward Equations for the
  Functions T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS
  Industrial Formulation 1997 for the Thermodynamic Properties of Water and
  Steam, http://www.iapws.org/relguide/Supp-Tv%28ph,ps%293-2014.pdf, Eq 2

  Examples
  --------
  >>> _Backward3a_T_Ph(20,1700)
  629.3083892
  >>> _Backward3a_T_Ph(100,2100)
  733.6163014
  """
  */
  I = [-12, -12, -12, -12, -12, -12, -12, -12, -10, -10, -10, -8, -8, -8, -8, -5, -3, -2, -2, -2, -1, -1, 0, 0, 1, 3, 3, 4, 4, 10, 12]
  J = [0, 1, 2, 6, 14, 16, 20, 22, 1, 5, 12, 0, 2, 4, 10, 2, 0, 1, 3, 4, 0, 2, 0, 1, 1, 0, 1, 0, 3, 4, 5]
  n = [-0.133645667811215e-6, 0.455912656802978e-5, -0.146294640700979e-4, 0.639341312970080e-2, 0.372783927268847e3, -0.718654377460447e4, 0.573494752103400e6, -0.267569329111439e7, -0.334066283302614e-4, -0.245479214069597e-1, 0.478087847764996e2, 0.764664131818904e-5, 0.128350627676972e-2, 0.171219081377331e-1, -0.851007304583213e1, -0.136513461629781e-1, -0.384460997596657e-5, 0.337423807911655e-2, -0.551624873066791, 0.729202277107470, -0.992522757376041e-2, -.119308831407288, .793929190615421, .454270731799386, .20999859125991, -0.642109823904738e-2, -0.235155868604540e-1, 0.252233108341612e-2, -0.764885133368119e-2, 0.136176427574291e-1, -0.133027883575669e-1]

  Pr = P/100.
  nu = h/2300.
  suma = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    n = n[iter];
    suma += n*(Pr+0.240)**i*(nu-0.615)**j
  }
  return 760*suma
}

void _Backward3b_T_Ph(P, h)
{
  /*
  """Backward equation for region 3b, T=f(P,h)

  Parameters
  ----------
  P : float
    Pressure, [MPa]
  h : float
    Specific enthalpy, [kJ/kg]

  Returns
  -------
  T : float
    Temperature, [K]

  References
  ----------
  IAPWS, Revised Supplementary Release on Backward Equations for the
  Functions T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS
  Industrial Formulation 1997 for the Thermodynamic Properties of Water and
  Steam, http://www.iapws.org/relguide/Supp-Tv%28ph,ps%293-2014.pdf, Eq 3

  Examples
  --------
  >>> _Backward3b_T_Ph(20,2500)
  641.8418053
  >>> _Backward3b_T_Ph(100,2700)
  842.0460876
  """
  */
  I = [-12, -12, -10, -10, -10, -10, -10, -8, -8, -8, -8, -8, -6, -6, -6, -4, -4, -3, -2, -2, -1, -1, -1, -1, -1, -1, 0, 0, 1, 3, 5, 6, 8]
  J = [0, 1, 0, 1, 5, 10, 12, 0, 1, 2, 4, 10, 0, 1, 2, 0, 1, 5, 0, 4, 2, 4, 6, 10, 14, 16, 0, 2, 1, 1, 1, 1, 1]
  n = [0.323254573644920e-4, -0.127575556587181e-3, -0.475851877356068e-3, 0.156183014181602e-2, 0.105724860113781, -0.858514221132534e2, 0.724140095480911e3, 0.296475810273257e-2, -0.592721983365988e-2, -0.126305422818666e-1, -0.115716196364853, 0.849000969739595e2, -0.108602260086615e-1, 0.154304475328851e-1, 0.750455441524466e-1, 0.252520973612982e-1, -0.602507901232996e-1, -0.307622221350501e1, -0.574011959864879e-1, 0.503471360939849e1, -0.925081888584834, 0.391733882917546e1, -0.773146007130190e2, 0.949308762098587e4, -0.141043719679409e7, 0.849166230819026e7, 0.861095729446704, 0.323346442811720, 0.873281936020439, -0.436653048526683, 0.286596714529479, -0.131778331276228, 0.676682064330275e-2]

  Pr = P/100.
  nu = h/2800.
  suma = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    n = n[iter];
    suma += n*(Pr+0.298)**i*(nu-0.72)**j
  }
  return 860*suma
}

void _Backward3_T_Ph(P, h)
{
  /*
  """Backward equation for region 3, T=f(P,h)

  Parameters
  ----------
  P : float
    Pressure, [MPa]
  h : float
    Specific enthalpy, [kJ/kg]

  Returns
  -------
  T : float
    Temperature, [K]
  """
  */
  hf = _h_3ab(P)
  if(h <= hf)
  {
    T = _Backward3a_T_Ph(P, h)
  }
  else
  {
    T = _Backward3b_T_Ph(P, h)
  }
  return T
}

void _Backward3a_v_Ps(P, s)
{
  /*
  """Backward equation for region 3a, v=f(P,s)

  Parameters
  ----------
  P : float
    Pressure, [MPa]
  s : float
    Specific entropy, [kJ/kgK]

  Returns
  -------
  v : float
    Specific volume, [m³/kg]

  References
  ----------
  IAPWS, Revised Supplementary Release on Backward Equations for the
  Functions T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS
  Industrial Formulation 1997 for the Thermodynamic Properties of Water and
  Steam, http://www.iapws.org/relguide/Supp-Tv%28ph,ps%293-2014.pdf, Eq 8

  Examples
  --------
  >>> _Backward3a_v_Ps(20,3.8)
  0.001733791463
  >>> _Backward3a_v_Ps(100,4)
  0.001555893131
  """
  */
  I = [-12, -12, -12, -10, -10, -10, -10, -8, -8, -8, -8, -6, -5, -4, -3, -3, -2, -2, -1, -1, 0, 0, 0, 1, 2, 4, 5, 6]
  J = [10, 12, 14, 4, 8, 10, 20, 5, 6, 14, 16, 28, 1, 5, 2, 4, 3, 8, 1, 2, 0, 1, 3, 0, 0, 2, 2, 0]
  n = [0.795544074093975e2, -0.238261242984590e4, 0.176813100617787e5, -0.110524727080379e-2, -0.153213833655326e2, 0.297544599376982e3, -0.350315206871242e8, 0.277513761062119, -0.523964271036888, -0.148011182995403e6, 0.160014899374266e7, 0.170802322663427e13, 0.246866996006494e-3, 0.165326084797980e1, -0.118008384666987, 0.253798642355900e1, 0.965127704669424, -0.282172420532826e2, 0.203224612353823, 0.110648186063513e1, 0.526127948451280, 0.277000018736321, 0.108153340501132e1, -0.744127885357893e-1, 0.164094443541384e-1, -0.680468275301065e-1, 0.257988576101640e-1, -0.145749861944416e-3]

  Pr = P/100
  sigma = s/4.4
  suma = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    ni = n[iter];
    suma += ni * (Pr+0.187)**i * (sigma-0.755)**j
  }
  return 0.0028*suma
}

void _Backward3b_v_Ps(P, s)
{
  /*
  """Backward equation for region 3b, v=f(P,s)

  Parameters
  ----------
  P : float
    Pressure, [MPa]
  s : float
    Specific entropy, [kJ/kgK]

  Returns
  -------
  v : float
    Specific volume, [m³/kg]

  References
  ----------
  IAPWS, Revised Supplementary Release on Backward Equations for the
  Functions T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS
  Industrial Formulation 1997 for the Thermodynamic Properties of Water and
  Steam, http://www.iapws.org/relguide/Supp-Tv%28ph,ps%293-2014.pdf, Eq 9

  Examples
  --------
  >>> _Backward3b_v_Ps(20,5)
  0.006262101987
  >>> _Backward3b_v_Ps(100,5)
  0.002449610757
  """
  */
  I = [-12, -12, -12, -12, -12, -12, -10, -10, -10, -10, -8, -5, -5, -5, -4, -4, -4, -4, -3, -2, -2, -2, -2, -2, -2, 0, 0, 0, 1, 1, 2]
  J = [0, 1, 2, 3, 5, 6, 0, 1, 2, 4, 0, 1, 2, 3, 0, 1, 2, 3, 1, 0, 1, 2, 3, 4, 12, 0, 1, 2, 0, 2, 2]
  n = [0.591599780322238e-4, -0.185465997137856e-2, 0.104190510480013e-1, 0.598647302038590e-2, -0.771391189901699, 0.172549765557036e1, -0.467076079846526e-3, 0.134533823384439e-1, -0.808094336805495e-1, 0.508139374365767, 0.128584643361683e-2, -0.163899353915435e1, 0.586938199318063e1, -0.292466667918613e1, -0.614076301499537e-2, 0.576199014049172e1, -0.121613320606788e2, 0.167637540957944e1, -0.744135838773463e1, 0.378168091437659e-1, 0.401432203027688e1, 0.160279837479185e2, 0.317848779347728e1, -0.358362310304853e1, -0.115995260446827e7, 0.199256573577909, -0.122270624794624, -0.191449143716586e2, -0.150448002905284e-1, 0.146407900162154e2, -0.327477787188230e1]

  Pr = P/100
  sigma = s/5.3
  suma = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    ni = n[iter];
    suma += ni * (Pr+0.298)**i * (sigma-0.816)**j
  }
  return 0.0088*suma
}

void _Backward3_v_Ps(P, s)
{
  /*
  """Backward equation for region 3, v=f(P,s)

  Parameters
  ----------
  P : float
    Pressure, [MPa]
  s : float
    Specific entropy, [kJ/kgK]

  Returns
  -------
  v : float
    Specific volume, [m³/kg]
  """
  */
  if(s <= sc)
  {
    return _Backward3a_v_Ps(P, s)
  }
  else
  {
    return _Backward3b_v_Ps(P, s)
  }
}

void _Backward3a_T_Ps(P, s)
{
  /*
  """Backward equation for region 3a, T=f(P,s)

  Parameters
  ----------
  P : float
    Pressure, [MPa]
  s : float
    Specific entropy, [kJ/kgK]

  Returns
  -------
  T : float
    Temperature, [K]

  References
  ----------
  IAPWS, Revised Supplementary Release on Backward Equations for the
  Functions T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS
  Industrial Formulation 1997 for the Thermodynamic Properties of Water and
  Steam, http://www.iapws.org/relguide/Supp-Tv%28ph,ps%293-2014.pdf, Eq 6

  Examples
  --------
  >>> _Backward3a_T_Ps(20,3.8)
  628.2959869
  >>> _Backward3a_T_Ps(100,4)
  705.6880237
  """
  */
  I = [-12, -12, -10, -10, -10, -10, -8, -8, -8, -8, -6, -6, -6, -5, -5, -5, -4, -4, -4, -2, -2, -1, -1, 0, 0, 0, 1, 2, 2, 3, 8, 8, 10]
  J = [28, 32, 4, 10, 12, 14, 5, 7, 8, 28, 2, 6, 32, 0, 14, 32, 6, 10, 36, 1, 4, 1, 6, 0, 1, 4, 0, 0, 3, 2, 0, 1, 2]
  n = [0.150042008263875e10, -0.159397258480424e12, 0.502181140217975e-3, -0.672057767855466e2, 0.145058545404456e4, -0.823889534888890e4, -0.154852214233853, 0.112305046746695e2, -0.297000213482822e2, 0.438565132635495e11, 0.137837838635464e-2, -0.297478527157462e1, 0.971777947349413e13, -0.571527767052398e-4, 0.288307949778420e5, -0.744428289262703e14, 0.128017324848921e2, -0.368275545889071e3, 0.664768904779177e16, 0.449359251958880e-1, -0.422897836099655e1, -0.240614376434179, -0.474341365254924e1, 0.724093999126110, 0.923874349695897, 0.399043655281015e1, 0.384066651868009e-1, -0.359344365571848e-2, -0.735196448821653, 0.188367048396131, 0.141064266818704e-3, -0.257418501496337e-2, 0.123220024851555e-2]

  Pr = P/100
  sigma = s/4.4
  suma = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    ni = n[iter];
    suma += ni * (Pr+0.240)**i * (sigma-0.703)**j
  }
  return 760*suma
}

void _Backward3b_T_Ps(P, s)
{
  /*
  """Backward equation for region 3b, T=f(P,s)

  Parameters
  ----------
  P : float
    Pressure, [MPa]
  s : float
    Specific entropy, [kJ/kgK]

  Returns
  -------
  T : float
    Temperature, [K]

  References
  ----------
  IAPWS, Revised Supplementary Release on Backward Equations for the
  Functions T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS
  Industrial Formulation 1997 for the Thermodynamic Properties of Water and
  Steam, http://www.iapws.org/relguide/Supp-Tv%28ph,ps%293-2014.pdf, Eq 7

  Examples
  --------
  >>> _Backward3b_T_Ps(20,5)
  640.1176443
  >>> _Backward3b_T_Ps(100,5)
  847.4332825
  """
  */
  I = [-12, -12, -12, -12, -8, -8, -8, -6, -6, -6, -5, -5, -5, -5, -5, -4, -3, -3, -2, 0, 2, 3, 4, 5, 6, 8, 12, 14]
  J = [1, 3, 4, 7, 0, 1, 3, 0, 2, 4, 0, 1, 2, 4, 6, 12, 1, 6, 2, 0, 1, 1, 0, 24, 0, 3, 1, 2]
  n = [0.527111701601660, -0.401317830052742e2, 0.153020073134484e3, -0.224799398218827e4, -0.193993484669048, -0.140467557893768e1, 0.426799878114024e2, 0.752810643416743, 0.226657238616417e2, -0.622873556909932e3, -0.660823667935396, 0.841267087271658, -0.253717501764397e2, 0.485708963532948e3, 0.880531517490555e3, 0.265015592794626e7, -0.359287150025783, -0.656991567673753e3, 0.241768149185367e1, 0.856873461222588, 0.655143675313458, -0.213535213206406, 0.562974957606348e-2, -0.316955725450471e15, -0.699997000152457e-3, 0.119845803210767e-1, 0.193848122022095e-4, -0.215095749182309e-4]

  Pr = P/100
  sigma = s/5.3
  suma = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    ni = n[iter];
    suma += ni * (Pr+0.760)**i * (sigma-0.818)**j
  }
  return 860*suma
}

void _Backward3_T_Ps(P, s)
{
  /*
  """Backward equation for region 3, T=f(P,s)

  Parameters
  ----------
  P : float
    Pressure, [MPa]
  s : float
    Specific entropy, [kJ/kgK]

  Returns
  -------
  T : float
    Temperature, [K]
  """
  */
  sc = 4.41202148223476
  if(s <= sc)
  {
    T = _Backward3a_T_Ps(P, s)
  }
  else
  {
    T = _Backward3b_T_Ps(P, s)
  }
  return T
}

void _Backward3a_P_hs(h, s)
{
  /*
  """Backward equation for region 3a, P=f(h,s)

  Parameters
  ----------
  h : float
    Specific enthalpy, [kJ/kg]
  s : float
    Specific entropy, [kJ/kgK]

  Returns
  -------
  P : float
    Pressure, [MPa]

  References
  ----------
  IAPWS, Revised Supplementary Release on Backward Equations p(h,s) for
  Region 3, Equations as a Function of h and s for the Region Boundaries, and
  an Equation Tsat(h,s) for Region 4 of the IAPWS Industrial Formulation 1997
  for the Thermodynamic Properties of Water and Steam,
  http://www.iapws.org/relguide/Supp-phs3-2014.pdf. Eq 1

  Examples
  --------
  >>> _Backward3a_P_hs(1700,3.8)
  25.55703246
  >>> _Backward3a_P_hs(2000,4.2)
  45.40873468
  >>> _Backward3a_P_hs(2100,4.3)
  60.78123340
  """
  */
  I = [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 6, 7, 8, 10, 10, 14, 18, 20, 22, 22, 24, 28, 28, 32, 32]
  J = [0, 1, 5, 0, 3, 4, 8, 14, 6, 16, 0, 2, 3, 0, 1, 4, 5, 28, 28, 24, 1, 32, 36, 22, 28, 36, 16, 28, 36, 16, 36, 10, 28]
  n = [0.770889828326934e1, -0.260835009128688e2, 0.267416218930389e3, 0.172221089496844e2, -0.293542332145970e3, 0.614135601882478e3, -0.610562757725674e5, -0.651272251118219e8, 0.735919313521937e5, -0.116646505914191e11, 0.355267086434461e2, -0.596144543825955e3, -0.475842430145708e3, 0.696781965359503e2, 0.335674250377312e3, 0.250526809130882e5, 0.146997380630766e6, 0.538069315091534e20, 0.143619827291346e22, 0.364985866165994e20, -0.254741561156775e4, 0.240120197096563e28, -0.393847464679496e30, 0.147073407024852e25, -0.426391250432059e32, 0.194509340621077e39, 0.666212132114896e24, 0.706777016552858e34, 0.175563621975576e42, 0.108408607429124e29, 0.730872705175151e44, 0.159145847398870e25, 0.377121605943324e41]

  nu = h/2300
  sigma = s/4.4
  suma = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    ni = n[iter];
    suma += ni * (nu-1.01)**i * (sigma-0.75)**j
  }
  return 99*suma
}

void _Backward3b_P_hs(h, s)
{
  /*
  """Backward equation for region 3b, P=f(h,s)

  Parameters
  ----------
  h : float
    Specific enthalpy, [kJ/kg]
  s : float
    Specific entropy, [kJ/kgK]

  Returns
  -------
  P : float
    Pressure, [MPa]

  References
  ----------
  IAPWS, Revised Supplementary Release on Backward Equations p(h,s) for
  Region 3, Equations as a Function of h and s for the Region Boundaries, and
  an Equation Tsat(h,s) for Region 4 of the IAPWS Industrial Formulation 1997
  for the Thermodynamic Properties of Water and Steam,
  http://www.iapws.org/relguide/Supp-phs3-2014.pdf. Eq 1

  Examples
  --------
  >>> _Backward3b_P_hs(2400,4.7)
  63.63924887
  >>> _Backward3b_P_hs(2600,5.1)
  34.34999263
  >>> _Backward3b_P_hs(2700,5.0)
  88.39043281
  """
  */
  I = [-12, -12, -12, -12, -12, -10, -10, -10, -10, -8, -8, -6, -6, -6, -6, -5, -4, -4, -4, -3, -3, -3, -3, -2, -2, -1, 0, 2, 2, 5, 6, 8, 10, 14, 14]
  J = [2, 10, 12, 14, 20, 2, 10, 14, 18, 2, 8, 2, 6, 7, 8, 10, 4, 5, 8, 1, 3, 5, 6, 0, 1, 0, 3, 0, 1, 0, 1, 1, 1, 3, 7]
  n = [0.125244360717979e-12, -0.126599322553713e-1, 0.506878030140626e1, 0.317847171154202e2, -0.391041161399932e6, -0.975733406392044e-10, -0.186312419488279e2, 0.510973543414101e3, 0.373847005822362e6, 0.299804024666572e-7, 0.200544393820342e2, -0.498030487662829e-5, -0.102301806360030e2, 0.552819126990325e2, -0.206211367510878e3, -0.794012232324823e4, 0.782248472028153e1, -0.586544326902468e2, 0.355073647696481e4, -0.115303107290162e-3, -0.175092403171802e1, 0.257981687748160e3, -0.727048374179467e3, 0.121644822609198e-3, 0.393137871762692e-1, 0.704181005909296e-2, -0.829108200698110e2, -0.265178818131250, 0.137531682453991e2, -0.522394090753046e2, 0.240556298941048e4, -0.227361631268929e5, 0.890746343932567e5, -0.239234565822486e8, 0.568795808129714e10]

  nu = h/2800
  sigma = s/5.3
  suma = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    ni = n[iter];
    suma += ni * (nu-0.681)**i * (sigma-0.792)**j
  }
  return 16.6/suma
}

void _Backward3_P_hs(h, s)
{
  /*
  """Backward equation for region 3, P=f(h,s)

  Parameters
  ----------
  h : float
    Specific enthalpy, [kJ/kg]
  s : float
    Specific entropy, [kJ/kgK]

  Returns
  -------
  P : float
    Pressure, [MPa]
  """
  */
  sc = 4.41202148223476
  if(s <= sc)
  {
    return _Backward3a_P_hs(h, s)
  }
  else
  {
    return _Backward3b_P_hs(h, s)
  }
}

void _Backward3_sat_v_P(P, T, x)
{
  /*
  """Backward equation for region 3 for saturated state, vs=f(P,x)

  Parameters
  ----------
  T : float
    Temperature, [K]
  P : float
    Pressure, [MPa]
  x : integer
    Vapor quality, [-]

  Returns
  -------
  v : float
    Specific volume, [m³/kg]

  Notes
  -----
  The vapor quality (x) can be 0 (saturated liquid) or 1 (saturated vapour)
  """
  */
  if(x == 0)
  {
    if(P < 19.00881189)
    {
      region = "c"
    }
    else if(P < 21.0434)
    {
      region = "s"
    }
    else if(P < 21.9316)
    {
      region = "u"
    }
    else
    {
      region = "y"
    }
  }
  else
  {
    if(P < 20.5)
    {
      region = "t"
    }
    else if(P < 21.0434)
    {
      region = "r"
    }
    else if(P < 21.9009)
    {
      region = "x"
    }
    else
    {
      region = "z"
    }
  }
  return _Backward3x_v_PT(T, P, region)
}

void _Backward3_v_PT(P, T)
{
  /*
  """Backward equation for region 3, v=f(P,T)

  Parameters
  ----------
  T : float
    Temperature, [K]
  P : float
    Pressure, [MPa]

  Returns
  -------
  v : float
    Specific volume, [m³/kg]

  References
  ----------
  IAPWS, Revised Supplementary Release on Backward Equations for Specific
  Volume as a Function of Pressure and Temperature v(p,T) for Region 3 of the
  IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water
  and Steam, http://www.iapws.org/relguide/Supp-VPT3-2016.pdf, Table 2 and 10
  """
  */
  if(P > 40)
  {
    if(T <= _tab_P(P))
    {
      region = "a"
    }
    else
    {
      region = "b"
    }
  }
  else if(25 < P <= 40)
  {
    tcd = _txx_P(P, "cd")
    tab = _tab_P(P)
    tef = _tef_P(P)
    if(T <= tcd)
    {
      region = "c"
    }
    else if(tcd < T <= tab)
    {
      region = "d"
    }
    else if(tab < T <= tef)
    {
      region = "e"
    }
    else
    {
      region = "f"
    }
  }
  else if(23.5 < P <= 25)
  {
    tcd = _txx_P(P, "cd")
    tgh = _txx_P(P, "gh")
    tef = _tef_P(P)
    tij = _txx_P(P, "ij")
    tjk = _txx_P(P, "jk")
    if(T <= tcd)
    {
      region = "c"
    }
    else if(tcd < T <= tgh)
    {
      region = "g"
    }
    else if(tgh < T <= tef)
    {
      region = "h"
    }
    else if(tef < T <= tij)
    {
      region = "i"
    }
    else if(tij < T <= tjk)
    {
      region = "j"
    }
    else
    {
      region = "k"
    }
  }
  else if(23 < P <= 23.5)
  {
    tcd = _txx_P(P, "cd")
    tgh = _txx_P(P, "gh")
    tef = _tef_P(P)
    tij = _txx_P(P, "ij")
    tjk = _txx_P(P, "jk")
    if(T <= tcd)
    {
      region = "c"
    }
    else if(tcd < T <= tgh)
    {
      region = "l"
    }
    else if(tgh < T <= tef)
    {
      region = "h"
    }
    else if(tef < T <= tij)
    {
      region = "i"
    }
    else if(tij < T <= tjk)
    {
      region = "j"
    }
    else
    {
      region = "k"
    }
  }
  else if(22.5 < P <= 23)
  {
    tcd = _txx_P(P, "cd")
    tgh = _txx_P(P, "gh")
    tmn = _txx_P(P, "mn")
    tef = _tef_P(P)
    top = _top_P(P)
    tij = _txx_P(P, "ij")
    tjk = _txx_P(P, "jk")
    if(T <= tcd)
    {
      region = "c"
    }
    else if(tcd < T <= tgh)
    {
      region = "l"
    }
    else if(tgh < T <= tmn)
    {
      region = "m"
    }
    else if(tmn < T <= tef)
    {
      region = "n"
    }
    else if(tef < T <= top)
    {
      region = "o"
    }
    else if(top < T <= tij)
    {
      region = "p"
    }
    else if(tij < T <= tjk)
    {
      region = "j"
    }
    else
    {
      region = "k"
    }
  }
  else if(_PSat_T(643.15) < P <= 22.5)
  {
    tcd = _txx_P(P, "cd")
    tqu = _txx_P(P, "qu")
    trx = _txx_P(P, "rx")
    tjk = _txx_P(P, "jk")
    if(T <= tcd)
    {
      region = "c"
    }
    else if(tcd < T <= tqu)
    {
      region = "q"
    }
    else if(tqu < T <= trx)
    {
      // Table 10
      tef = _tef_P(P)
      twx = _twx_P(P)
      tuv = _txx_P(P, "uv")
      if(22.11 < P <= 22.5)
      {
        if(T <= tuv)
        {
          region = "u"
        }
        else if(tuv <= T <= tef)
        {
          region = "v"
        }
        else if(tef <= T <= twx)
        {
          region = "w"
        }
        else
        {
          region = "x"
        }
      }
      else if(22.064 < P <= 22.11)
      {
        if(T <= tuv)
        {
          region = "u"
        }
        else if(tuv <= T <= tef)
        {
          region = "y"
        }
        else if(tef <= T <= twx)
        {
          region = "z"
        }
        else
        {
          region = "x"
        }
      }
      else if(T > _TSat_P(P))
      {
        if(_PSat_T(643.15) < P <= 21.90096265)
        {
          region = "x"
        }
        else if(21.90096265 < P <= 22.064)
        {
          if(T <= twx)
          {
            region = "z"
          }
          else
          {
            region = "x"
          }
      }
      else if(T <= _TSat_P(P))
      {
        if(_PSat_T(643.15) < P <= 21.93161551)
        {
          region = "u"
        }
        else if(21.93161551 < P <= 22.064)
        {
          if(T <= tuv)
          {
            region = "u"
          }
          else
          {
            region = "y"
          }
        }
      }
    }
    else if(trx < T <= tjk)
    {
      region = "r"
    }
    else
    {
      region = "k"
    }
  }
  else if(20.5 < P <= _PSat_T(643.15))
  {
    tcd = _txx_P(P, "cd")
    Ts = _TSat_P(P)
    tjk = _txx_P(P, "jk")
    if(T <= tcd)
    {
      region = "c"
    }
    else if(tcd < T <= Ts)
    {
      region = "s"
    }
    else if(Ts < T <= tjk)
    {
      region = "r"
    }
    else
    {
      region = "k"
    }
  }
  else if(19.00881189173929 < P <= 20.5)
  {
    tcd = _txx_P(P, "cd")
    Ts = _TSat_P(P)
    if(T <= tcd)
    {
      region = "c"
    }
    else if(tcd < T <= Ts)
    {
      region = "s"
    }
    else
    {
      region = "t"
    }
  }
  else if(Ps_623 < P <= 19.00881189173929)
  {
    Ts = _TSat_P(P)
    if(T <= Ts)
    {
      region = "c"
    }
    else
    {
      region = "t"
    }
  }
  return _Backward3x_v_PT(T, P, region)
}

void _Backward3x_v_PT(T, P, x)
{
  /*
  """Backward equation for region 3x, v=f(P,T)

  Parameters
  ----------
  T : float
    Temperature, [K]
  P : float
    Pressure, [MPa]
  x : char
    Region 3 subregion code

  Returns
  -------
  v : float
    Specific volume, [m³/kg]

  References
  ----------
  IAPWS, Revised Supplementary Release on Backward Equations for Specific
  Volume as a Function of Pressure and Temperature v(p,T) for Region 3 of the
  IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water
  and Steam, http://www.iapws.org/relguide/Supp-VPT3-2016.pdf, Eq. 4-5

  Examples
  --------
  >>> _Backward3x_v_PT(630,50,"a")
  0.001470853100
  >>> _Backward3x_v_PT(670,80,"a")
  0.001503831359
  >>> _Backward3x_v_PT(710,50,"b")
  0.002204728587
  >>> _Backward3x_v_PT(750,80,"b")
  0.001973692940
  >>> _Backward3x_v_PT(630,20,"c")
  0.001761696406
  >>> _Backward3x_v_PT(650,30,"c")
  0.001819560617
  >>> _Backward3x_v_PT(656,26,"d")
  0.002245587720
  >>> _Backward3x_v_PT(670,30,"d")
  0.002506897702
  >>> _Backward3x_v_PT(661,26,"e")
  0.002970225962
  >>> _Backward3x_v_PT(675,30,"e")
  0.003004627086
  >>> _Backward3x_v_PT(671,26,"f")
  0.005019029401
  >>> _Backward3x_v_PT(690,30,"f")
  0.004656470142
  >>> _Backward3x_v_PT(649,23.6,"g")
  0.002163198378
  >>> _Backward3x_v_PT(650,24,"g")
  0.002166044161
  >>> _Backward3x_v_PT(652,23.6,"h")
  0.002651081407
  >>> _Backward3x_v_PT(654,24,"h")
  0.002967802335
  >>> _Backward3x_v_PT(653,23.6,"i")
  0.003273916816
  >>> _Backward3x_v_PT(655,24,"i")
  0.003550329864
  >>> _Backward3x_v_PT(655,23.5,"j")
  0.004545001142
  >>> _Backward3x_v_PT(660,24,"j")
  0.005100267704
  >>> _Backward3x_v_PT(660,23,"k")
  0.006109525997
  >>> _Backward3x_v_PT(670,24,"k")
  0.006427325645
  >>> _Backward3x_v_PT(646,22.6,"l")
  0.002117860851
  >>> _Backward3x_v_PT(646,23,"l")
  0.002062374674
  >>> _Backward3x_v_PT(648.6,22.6,"m")
  0.002533063780
  >>> _Backward3x_v_PT(649.3,22.8,"m")
  0.002572971781
  >>> _Backward3x_v_PT(649,22.6,"n")
  0.002923432711
  >>> _Backward3x_v_PT(649.7,22.8,"n")
  0.002913311494
  >>> _Backward3x_v_PT(649.1,22.6,"o")
  0.003131208996
  >>> _Backward3x_v_PT(649.9,22.8,"o")
  0.003221160278
  >>> _Backward3x_v_PT(649.4,22.6,"p")
  0.003715596186
  >>> _Backward3x_v_PT(650.2,22.8,"p")
  0.003664754790
  >>> _Backward3x_v_PT(640,21.1,"q")
  0.001970999272
  >>> _Backward3x_v_PT(643,21.8,"q")
  0.002043919161
  >>> _Backward3x_v_PT(644,21.1,"r")
  0.005251009921
  >>> _Backward3x_v_PT(648,21.8,"r")
  0.005256844741
  >>> _Backward3x_v_PT(635,19.1,"s")
  0.001932829079
  >>> _Backward3x_v_PT(638,20,"s")
  0.001985387227
  >>> _Backward3x_v_PT(626,17,"t")
  0.008483262001
  >>> _Backward3x_v_PT(640,20,"t")
  0.006227528101
  >>> _Backward3x_v_PT(644.6,21.5,"u")
  0.002268366647
  >>> _Backward3x_v_PT(646.1,22,"u")
  0.002296350553
  >>> _Backward3x_v_PT(648.6,22.5,"v")
  0.002832373260
  >>> _Backward3x_v_PT(647.9,22.3,"v")
  0.002811424405
  >>> _Backward3x_v_PT(647.5,22.15,"w")
  0.003694032281
  >>> _Backward3x_v_PT(648.1,22.3,"w")
  0.003622226305
  >>> _Backward3x_v_PT(648,22.11,"x")
  0.004528072649
  >>> _Backward3x_v_PT(649,22.3,"x")
  0.004556905799
  >>> _Backward3x_v_PT(646.84,22,"y")
  0.002698354719
  >>> _Backward3x_v_PT(647.05,22.064,"y")
  0.002717655648
  >>> _Backward3x_v_PT(646.89,22,"z")
  0.003798732962
  >>> _Backward3x_v_PT(647.15,22.064,"z")
  0.003701940009
  """
  */

  par =
  {
    "a": [0.0024, 100, 760, 0.085, 0.817, 1, 1, 1],
    "b": [0.0041, 100, 860, 0.280, 0.779, 1, 1, 1],
    "c": [0.0022, 40, 690, 0.259, 0.903, 1, 1, 1],
    "d": [0.0029, 40, 690, 0.559, 0.939, 1, 1, 4],
    "e": [0.0032, 40, 710, 0.587, 0.918, 1, 1, 1],
    "f": [0.0064, 40, 730, 0.587, 0.891, 0.5, 1, 4],
    "g": [0.0027, 25, 660, 0.872, 0.971, 1, 1, 4],
    "h": [0.0032, 25, 660, 0.898, 0.983, 1, 1, 4],
    "i": [0.0041, 25, 660, 0.910, 0.984, 0.5, 1, 4],
    "j": [0.0054, 25, 670, 0.875, 0.964, 0.5, 1, 4],
    "k": [0.0077, 25, 680, 0.802, 0.935, 1, 1, 1],
    "l": [0.0026, 24, 650, 0.908, 0.989, 1, 1, 4],
    "m": [0.0028, 23, 650, 1.000, 0.997, 1, 0.25, 1],
    "n": [0.0031, 23, 650, 0.976, 0.997, None, None, None],
    "o": [0.0034, 23, 650, 0.974, 0.996, 0.5, 1, 1],
    "p": [0.0041, 23, 650, 0.972, 0.997, 0.5, 1, 1],
    "q": [0.0022, 23, 650, 0.848, 0.983, 1, 1, 4],
    "r": [0.0054, 23, 650, 0.874, 0.982, 1, 1, 1],
    "s": [0.0022, 21, 640, 0.886, 0.990, 1, 1, 4],
    "t": [0.0088, 20, 650, 0.803, 1.020, 1, 1, 1],
    "u": [0.0026, 23, 650, 0.902, 0.988, 1, 1, 1],
    "v": [0.0031, 23, 650, 0.960, 0.995, 1, 1, 1],
    "w": [0.0039, 23, 650, 0.959, 0.995, 1, 1, 4],
    "x": [0.0049, 23, 650, 0.910, 0.988, 1, 1, 1],
    "y": [0.0031, 22, 650, 0.996, 0.994, 1, 1, 4],
    "z": [0.0038, 22, 650, 0.993, 0.994, 1, 1, 4],
  }

  I =
  {
    "a": [-12, -12, -12, -10, -10, -10, -8, -8, -8, -6, -5, -5, -5, -4, -3, -3, -3, -3, -2, -2, -2, -1, -1, -1, 0, 0, 1, 1, 2, 2],
    "b": [-12, -12, -10, -10, -8, -6, -6, -6, -5, -5, -5, -4, -4, -4, -3, -3, -3, -3, -3, -2, -2, -2, -1, -1, 0, 0, 1, 1, 2, 3, 4, 4],
    "c": [-12, -12, -12, -10, -10, -10, -8, -8, -8, -6, -5, -5, -5, -4, -4, -3, -3, -2, -2, -2, -1, -1, -1, 0, 0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 8],
    "d": [-12, -12, -12, -12, -12, -12, -10, -10, -10, -10, -10, -10, -10, -8, -8, -8, -8, -6, -6, -5, -5, -5, -5, -4, -4, -4, -3, -3, -2, -2, -1, -1, -1, 0, 0, 1, 1, 3],
    "e": [-12, -12, -10, -10, -10, -10, -10, -8, -8, -8, -6, -5, -4, -4, -3, -3, -3, -2, -2, -2, -2, -1, 0, 0, 1, 1, 1, 2, 2],
    "f": [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 3, 3, 4, 5, 5, 6, 7, 7, 10, 12, 12, 12, 14, 14, 14, 14, 14, 16, 16, 18, 18, 20, 20, 20, 22, 24, 24, 28, 32],
    "g": [-12, -12, -12, -12, -12, -12, -10, -10, -10, -8, -8, -8, -8, -6, -6, -5, -5, -4, -3, -2, -2, -2, -2, -1, -1, -1, 0, 0, 0, 1, 1, 1, 3, 5, 6, 8, 10, 10],
    "h": [-12, -12, -10, -10, -10, -10, -10, -10, -8, -8, -8, -8, -8, -6, -6, -6, -5, -5, -5, -4, -4, -3, -3, -2, -1, -1, 0, 1, 1],
    "i": [0, 0, 0, 1, 1, 1, 1, 2, 3, 3, 4, 4, 4, 5, 5, 5, 7, 7, 8, 8, 10, 12, 12, 12, 14, 14, 14, 14, 18, 18, 18, 18, 18, 20, 20, 22, 24, 24, 32, 32, 36, 36],
    "j": [0, 0, 0, 1, 1, 1, 2, 2, 3, 4, 4, 5, 5, 5, 6, 10, 12, 12, 14, 14, 14, 16, 18, 20, 20, 24, 24, 28, 28],
    "k": [-2, -2, -1, -1, 0, -0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 5, 5, 5, 6, 6, 6, 6, 8, 10, 12],
    "l": [-12, -12, -12, -12, -12, -10, -10, -8, -8, -8, -8, -8, -8, -8, -6, -5, -5, -4, -4, -3, -3, -3, -3, -2, -2, -2, -1, -1, -1, 0, 0, 0, 0, 1, 1, 2, 4, 5, 5, 6, 10, 10, 14],
    "m": [0, 3, 8, 20, 1, 3, 4, 5, 1, 6, 2, 4, 14, 2, 5, 3, 0, 1, 1, 1, 28, 2, 16, 0, 5, 0, 3, 4, 12, 16, 1, 8, 14, 0, 2, 3, 4, 8, 14, 24],
    "n": [0, 3, 4, 6, 7, 10, 12, 14, 18, 0, 3, 5, 6, 8, 12, 0, 3, 7, 12, 2, 3, 4, 2, 4, 7, 4, 3, 5, 6, 0, 0, 3, 1, 0, 1, 0, 1, 0, 1],
    "o": [0, 0, 0, 2, 3, 4, 4, 4, 4, 4, 5, 5, 6, 7, 8, 8, 8, 10, 10, 14, 14, 20, 20, 24],
    "p": [0, 0, 0, 0, 1, 2, 3, 3, 4, 6, 7, 7, 8, 10, 12, 12, 12, 14, 14, 14, 16, 18, 20, 22, 24, 24, 36],
    "q": [-12, -12, -10, -10, -10, -10, -8, -6, -5, -5, -4, -4, -3, -2, -2, -2, -2, -1, -1, -1, 0, 1, 1, 1],
    "r": [-8, -8, -3, -3, -3, -3, -3, 0, 0, 0, 0, 3, 3, 8, 8, 8, 8, 10, 10, 10, 10, 10, 10, 10, 10, 12, 14],
    "s": [-12, -12, -10, -8, -6, -5, -5, -4, -4, -3, -3, -2, -1, -1, -1, 0, 0, 0, 0, 1, 1, 3, 3, 3, 4, 4, 4, 5, 14],
    "t": [0, 0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 4, 4, 7, 7, 7, 7, 7, 10, 10, 10, 10, 10, 18, 20, 22, 22, 24, 28, 32, 32, 32, 36],
    "u": [-12, -10, -10, -10, -8, -8, -8, -6, -6, -5, -5, -5, -3, -1, -1, -1, -1, 0, 0, 1, 2, 2, 3, 5, 5, 5, 6, 6, 8, 8, 10, 12, 12, 12, 14, 14, 14, 14],
    "v": [-10, -8, -6, -6, -6, -6, -6, -6, -5, -5, -5, -5, -5, -5, -4, -4, -4, -4, -3, -3, -3, -2, -2, -1, -1, 0, 0, 0, 1, 1, 3, 4, 4, 4, 5, 8, 10, 12, 14],
    "w": [-12, -12, -10, -10, -8, -8, -8, -6, -6, -6, -6, -5, -4, -4, -3, -3, -2, -2, -1, -1, -1, 0, 0, 1, 2, 2, 3, 3, 5, 5, 5, 8, 8, 10, 10],
    "x": [-8, -6, -5, -4, -4, -4, -3, -3, -1, 0, 0, 0, 1, 1, 2, 3, 3, 3, 4, 5, 5, 5, 6, 8, 8, 8, 8, 10, 12, 12, 12, 12, 14, 14, 14, 14],
    "y": [0, 0, 0, 0, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 8, 8, 10, 12],
    "z": [-8, -6, -5, -5, -4, -4, -4, -3, -3, -3, -2, -1, 0, 1, 2, 3, 3, 6, 6, 6, 6, 8, 8]
  }

  J =
  {
    "a": [5, 10, 12, 5, 10, 12, 5, 8, 10, 1, 1, 5, 10, 8, 0, 1, 3, 6, 0, 2, 3, 0, 1, 2, 0, 1, 0, 2, 0, 2],
    "b": [10, 12, 8, 14, 8, 5, 6, 8, 5, 8, 10, 2, 4, 5, 0, 1, 2, 3, 5, 0, 2, 5, 0, 2, 0, 1, 0, 2, 0, 2, 0, 1],
    "c": [6, 8, 10, 6, 8, 10, 5, 6, 7, 8, 1, 4, 7, 2, 8, 0, 3, 0, 4, 5, 0, 1, 2, 0, 1, 2, 0, 2, 0, 1, 3, 7, 0, 7, 1],
    "d": [4, 6, 7, 10, 12, 16, 0, 2, 4, 6, 8, 10, 14, 3, 7, 8, 10, 6, 8, 1, 2, 5, 7, 0, 1, 7, 2, 4, 0, 1, 0, 1, 5, 0, 2, 0, 6, 0],
    "e": [14, 16, 3, 6, 10, 14, 16, 7, 8, 10, 6, 6, 2, 4, 2, 6, 7, 0, 1, 3, 4, 0, 0, 1, 0, 4, 6, 0, 2],
    "f": [-3, -2, -1, 0, 1, 2, -1, 1, 2, 3, 0, 1, -5, -2, 0, -3, -8, 1, -6, -4, 1, -6, -10, -8, -4, -12, -10, -8, -6, -4, -10, -8, -12, -10, -12, -10, -6, -12, -12, -4, -12, -12],
    "g": [7, 12, 14, 18, 22, 24, 14, 20, 24, 7, 8, 10, 12, 8, 22, 7, 20, 22, 7, 3, 5, 14, 24, 2, 8, 18, 0, 1, 2, 0, 1, 3, 24, 22, 12, 3, 0, 6],
    "h": [8, 12, 4, 6, 8, 10, 14, 16, 0, 1, 6, 7, 8, 4, 6, 8, 2, 3, 4, 2, 4, 1, 2, 0, 0, 2, 0, 0, 2],
    "i": [0, 1, 10, -4, -2, -1, 0, 0, -5, 0, -3, -2, -1, -6, -1, 12, -4, -3, -6, 10, -8, -12, -6, -4, -10, -8, -4, 5, -12, -10, -8, -6, 2, -12, -10, -12, -12, -8, -10, -5, -10, -8],
    "j": [-1, 0, 1, -2, -1, 1, -1, 1, -2, -2, 2, -3, -2, 0, 3, -6, -8, -3, -10, -8, -5, -10, -12, -12, -10, -12, -6, -12, -5],
    "k": [10, 12, -5, 6, -12, -6, -2, -1, 0, 1, 2, 3, 14, -3, -2, 0, 1, 2, -8, -6, -3, -2, 0, 4, -12, -6, -3, -12, -10, -8, -5, -12, -12, -10],
    "l": [14, 16, 18, 20, 22, 14, 24, 6, 10, 12, 14, 18, 24, 36, 8, 4, 5, 7, 16, 1, 3, 18, 20, 2, 3, 10, 0, 1, 3, 0, 1, 2, 12, 0, 16, 1, 0, 0, 1, 14, 4, 12, 10],
    "m": [0, 0, 0, 2, 5, 5, 5, 5, 6, 6, 7, 8, 8, 10, 10, 12, 14, 14, 18, 20, 20, 22, 22, 24, 24, 28, 28, 28, 28, 28, 32, 32, 32, 36, 36, 36, 36, 36, 36, 36],
    "n": [-12, -12, -12, -12, -12, -12, -12, -12, -12, -10, -10, -10, -10, -10, -10, -8, -8, -8, -8, -6, -6, -6, -5, -5, -5, -4, -3, -3, -3, -2, -1, -1, 0, 1, 1, 2, 4, 5, 6],
    "o": [-12, -4, -1, -1, -10, -12, -8, -5, -4, -1, -4, -3, -8, -12, -10, -8, -4, -12, -8, -12, -8, -12, -10, -12],
    "p": [-1, 0, 1, 2, 1, -1, -3, 0, -2, -2, -5, -4, -2, -3, -12, -6, -5, -10, -8, -3, -8, -8, -10, -10, -12, -8, -12],
    "q": [10, 12, 6, 7, 8, 10, 8, 6, 2, 5, 3, 4, 3, 0, 1, 2, 4, 0, 1, 2, 0, 0, 1, 3],
    "r": [6, 14, -3, 3, 4, 5, 8, -1, 0, 1, 5, -6, -2, -12, -10, -8, -5, -12, -10, -8, -6, -5, -4, -3, -2, -12, -12],
    "s": [20, 24, 22, 14, 36, 8, 16, 6, 32, 3, 8, 4, 1, 2, 3, 0, 1, 4, 28, 0, 32, 0, 1, 2, 3, 18, 24, 4, 24],
    "t": [0, 1, 4, 12, 0, 10, 0, 6, 14, 3, 8, 0, 10, 3, 4, 7, 20, 36, 10, 12, 14, 16, 22, 18, 32, 22, 36, 24, 28, 22, 32, 36, 36],
    "u": [14, 10, 12, 14, 10, 12, 14, 8, 12, 4, 8, 12, 2, -1, 1, 12, 14, -3, 1, -2, 5, 10, -5, -4, 2, 3, -5, 2, -8, 8, -4, -12, -4, 4, -12, -10, -6, 6],
    "v": [-8, -12, -12, -3, 5, 6, 8, 10, 1, 2, 6, 8, 10, 14, -12, -10, -6, 10, -3, 10, 12, 2, 4, -2, 0, -2, 6, 10, -12, -10, 3, -6, 3, 10, 2, -12, -2, -3, 1],
    "w": [8, 14, -1, 8, 6, 8, 14, -4, -3, 2, 8, -10, -1, 3, -10, 3, 1, 2, -8, -4, 1, -12, 1, -1, -1, 2, -12, -5, -10, -8, -6, -12, -10, -12, -8],
    "x": [14, 10, 10, 1, 2, 14, -2, 12, 5, 0, 4, 10, -10, -1, 6, -12, 0, 8, 3, -6, -2, 1, 1, -6, -3, 1, 8, -8, -10, -8, -5, -4, -12, -10, -8, -6],
    "y": [-3, 1, 5, 8, 8, -4, -1, 4, 5, -8, 4, 8, -6, 6, -2, 1, -8, -2, -5, -8],
    "z": [3, 6, 6, 8, 5, 6, 8, -2, 5, 6, 2, -6, 3, 1, 6, -6, -2, -6, -5, -4, -1, -8, -4]
  }

  n =
  {
    "a": [0.110879558823853e-2, 0.572616740810616e3, -0.767051948380852e5, -0.253321069529674e-1, 0.628008049345689e4, 0.234105654131876e6, 0.216867826045856, -0.156237904341963e3, -0.269893956176613e5, -0.180407100085505e-3, 0.116732227668261e-2, 0.266987040856040e2, 0.282776617243286e5, -0.242431520029523e4, 0.435217323022733e-3, -0.122494831387441e-1, 0.179357604019989e1, 0.442729521058314e2, -0.593223489018342e-2, 0.453186261685774, 0.135825703129140e1, 0.408748415856745e-1, 0.474686397863312, 0.118646814997915e1, 0.546987265727549, 0.195266770452643, -0.502268790869663e-1, -0.369645308193377, 0.633828037528420e-2, 0.797441793901017e-1],
    "b": [-0.827670470003621e-1, 0.416887126010565e2, 0.483651982197059e-1, -0.291032084950276e5, -0.111422582236948e3, -.202300083904014e-1, 0.294002509338515e3, 0.140244997609658e3, -0.344384158811459e3, 0.361182452612149e3, -0.140699677420738e4, -0.202023902676481e-2, 0.171346792457471e3, -0.425597804058632e1, 0.691346085000334e-5, 0.151140509678925e-2, -0.416375290166236e-1, -.413754957011042e2, -0.506673295721637e2, -0.572212965569023e-3, 0.608817368401785e1, 0.239600660256161e2, 0.122261479925384e-1, 0.216356057692938e1, 0.398198903368642, -0.116892827834085, -0.102845919373532, -0.492676637589284, 0.655540456406790e-1, -0.240462535078530, -0.269798180310075e-1, 0.128369435967012],
    "c": [0.311967788763030e1, 0.276713458847564e5, 0.322583103403269e8, -0.342416065095363e3, -0.899732529907377e6, -0.793892049821251e8, 0.953193003217388e2, 0.229784742345072e4, 0.175336675322499e6, 0.791214365222792e7, 0.319933345844209e-4, -0.659508863555767e2, -0.833426563212851e6, 0.645734680583292e-1, -0.382031020570813e7, 0.406398848470079e-4, 0.310327498492008e2, -0.892996718483724e-3, 0.234604891591616e3, 0.377515668966951e4, 0.158646812591361e-1, 0.707906336241843, 0.126016225146570e2, 0.736143655772152, 0.676544268999101, -0.178100588189137e2, -0.156531975531713, 0.117707430048158e2, 0.840143653860447e-1, -0.186442467471949, -0.440170203949645e2, 0.123290423502494e7, -0.240650039730845e-1, -0.107077716660869e7, 0.438319858566475e-1],
    "d": [-0.452484847171645e-9, .315210389538801e-4, -.214991352047545e-2, 0.508058874808345e3, -0.127123036845932e8, 0.115371133120497e13, -.197805728776273e-15, .241554806033972e-10, -.156481703640525e-5, 0.277211346836625e-2, -0.203578994462286e2, 0.144369489909053e7, -0.411254217946539e11, 0.623449786243773e-5, -.221774281146038e2, -0.689315087933158e5, -0.195419525060713e8, 0.316373510564015e4, 0.224040754426988e7, -0.436701347922356e-5, -.404213852833996e-3, -0.348153203414663e3, -0.385294213555289e6, 0.135203700099403e-6, 0.134648383271089e-3, 0.125031835351736e6, 0.968123678455841e-1, 0.225660517512438e3, -0.190102435341872e-3, -.299628410819229e-1, 0.500833915372121e-2, 0.387842482998411, -0.138535367777182e4, 0.870745245971773, 0.171946252068742e1, -.326650121426383e-1, 0.498044171727877e4, 0.551478022765087e-2],
    "e": [0.715815808404721e9, -0.114328360753449e12, .376531002015720e-11, -0.903983668691157e-4, 0.665695908836252e6, 0.535364174960127e10, 0.794977402335603e11, 0.922230563421437e2, -0.142586073991215e6, -0.111796381424162e7, 0.896121629640760e4, -0.669989239070491e4, 0.451242538486834e-2, -0.339731325977713e2, -0.120523111552278e1, 0.475992667717124e5, -0.266627750390341e6, -0.153314954386524e-3, 0.305638404828265, 0.123654999499486e3, -0.104390794213011e4, -0.157496516174308e-1, 0.685331118940253, 0.178373462873903e1, -0.544674124878910, 0.204529931318843e4, -0.228342359328752e5, 0.413197481515899, -0.341931835910405e2],
    "f": [-0.251756547792325e-7, .601307193668763e-5, -.100615977450049e-2, 0.999969140252192, 0.214107759236486e1, -0.165175571959086e2, -0.141987303638727e-2, 0.269251915156554e1, 0.349741815858722e2, -0.300208695771783e2, -0.131546288252539e1, -0.839091277286169e1, 0.181545608337015e-9, -0.591099206478909e-3, 0.152115067087106e1, 0.252956470663225e-4, 0.100726265203786e-14, -0.14977453386065e1, -0.793940970562969e-9, -0.150290891264717e-3, .151205531275133e1, 0.470942606221652e-5, .195049710391712e-12, -.911627886266077e-8, .604374640201265e-3, -.225132933900136e-15, .610916973582981e-11, -.303063908043404e-6, -.137796070798409e-4, -.919296736666106e-3, .639288223132545e-9, .753259479898699e-6, -0.400321478682929e-12, .756140294351614e-8, -.912082054034891e-11, -.237612381140539e-7, 0.269586010591874e-4, -.732828135157839e-10, .241995578306660e-9, -.405735532730322e-3, .189424143498011e-9, -.486632965074563e-9],
    "g": [0.412209020652996e-4, -0.114987238280587e7, 0.948180885032080e10, -0.195788865718971e18, 0.4962507048713e25, -0.105549884548496e29, -0.758642165988278e12, -.922172769596101e23, .725379072059348e30, -0.617718249205859e2, 0.107555033344858e5, -0.379545802336487e8, 0.228646846221831e12, -0.499741093010619e7, -.280214310054101e31, 0.104915406769586e7, 0.613754229168619e28, 0.802056715528378e32, -0.298617819828065e8, -0.910782540134681e2, 0.135033227281565e6, -0.712949383408211e19, -0.104578785289542e37, .304331584444093e2, 0.593250797959445e10, -0.364174062110798e28, 0.921791403532461, -0.337693609657471, -0.724644143758508e2, -0.110480239272601, 0.536516031875059e1, -0.291441872156205e4, 0.616338176535305e40, -0.120889175861180e39, 0.818396024524612e23, 0.940781944835829e9, -0.367279669545448e5, -0.837513931798655e16],
    "h": [0.561379678887577e-1, 0.774135421587083e10, 0.111482975877938e-8, -0.143987128208183e-2, 0.193696558764920e4, -0.605971823585005e9, 0.171951568124337e14, -.185461154985145e17, 0.38785116807801e-16, -.395464327846105e-13, -0.170875935679023e3, -0.21201062070122e4, 0.177683337348191e8, 0.110177443629575e2, -0.234396091693313e6, -0.656174421999594e7, 0.156362212977396e-4, -0.212946257021400e1, 0.135249306374858e2, 0.177189164145813, 0.139499167345464e4, -0.703670932036388e-2, -0.152011044389648, 0.981916922991113e-4, 0.147199658618076e-2, 0.202618487025578e2, 0.899345518944240, -0.211346402240858, 0.249971752957491e2],
    "i": [0.106905684359136e1, -0.148620857922333e1, 0.259862256980408e15, -.446352055678749e-11, -.566620757170032e-6, -.235302885736849e-2, -0.269226321968839, 0.922024992944392e1, 0.357633505503772e-11, -.173942565562222e2, 0.700681785556229e-5, -.267050351075768e-3, -.231779669675624e1, -.753533046979752e-12, .481337131452891e1, -0.223286270422356e22, -.118746004987383e-4, .646412934136496e-2, -0.410588536330937e-9, .422739537057241e20, .313698180473812e-12, 0.16439533434504e-23, -.339823323754373e-5, -.135268639905021e-1, -.723252514211625e-14, .184386437538366e-8, -.463959533752385e-1, -.99226310037675e14, .688169154439335e-16, -.222620998452197e-10, -.540843018624083e-7, .345570606200257e-2, .422275800304086e11, -.126974478770487e-14, .927237985153679e-9, .612670812016489e-13, -.722693924063497e-11, -.383669502636822e-3, .374684572410204e-3, -0.931976897511086e5, -0.247690616026922e-1, .658110546759474e2],
    "j": [-0.111371317395540e-3, 0.100342892423685e1, 0.530615581928979e1, 0.179058760078792e-5, -0.728541958464774e-3, -.187576133371704e2, 0.199060874071849e-2, 0.243574755377290e2, -0.177040785499444e-3, -0.25968038522713e-2, -0.198704578406823e3, 0.738627790224287e-4, -0.236264692844138e-2, -0.161023121314333e1, 0.622322971786473e4, -.960754116701669e-8, -.510572269720488e-10, .767373781404211e-2, .663855469485254e-14, -.717590735526745e-9, 0.146564542926508e-4, .309029474277013e-11, -.464216300971708e-15, -.390499637961161e-13, -.236716126781431e-9, .454652854268717e-11, -.422271787482497e-2, 0.283911742354706e-10, 0.270929002720228e1],
    "k": [-0.401215699576099e9, 0.484501478318406e11, .394721471363678e-14, .372629967374147e5, -.369794374168666e-29, -.380436407012452e-14, 0.475361629970233e-6, -0.879148916140706e-3, 0.844317863844331, 0.122433162656600e2, -0.104529634830279e3, 0.589702771277429e3, -.291026851164444e14, .170343072841850e-5, -0.277617606975748e-3, -0.344709605486686e1, 0.221333862447095e2, -0.194646110037079e3, .808354639772825e-15, -.18084520914547e-10, -.696664158132412e-5, -0.181057560300994e-2, 0.255830298579027e1, 0.328913873658481e4, -.173270241249904e-18, -.661876792558034e-6, -.39568892342125e-2, .604203299819132e-17, -.400879935920517e-13, .160751107464958e-8, .383719409025556e-4, -.649565446702457e-14, -.149095328506e-11, 0.541449377329581e-8],
    "l": [0.260702058647537e10, -.188277213604704e15, 0.554923870289667e19, -.758966946387758e23, .413865186848908e27, -.81503800073806e12, -.381458260489955e33, -.123239564600519e-1, 0.226095631437174e8, -.49501780950672e12, 0.529482996422863e16, -0.444359478746295e23, .521635864527315e35, -0.487095672740742e55, -0.714430209937547e6, 0.127868634615495, -0.100752127917598e2, 0.777451437960990e7, -.108105480796471e25, -.357578581169659e-5, -0.212857169423484e1, 0.270706111085238e30, -0.695953622348829e33, 0.110609027472280, 0.721559163361354e2, -0.306367307532219e15, 0.265839618885530e-4, 0.253392392889754e-1, -0.214443041836579e3, 0.937846601489667, 0.223184043101700e1, 0.338401222509191e2, 0.494237237179718e21, -0.198068404154428, -0.141415349881140e31, -0.993862421613651e2, 0.125070534142731e3, -0.996473529004439e3, 0.473137909872765e5, 0.116662121219322e33, -0.315874976271533e16, -0.445703369196945e33, 0.642794932373694e33],
    "m": [0.811384363481847, -0.568199310990094e4, -0.178657198172556e11, 0.795537657613427e32, -0.814568209346872e5, -0.659774567602874e8, -.152861148659302e11, -0.560165667510446e12, 0.458384828593949e6, -0.385754000383848e14, 0.453735800004273e8, 0.939454935735563e12, .266572856432938e28, -0.547578313899097e10, 0.200725701112386e15, 0.185007245563239e13, 0.185135446828337e9, -0.170451090076385e12, 0.157890366037614e15, -0.202530509748774e16, 0.36819392618357e60, 0.170215539458936e18, 0.639234909918741e42, -.821698160721956e15, -.795260241872306e24, 0.23341586947851e18, -0.600079934586803e23, 0.594584382273384e25, 0.189461279349492e40, -.810093428842645e46, 0.188813911076809e22, 0.111052244098768e36, 0.291133958602503e46, -.329421923951460e22, -.137570282536696e26, 0.181508996303902e28, -.346865122768353e30, -.21196114877426e38, -0.128617899887675e49, 0.479817895699239e65],
    "n": [.280967799943151e-38, .614869006573609e-30, .582238667048942e-27, .390628369238462e-22, .821445758255119e-20, .402137961842776e-14, .651718171878301e-12, -.211773355803058e-7, 0.264953354380072e-2, -.135031446451331e-31, -.607246643970893e-23, -.402352115234494e-18, -.744938506925544e-16, .189917206526237e-12, .364975183508473e-5, .177274872361946e-25, -.334952758812999e-18, -.421537726098389e-8, -.391048167929649e-1, .541276911564176e-13, .705412100773699e-11, .258585887897486e-8, -.493111362030162e-10, -.158649699894543e-5, -0.525037427886100, 0.220019901729615e-2, -0.643064132636925e-2, 0.629154149015048e2, 0.135147318617061e3, 0.240560808321713e-6, -.890763306701305e-3, -0.440209599407714e4, -0.302807107747776e3, 0.159158748314599e4, 0.232534272709876e6, -0.792681207132600e6, -.869871364662769e11, .354542769185671e12, 0.400849240129329e15],
    "o": [.128746023979718e-34, -.735234770382342e-11, .28907869214915e-2, 0.244482731907223, 0.141733492030985e-23, -0.354533853059476e-28, -.594539202901431e-17, -.585188401782779e-8, .201377325411803e-5, 0.138647388209306e1, -0.173959365084772e-4, 0.137680878349369e-2, .814897605805513e-14, .425596631351839e-25, -.387449113787755e-17, .13981474793024e-12, -.171849638951521e-2, 0.641890529513296e-21, .118960578072018e-10, -.155282762571611e-17, .233907907347507e-7, -.174093247766213e-12, .377682649089149e-8, -.516720236575302e-10],
    "p": [-0.982825342010366e-4, 0.105145700850612e1, 0.116033094095084e3, 0.324664750281543e4, -0.123592348610137e4, -0.561403450013495e-1, 0.856677401640869e-7, 0.236313425393924e3, 0.972503292350109e-2, -.103001994531927e1, -0.149653706199162e-8, -.215743778861592e-4, -0.834452198291445e1, 0.586602660564988, 0.343480022104968e-25, .816256095947021e-5, .294985697916798e-2, 0.711730466276584e-16, 0.400954763806941e-9, 0.107766027032853e2, -0.409449599138182e-6, -.729121307758902e-5, 0.677107970938909e-8, 0.602745973022975e-7, -.382323011855257e-10, .179946628317437e-2, -.345042834640005e-3],
    "q": [-0.820433843259950e5, 0.473271518461586e11, -.805950021005413e-1, 0.328600025435980e2, -0.35661702998249e4, -0.172985781433335e10, 0.351769232729192e8, -0.775489259985144e6, 0.710346691966018e-4, 0.993499883820274e5, -0.642094171904570, -0.612842816820083e4, .232808472983776e3, -0.142808220416837e-4, -0.643596060678456e-2, -0.428577227475614e1, 0.225689939161918e4, 0.100355651721510e-2, 0.333491455143516, 0.109697576888873e1, 0.961917379376452, -0.838165632204598e-1, 0.247795908411492e1, -.319114969006533e4],
    "r": [.144165955660863e-2, -.701438599628258e13, -.830946716459219e-16, 0.261975135368109, 0.393097214706245e3, -0.104334030654021e5, 0.490112654154211e9, -0.147104222772069e-3, 0.103602748043408e1, 0.305308890065089e1, -0.399745276971264e7, 0.569233719593750e-11, -.464923504407778e-1, -.535400396512906e-17, .399988795693162e-12, -.536479560201811e-6, .159536722411202e-1, .270303248860217e-14, .244247453858506e-7, -0.983430636716454e-5, 0.663513144224454e-1, -0.993456957845006e1, 0.546491323528491e3, -0.143365406393758e5, 0.150764974125511e6, -.337209709340105e-9, 0.377501980025469e-8],
    "s": [-0.532466612140254e23, .100415480000824e32, -.191540001821367e30, 0.105618377808847e17, 0.202281884477061e59, 0.884585472596134e8, 0.166540181638363e23, -0.313563197669111e6, -.185662327545324e54, -.624942093918942e-1, -0.50416072413259e10, 0.187514491833092e5, 0.121399979993217e-2, 0.188317043049455e1, -0.167073503962060e4, 0.965961650599775, 0.294885696802488e1, -0.653915627346115e5, 0.604012200163444e50, -0.198339358557937, -0.175984090163501e58, 0.356314881403987e1, -0.575991255144384e3, 0.456213415338071e5, -.109174044987829e8, 0.437796099975134e34, -0.616552611135792e46, 0.193568768917797e10, 0.950898170425042e54],
    "t": [0.155287249586268e1, 0.664235115009031e1, -0.289366236727210e4, -0.385923202309848e13, -.291002915783761e1, -.829088246858083e12, 0.176814899675218e1, -0.534686695713469e9, 0.160464608687834e18, 0.196435366560186e6, 0.156637427541729e13, -0.178154560260006e1, -0.229746237623692e16, 0.385659001648006e8, 0.110554446790543e10, -.677073830687349e14, -.327910592086523e31, -.341552040860644e51, -.527251339709047e21, .245375640937055e24, -0.168776617209269e27, .358958955867578e29, -0.656475280339411e36, 0.355286045512301e39, .569021454413270e58, -.700584546433113e48, -0.705772623326374e65, 0.166861176200148e53, -.300475129680486e61, -.668481295196808e51, .428432338620678e69, -.444227367758304e72, -.281396013562745e77],
    "u": [0.122088349258355e18, 0.104216468608488e10, -.882666931564652e16, .259929510849499e20, 0.222612779142211e15, -0.878473585050085e18, -0.314432577551552e22, -.216934916996285e13, .159079648196849e21, -.339567617303423e3, 0.884387651337836e13, -0.843405926846418e21, 0.114178193518022e2, -0.122708229235641e-3, -0.106201671767107e3, .903443213959313e25, -0.693996270370852e28, 0.648916718965575e-8, 0.718957567127851e4, 0.105581745346187e-2, -0.651903203602581e15, -0.160116813274676e25, -0.510254294237837e-8, -0.152355388953402, 0.677143292290144e12, 0.276378438378930e15, 0.116862983141686e-1, -.301426947980171e14, 0.169719813884840e-7, 0.104674840020929e27, -0.10801690456014e5, -0.990623601934295e-12, 0.536116483602738e7, .226145963747881e22, -0.488731565776210e-9, 0.151001548880670e-4, -0.227700464643920e5, -0.781754507698846e28],
    "v": [-.415652812061591e-54, .177441742924043e-60, -.357078668203377e-54, 0.359252213604114e-25, -0.259123736380269e2, 0.594619766193460e5, -0.624184007103158e11, 0.313080299915944e17, .105006446192036e-8, -0.192824336984852e-5, 0.654144373749937e6, 0.513117462865044e13, -.697595750347391e19, -.103977184454767e29, .119563135540666e-47, -.436677034051655e-41, .926990036530639e-29, .587793105620748e21, .280375725094731e-17, -0.192359972440634e23, .742705723302738e27, -0.517429682450605e2, 0.820612048645469e7, -0.188214882341448e-8, .184587261114837e-1, -0.135830407782663e-5, -.723681885626348e17, -.223449194054124e27, -.111526741826431e-34, .276032601145151e-28, 0.134856491567853e15, 0.652440293345860e-9, 0.510655119774360e17, -.468138358908732e32, -.760667491183279e16, -.417247986986821e-18, 0.312545677756104e14, -.100375333864186e15, .247761392329058e27],
    "w": [-.586219133817016e-7, -.894460355005526e11, .531168037519774e-30, 0.109892402329239, -0.575368389425212e-1, 0.228276853990249e5, -.158548609655002e19, .329865748576503e-27, -.634987981190669e-24, 0.615762068640611e-8, -.961109240985747e8, -.406274286652625e-44, -0.471103725498077e-12, 0.725937724828145, .187768525763682e-38, -.103308436323771e4, -0.662552816342168e-1, 0.579514041765710e3, .237416732616644e-26, .271700235739893e-14, -0.9078862134836e2, -0.171242509570207e-36, 0.156792067854621e3, 0.923261357901470, -0.597865988422577e1, 0.321988767636389e7, -.399441390042203e-29, .493429086046981e-7, .812036983370565e-19, -.207610284654137e-11, -.340821291419719e-6, .542000573372233e-17, -.856711586510214e-12, 0.266170454405981e-13, 0.858133791857099e-5],
    "x": [.377373741298151e19, -.507100883722913e13, -0.10336322559886e16, .184790814320773e-5, -.924729378390945e-3, -0.425999562292738e24, -.462307771873973e-12, .107319065855767e22, 0.648662492280682e11, 0.244200600688281e1, -0.851535733484258e10, 0.169894481433592e22, 0.215780222509020e-26, -0.320850551367334, -0.382642448458610e17, -.275386077674421e-28, -.563199253391666e6, -.326068646279314e21, 0.397949001553184e14, 0.100824008584757e-6, 0.162234569738433e5, -0.432355225319745e11, -.59287424559861e12, 0.133061647281106e1, 0.157338197797544e7, 0.258189614270853e14, 0.262413209706358e25, -.920011937431142e-1, 0.220213765905426e-2, -0.110433759109547e2, 0.847004870612087e7, -0.592910695762536e9, -0.183027173269660e-4, 0.181339603516302, -0.119228759669889e4, 0.430867658061468e7],
    "y": [-0.525597995024633e-9, 0.583441305228407e4, -.134778968457925e17, .118973500934212e26, -0.159096490904708e27, -.315839902302021e-6, 0.496212197158239e3, 0.327777227273171e19, -0.527114657850696e22, .210017506281863e-16, 0.705106224399834e21, -.266713136106469e31, -0.145370512554562e-7, 0.149333917053130e28, -.149795620287641e8, -.3818819062711e16, 0.724660165585797e-4, -0.937808169550193e14, 0.514411468376383e10, -0.828198594040141e5],
    "z": [0.24400789229065e-10, -0.463057430331242e7, 0.728803274777712e10, .327776302858856e16, -.110598170118409e10, -0.323899915729957e13, .923814007023245e16, 0.842250080413712e-12, 0.663221436245506e12, -.167170186672139e15, .253749358701391e4, -0.819731559610523e-20, 0.328380587890663e12, -0.625004791171543e8, 0.803197957462023e21, -.204397011338353e-10, -.378391047055938e4, 0.97287654593862e-2, 0.154355721681459e2, -0.373962862928643e4, -0.682859011374572e11, -0.248488015614543e-3, 0.394536049497068e7]
  }

  I = I[x]
  J = J[x]
  n = n[x]
  v_, P_, T_, a, b, c, d, e = par[x]

  Pr = P/P_
  Tr = T/T_
  suma = 0
  if(x == "n")
  {
    for(int iter = 0; iter < I.length; iter++)
    {
      i = I[iter];
      j = J[iter];
      ni = n[iter];
      suma += ni * (Pr-a)**i * (Tr-b)**j
    }
    return v_*exp(suma)
  }
  else
  {
    for(int iter = 0; iter < I.length; iter++)
    {
      i = I[iter];
      j = J[iter];
      ni = n[iter];
      suma += ni * (Pr-a)**(c*i) * (Tr-b)**(j*d)
    }
    return v_*suma**e
  }
}

// Region 4
void _Region4(P, x)
{
  /*
  """Basic equation for region 4

  Parameters
  ----------
  P : float
    Pressure, [MPa]
  x : float
    Vapor quality, [-]

  Returns
  -------
  prop : dict
    Dict with calculated properties. The available properties are:

      * T: Saturated temperature, [K]
      * P: Saturated pressure, [MPa]
      * x: Vapor quality, [-]
      * v: Specific volume, [m³/kg]
      * h: Specific enthalpy, [kJ/kg]
      * s: Specific entropy, [kJ/kgK]
  """
  */

  T = _TSat_P(P)
  if(T > 623.15)
  {
    rhol = 1./_Backward3_sat_v_P(P, T, 0)
    P1 = _Region3(rhol, T)
    rhov = 1./_Backward3_sat_v_P(P, T, 1)
    P2 = _Region3(rhov, T)
  }
  else
  {
    P1 = _Region1(T, P)
    P2 = _Region2(T, P)
  }

  propiedades = {}
  propiedades["T"] = T
  propiedades["P"] = P
  propiedades["v"] = P1["v"]+x*(P2["v"]-P1["v"])
  propiedades["h"] = P1["h"]+x*(P2["h"]-P1["h"])
  propiedades["s"] = P1["s"]+x*(P2["s"]-P1["s"])
  propiedades["cp"] = None
  propiedades["cv"] = None
  propiedades["w"] = None
  propiedades["alfav"] = None
  propiedades["kt"] = None
  propiedades["region"] = 4
  propiedades["x"] = x
  return propiedades
}

void _Backward4_T_hs(h, s)
{
  /*
  """Backward equation for region 4, T=f(h,s)

  Parameters
  ----------
  h : float
    Specific enthalpy, [kJ/kg]
  s : float
    Specific entropy, [kJ/kgK]

  Returns
  -------
  T : float
    Temperature, [K]

  References
  ----------
  IAPWS, Revised Supplementary Release on Backward Equations p(h,s) for
  Region 3, Equations as a Function of h and s for the Region Boundaries, and
  an Equation Tsat(h,s) for Region 4 of the IAPWS Industrial Formulation 1997
  for the Thermodynamic Properties of Water and Steam,
  http://www.iapws.org/relguide/Supp-phs3-2014.pdf. Eq 9

  Examples
  --------
  >>> _Backward4_T_hs(1800,5.3)
  346.8475498
  >>> _Backward4_T_hs(2400,6.0)
  425.1373305
  >>> _Backward4_T_hs(2500,5.5)
  522.5579013
  """
  */

  I = [0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 5, 5, 6, 6, 6, 8, 10, 10, 12, 14, 14, 16, 16, 18, 18, 18, 20, 28]
  J = [0, 3, 12, 0, 1, 2, 5, 0, 5, 8, 0, 2, 3, 4, 0, 1, 1, 2, 4, 16, 6, 8, 22, 1, 20, 36, 24, 1, 28, 12, 32, 14, 22, 36, 24, 36]
  n = [0.179882673606601, -0.267507455199603, 0.116276722612600e1, 0.147545428713616, -0.512871635973248, 0.421333567697984, 0.563749522189870, 0.429274443819153, -0.335704552142140e1, 0.108890916499278e2, -0.248483390456012, 0.304153221906390, -0.494819763939905, 0.107551674933261e1, 0.733888415457688e-1, 0.140170545411085e-1, -0.106110975998808, 0.168324361811875e-1, 0.125028363714877e1, 0.101316840309509e4, -0.151791558000712e1, 0.524277865990866e2, 0.230495545563912e5, 0.249459806365456e-1, 0.210796467412137e7, 0.366836848613065e9, -0.144814105365163e9, -0.179276373003590e-2, 0.489955602100459e10, 0.471262212070518e3, -0.829294390198652e11, -0.171545662263191e4, 0.355777682973575e7, 0.586062760258436e12, -0.129887635078195e8, 0.317247449371057e11]

  nu = h/2800
  sigma = s/9.2
  suma = 0
  for(int iter = 0; iter < I.length; iter++)
  {
    i = I[iter];
    j = J[iter];
    ni = n[iter];
    suma += ni * (nu-0.119)**i * (sigma-1.07)**j
  }
  return 550*suma
}

// Region 5
void _Region5(T, P)
{
  /*
  """Basic equation for region 5

  Parameters
  ----------
  T : float
    Temperature, [K]
  P : float
    Pressure, [MPa]

  Returns
  -------
  prop : dict
    Dict with calculated properties. The available properties are:

      * v: Specific volume, [m³/kg]
      * h: Specific enthalpy, [kJ/kg]
      * s: Specific entropy, [kJ/kgK]
      * cp: Specific isobaric heat capacity, [kJ/kgK]
      * cv: Specific isocoric heat capacity, [kJ/kgK]
      * w: Speed of sound, [m/s]
      * alfav: Cubic expansion coefficient, [1/K]
      * kt: Isothermal compressibility, [1/MPa]

  References
  ----------
  IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  Thermodynamic Properties of Water and Steam August 2007,
  http://www.iapws.org/relguide/IF97-Rev.html, Eq 32-34

  Examples
  --------
  >>> _Region5(1500,0.5)["v"]
  1.38455090
  >>> _Region5(1500,0.5)["h"]
  5219.76855
  >>> _Region5(1500,0.5)["h"]-500*_Region5(1500,0.5)["v"]
  4527.49310
  >>> _Region5(1500,30)["s"]
  7.72970133
  >>> _Region5(1500,30)["cp"]
  2.72724317
  >>> _Region5(1500,30)["cv"]
  2.19274829
  >>> _Region5(2000,30)["w"]
  1067.36948
  >>> _Region5(2000,30)["alfav"]
  0.000508830641
  >>> _Region5(2000,30)["kt"]
  0.0329193892
  """
  */

  if(P < 0)
  {
    P = Pmin

  }
  Tr = 1000/T
  Pr = P/1

  go, gop, gopp, got, gott, gopt = Region5_cp0(Tr, Pr)

  Ir = [1, 1, 1, 2, 2, 3]
  Jr = [1, 2, 3, 3, 9, 7]
  nr = [0.15736404855259e-2, 0.90153761673944e-3, -0.50270077677648e-2, 0.22440037409485e-5, -0.41163275453471e-5, 0.37919454822955e-7]
  gr = grp = grpp = grt = grtt = grpt = 0
  for(int iter = 0; iter < Ir.length; iter++)
  {
    i = Ir[iter];
    j = Jr[iter];
    ni = nr[iter];
    gr += ni * Pr**i * Tr**j
    grp += ni*i * Pr**(i-1) * Tr**j
    grpp += ni*i*(i-1) * Pr**(i-2) * Tr**j
    grt += ni*j * Pr**i * Tr**(j-1)
    grtt += ni*j*(j-1) * Pr**i * Tr**(j-2)
    grpt += ni*i*j * Pr**(i-1) * Tr**(j-1)
  }

  propiedades = {}
  propiedades["T"] = T
  propiedades["P"] = P
  propiedades["v"] = Pr*(gop+grp)*R*T/P/1000
  propiedades["h"] = Tr*(got+grt)*R*T
  propiedades["s"] = R*(Tr*(got+grt)-(go+gr))
  propiedades["cp"] = -R*Tr**2*(gott+grtt)
  propiedades["cv"] = R*(-Tr**2*(gott+grtt)+((gop+grp)-Tr*(gopt+grpt))**2 /
               (gopp+grpp))
  propiedades["w"] = (R*T*1000*(1+2*Pr*grp+Pr**2*grp**2)/(1-Pr**2*grpp+(
    1+Pr*grp-Tr*Pr*grpt)**2/Tr**2/(gott+grtt)))**0.5
  propiedades["alfav"] = (1+Pr*grp-Tr*Pr*grpt)/(1+Pr*grp)/T
  propiedades["kt"] = (1-Pr**2*grpp)/(1+Pr*grp)/P
  propiedades["region"] = 5
  propiedades["x"] = 1
  return propiedades
}

void Region5_cp0(Tr, Pr)
{
  /*
  """Ideal properties for Region 5

  Parameters
  ----------
  Tr : float
    Reduced temperature, [-]
  Pr : float
    Reduced pressure, [-]

  Returns
  -------
  prop : array
    Array with ideal Gibbs energy partial derivatives:

      * g: Ideal Specific Gibbs energy, [kJ/kg]
      * gp: [∂g/∂P]T
      * gpp: [∂²g/∂P²]T
      * gt: [∂g/∂T]P
      * gtt: [∂²g/∂T²]P
      * gpt: [∂²g/∂T∂P]

  References
  ----------
  IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  Thermodynamic Properties of Water and Steam August 2007,
  http://www.iapws.org/relguide/IF97-Rev.html, Eq 33
  """
  */
  Jo = [0, 1, -3, -2, -1, 2]
  no = [-0.13179983674201e2, 0.68540841634434e1, -0.24805148933466e-1, 0.36901534980333, -0.31161318213925e1, -0.32961626538917]
  go = log(Pr)
  gop = Pr**-1
  gopp = -Pr**-2
  got = gott = gopt = 0
  for(int i = 0; i < Jo.length; i++)
  {
    j = Jo[i];
    ni= no[i];
    go += ni * Tr**j
    got += ni*j * Tr**(j-1)
    gott += ni*j*(j-1) * Tr**(j-2)
  }

  return go, gop, gopp, got, gott, gopt
}

// Region definitions
void _Bound_TP(T, P)
{
  /*
  """Region definition for input T and P

  Parameters
  ----------
  T : float
    Temperature, [K]
  P : float
    Pressure, [MPa]

  Returns
  -------
  region : float
    IAPWS-97 region code

  References
  ----------
  Wagner, W; Kretzschmar, H-J: International Steam Tables: Properties of
  Water and Steam Based on the Industrial Formulation IAPWS-IF97; Springer,
  2008; doi: 10.1007/978-3-540-74234-0. Fig. 2.3
  """
  */
  region = None
  if(1073.15 < T <= 2273.15 && Pmin <= P <= 50)
  {
    region = 5
  }
  else if(Pmin <= P <= Ps_623)
  {
    Tsat = _TSat_P(P)
    if(273.15 <= T <= Tsat)
    {
      region = 1
    }
    else if(Tsat < T <= 1073.15)
    {
      region = 2
    }
  }
  else if(Ps_623 < P <= 100)
  {
    T_b23 = _t_P(P)
    if(273.15 <= T <= 623.15)
    {
      region = 1
    }
    else if(623.15 < T < T_b23)
    {
      region = 3
    }
    else if(T_b23 <= T <= 1073.15)
    {
      region = 2
    }
  }
  return region
}

void _Bound_Ph(P, h)
{
  /*
  """Region definition for input P y h

  Parameters
  ----------
  P : float
    Pressure, [MPa]
  h : float
    Specific enthalpy, [kJ/kg]

  Returns
  -------
  region : float
    IAPWS-97 region code

  References
  ----------
  Wagner, W; Kretzschmar, H-J: International Steam Tables: Properties of
  Water and Steam Based on the Industrial Formulation IAPWS-IF97; Springer,
  2008; doi: 10.1007/978-3-540-74234-0. Fig. 2.5
  """
  */

  region = None
  if(Pmin <= P <= Ps_623)
  {
    h14 = _Region1(_TSat_P(P), P)["h"]
    h24 = _Region2(_TSat_P(P), P)["h"]
    h25 = _Region2(1073.15, P)["h"]
    hmin = _Region1(273.15, P)["h"]
    hmax = _Region5(2273.15, P)["h"]
    if(hmin <= h <= h14)
    {
      region = 1
    }
    else if(h14 < h < h24)
    {
      region = 4
    }
    else if(h24 <= h <= h25)
    {
      region = 2
    }
    else if(h25 < h <= hmax)
    {
      region = 5
    }
  }
  else if(Ps_623 < P < Pc)
  {
    hmin = _Region1(273.15, P)["h"]
    h13 = _Region1(623.15, P)["h"]
    h32 = _Region2(_t_P(P), P)["h"]
    h25 = _Region2(1073.15, P)["h"]
    hmax = _Region5(2273.15, P)["h"]
    if(hmin <= h <= h13)
    {
      region = 1
    }
    else if(h13 < h < h32)
    {
      try:
      {
        p34 = _PSat_h(h)
      }
      except NotImplementedError:
      {
        p34 = Pc
      }
      if(P < p34)
      {
        region = 4
      }
      else
      {
        region = 3
      }
    }
    else if(h32 <= h <= h25)
    {
      region = 2
    }
    else if(h25 < h <= hmax)
    {
      region = 5
    }
  }
  else if(Pc <= P <= 100)
  {
    hmin = _Region1(273.15, P)["h"]
    h13 = _Region1(623.15, P)["h"]
    h32 = _Region2(_t_P(P), P)["h"]
    h25 = _Region2(1073.15, P)["h"]
    hmax = _Region5(2273.15, P)["h"]
    if(hmin <= h <= h13)
    {
      region = 1
    }
    else if(h13 < h < h32)
    {
      region = 3
    }
    else if(h32 <= h <= h25)
    {
      region = 2
    }
    else if(P <= 50 && h25 <= h <= hmax)
    {
      region = 5
    }
  }
  return region
}

void _Bound_Ps(P, s)
{
  /*
  """Region definition for input P and s

  Parameters
  ----------
  P : float
    Pressure, [MPa]
  s : float
    Specific entropy, [kJ/kgK]

  Returns
  -------
  region : float
    IAPWS-97 region code

  References
  ----------
  Wagner, W; Kretzschmar, H-J: International Steam Tables: Properties of
  Water and Steam Based on the Industrial Formulation IAPWS-IF97; Springer,
  2008; doi: 10.1007/978-3-540-74234-0. Fig. 2.9
  """
  */
  region = None
  if(Pmin <= P <= Ps_623)
  {
    smin = _Region1(273.15, P)["s"]
    s14 = _Region1(_TSat_P(P), P)["s"]
    s24 = _Region2(_TSat_P(P), P)["s"]
    s25 = _Region2(1073.15, P)["s"]
    smax = _Region5(2273.15, P)["s"]
    if(smin <= s <= s14)
    {
      region = 1
    }
    else if(s14 < s < s24)
    {
      region = 4
    }
    else if(s24 <= s <= s25)
    {
      region = 2
    }
    else if(s25 < s <= smax)
    {
      region = 5
    }
  }
  else if(Ps_623 < P < Pc)
  {
    smin = _Region1(273.15, P)["s"]
    s13 = _Region1(623.15, P)["s"]
    s32 = _Region2(_t_P(P), P)["s"]
    s25 = _Region2(1073.15, P)["s"]
    smax = _Region5(2273.15, P)["s"]
    if(smin <= s <= s13)
    {
      region = 1
    }
    else if(s13 < s < s32)
    {
      try:
      {
        p34 = _PSat_s(s)
      }
      except NotImplementedError:
      {
        p34 = Pc
      }
      if(P < p34)
      {
        region = 4
      }
      else
      {
        region = 3
      }
    }
    else if(s32 <= s <= s25)
    {
      region = 2
    }
    else if(s25 < s <= smax)
    {
      region = 5
    }
  }
  else if(Pc <= P <= 100)
  {
    smin = _Region1(273.15, P)["s"]
    s13 = _Region1(623.15, P)["s"]
    s32 = _Region2(_t_P(P), P)["s"]
    s25 = _Region2(1073.15, P)["s"]
    smax = _Region5(2273.15, P)["s"]
    if(smin <= s <= s13)
    {
      region = 1
    }
    else if(s13 < s < s32)
    {
      region = 3
    }
    else if(s32 <= s <= s25)
    {
      region = 2
    }
    else if(P <= 50 && s25 <= s <= smax)
    {
      region = 5
    }
  }
  return region
}

void _Bound_hs(h, s)
{
  /*
  """Region definition for input h and s

  Parameters
  ----------
  h : float
    Specific enthalpy, [kJ/kg]
  s : float
    Specific entropy, [kJ/kgK]

  Returns
  -------
  region : float
    IAPWS-97 region code

  References
  ----------
  Wagner, W; Kretzschmar, H-J: International Steam Tables: Properties of
  Water and Steam Based on the Industrial Formulation IAPWS-IF97; Springer,
  2008; doi: 10.1007/978-3-540-74234-0. Fig. 2.14
  """
  */
  region = None
  s13 = _Region1(623.15, 100)["s"]
  s13s = _Region1(623.15,  Ps_623)["s"]
  sTPmax = _Region2(1073.15, 100)["s"]
  s2ab = _Region2(1073.15, 4)["s"]

  // Left point in h-s plot
  smin = _Region1(273.15, 100)["s"]
  hmin = _Region1(273.15, Pmin)["h"]

  // Right point in h-s plot
  _Pmax = _Region2(1073.15, Pmin)
  hmax = _Pmax["h"]
  smax = _Pmax["s"]

  // Region 4 left and right point
  _sL = _Region1(273.15, Pmin)
  h4l = _sL["h"]
  s4l = _sL["s"]
  _sV = _Region2(273.15, Pmin)
  h4v = _sV["h"]
  s4v = _sV["s"]

  if(smin <= s <= s13)
  {
    hmin = h4l+(s-s4l)/(s4v-s4l)*(h4v-h4l)
    hs = _h1_s(s)
    T = _Backward1_T_Ps(100, s)-0.0218
    hmax = _Region1(T, 100)["h"]
    if(hmin <= h < hs)
    {
      region = 4
    }
    else if(hs <= h <= hmax)
    {
      region = 1
    }
  }
  else if(s13 < s <= s13s)
  {
    hmin = h4l+(s-s4l)/(s4v-s4l)*(h4v-h4l)
    hs = _h1_s(s)
    h13 = _h13_s(s)
    v = _Backward3_v_Ps(100, s)*(1+9.6e-5)
    T = _Backward3_T_Ps(100, s)-0.0248
    hmax = _Region3(1/v, T)["h"]
    if(hmin <= h < hs)
    {
      region = 4
    }
    else if(hs <= h < h13)
    {
      region = 1
    }
    else if(h13 <= h <= hmax)
    {
      region = 3
    }
  }
  else if(s13s < s <= sc)
  {
    hmin = h4l+(s-s4l)/(s4v-s4l)*(h4v-h4l)
    hs = _h3a_s(s)
    v = _Backward3_v_Ps(100, s)*(1+9.6e-5)
    T = _Backward3_T_Ps(100, s)-0.0248
    hmax = _Region3(1/v, T)["h"]
    if(hmin <= h < hs)
    {
      region = 4
    }
    else if(hs <= h <= hmax)
    {
      region = 3
    }
  }
  else if(sc < s < 5.049096828)
  {
    hmin = h4l+(s-s4l)/(s4v-s4l)*(h4v-h4l)
    hs = _h2c3b_s(s)
    v = _Backward3_v_Ps(100, s)*(1+9.6e-5)
    T = _Backward3_T_Ps(100, s)-0.0248
    hmax = _Region3(1/v, T)["h"]
    if(hmin <= h < hs)
    {
      region = 4
    }
    else if(hs <= h <= hmax)
    {
      region = 3
    }
  }
  else if(5.049096828 <= s < 5.260578707)
  {
    // Specific zone with 2-3 boundary in s shape
    hmin = h4l+(s-s4l)/(s4v-s4l)*(h4v-h4l)
    hs = _h2c3b_s(s)
    h23max = _Region2(863.15, 100)["h"]
    h23min = _Region2(623.15, Ps_623)["h"]
    T = _Backward2_T_Ps(100, s)-0.019
    hmax = _Region2(T, 100)["h"]

    if(hmin <= h < hs)
    {
      region = 4
    }
    else if(hs <= h < h23min)
    {
      region = 3
    }
    else if(h23min <= h < h23max)
    {
      if(_Backward2c_P_hs(h, s) <= _P23_T(_t_hs(h, s)))
      {
        region = 2
      }
      else
      {
        region = 3
      }
    }
    else if(h23max <= h <= hmax)
    {
      region = 2
    }
  }
  else if(5.260578707 <= s < 5.85)
  {
    hmin = h4l+(s-s4l)/(s4v-s4l)*(h4v-h4l)
    hs = _h2c3b_s(s)
    T = _Backward2_T_Ps(100, s)-0.019
    hmax = _Region2(T, 100)["h"]
    if(hmin <= h < hs)
    {
      region = 4
    }
    else if(hs <= h <= hmax)
    {
      region = 2
    }
  }
  else if(5.85 <= s < sTPmax)
  {
    hmin = h4l+(s-s4l)/(s4v-s4l)*(h4v-h4l)
    hs = _h2ab_s(s)
    T = _Backward2_T_Ps(100, s)-0.019
    hmax = _Region2(T, 100)["h"]
    if(hmin <= h < hs)
    {
      region = 4
    }
    else if(hs <= h <= hmax)
    {
      region = 2
    }
  }
  else if(sTPmax <= s < s2ab)
  {
    hmin = h4l+(s-s4l)/(s4v-s4l)*(h4v-h4l)
    hs = _h2ab_s(s)
    P = _Backward2_P_hs(h, s)
    hmax = _Region2(1073.15, P)["h"]
    if(hmin <= h < hs)
    {
      region = 4
    }
    else if(hs <= h <= hmax)
    {
      region = 2
    }
  }
  else if(s2ab <= s < s4v)
  {
    hmin = h4l+(s-s4l)/(s4v-s4l)*(h4v-h4l)
    hs = _h2ab_s(s)
    P = _Backward2_P_hs(h, s)
    hmax = _Region2(1073.15, P)["h"]
    if(hmin <= h < hs)
    {
      region = 4
    }
    else if(hs <= h <= hmax)
    {
      region = 2
    }
  }
  else if(s4v <= s <= smax)
  {
    hmin = _Region2(273.15, Pmin)["h"]
    P = _Backward2a_P_hs(h, s)
    hmax = _Region2(1073.15, P)["h"]
    if(Pmin <= P <= 100 && hmin <= h <= hmax)
    {
      region = 2
    }
  }
  // Check region 5
  if(!region && _Region5(1073.15, 50)["s"] < s <= _Region5(2273.15, Pmin)["s"] && _Region5(1073.15, 50)["h"] < h <= _Region5(2273.15, Pmin)["h"])
  {
    void funcion(par)
    {
      return (_Region5(par[0], par[1])["h"]-h, _Region5(par[0], par[1])["s"]-s)
    }

    T, P = fsolve(funcion, [1400, 1])
    if(1073.15 < T <= 2273.15 && Pmin <= P <= 50)
    {
      region = 5
    }
  }
  return region
}

void prop0(T, P)
{
  /*
  """Ideal gas properties

  Parameters
  ----------
  T : float
    Temperature, [K]
  P : float
    Pressure, [MPa]

  Returns
  -------
  prop : dict
    Dict with calculated properties. The available properties are:

      * v: Specific volume, [m³/kg]
      * h: Specific enthalpy, [kJ/kg]
      * s: Specific entropy, [kJ/kgK]
      * cp: Specific isobaric heat capacity, [kJ/kgK]
      * cv: Specific isocoric heat capacity, [kJ/kgK]
      * w: Speed of sound, [m/s]
      * alfav: Cubic expansion coefficient, [1/K]
      * kt: Isothermal compressibility, [1/MPa]
  """
  */
  if(T <= 1073.15)
  {
    Tr = 540/T
    Pr = P/1.
    go, gop, gopp, got, gott, gopt = Region2_cp0(Tr, Pr)
  }
  else
  {
    Tr = 1000/T
    Pr = P/1.
    go, gop, gopp, got, gott, gopt = Region5_cp0(Tr, Pr)
  }
  prop0 = {}
  prop0["v"] = Pr*gop*R*T/P/1000
  prop0["h"] = Tr*got*R*T
  prop0["s"] = R*(Tr*got-go)
  prop0["cp"] = -R*Tr**2*gott
  prop0["cv"] = R*(-Tr**2*gott-1)

  prop0["w"] = (R*T*1000/(1+1/Tr**2/gott))**0.5
  prop0["alfav"] = 1/T
  prop0["xkappa"] = 1/P
  return prop0
}

class IAPWS97(object):
{
  /*
  """Class to model a state of liquid water or steam with the IAPWS-IF97

  Parameters
  ----------
  T : float
    Temperature, [K]
  P : float
    Pressure, [MPa]
  h : float
    Specific enthalpy, [kJ/kg]
  s : float
    Specific entropy, [kJ/kgK]
  x : float
    Vapor quality, [-]
  l : float, optional
    Wavelength of light, for refractive index, [nm]

  Notes
  -----
  Definitions options:

    * T, P: Not valid for two-phases region
    * P, h
    * P, s
    * h, s
    * T, x: Only for two-phases region
    * P, x: Only for two-phases region

  Returns
  -------
  prop : dict
    The calculated instance has the following properties:

      * P: Pressure, [MPa]
      * T: Temperature, [K]
      * g: Specific Gibbs free energy, [kJ/kg]
      * a: Specific Helmholtz free energy, [kJ/kg]
      * v: Specific volume, [m³/kg]
      * rho: Density, [kg/m³]
      * h: Specific enthalpy, [kJ/kg]
      * u: Specific internal energy, [kJ/kg]
      * s: Specific entropy, [kJ/kg·K]
      * cp: Specific isobaric heat capacity, [kJ/kg·K]
      * cv: Specific isochoric heat capacity, [kJ/kg·K]
      * Z: Compression factor, [-]
      * fi: Fugacity coefficient, [-]
      * f: Fugacity, [MPa]

      * gamma: Isoentropic exponent, [-]
      * alfav: Isobaric cubic expansion coefficient, [1/K]
      * xkappa: Isothermal compressibility, [1/MPa]
      * kappas: Adiabatic compresibility, [1/MPa]
      * alfap: Relative pressure coefficient, [1/K]
      * betap: Isothermal stress coefficient, [kg/m³]
      * joule: Joule-Thomson coefficient, [K/MPa]
      * deltat: Isothermal throttling coefficient, [kJ/kg·MPa]
      * region: Region

      * v0: Ideal specific volume, [m³/kg]
      * u0: Ideal specific internal energy, [kJ/kg]
      * h0: Ideal specific enthalpy, [kJ/kg]
      * s0: Ideal specific entropy, [kJ/kg·K]
      * a0: Ideal specific Helmholtz free energy, [kJ/kg]
      * g0: Ideal specific Gibbs free energy, [kJ/kg]
      * cp0: Ideal specific isobaric heat capacity, [kJ/kg·K]
      * cv0: Ideal specific isochoric heat capacity [kJ/kg·K]
      * w0: Ideal speed of sound, [m/s]
      * gamma0: Ideal isoentropic exponent, [-]

      * w: Speed of sound, [m/s]
      * mu: Dynamic viscosity, [Pa·s]
      * nu: Kinematic viscosity, [m²/s]
      * k: Thermal conductivity, [W/m·K]
      * alfa: Thermal diffusivity, [m²/s]
      * sigma: Surface tension, [N/m]
      * epsilon: Dielectric constant, [-]
      * n: Refractive index, [-]
      * Prandt: Prandtl number, [-]
      * Pr: Reduced Pressure, [-]
      * Tr: Reduced Temperature, [-]
      * Hvap: Vaporization heat, [kJ/kg]
      * Svap: Vaporization entropy, [kJ/kg·K]

  Examples
  --------
  >>> water=IAPWS97(T=170+273.15, x=0.5)
  >>> water.Liquid.cp, water.Vapor.cp, water.Liquid.w, water.Vapor.w
  4.3695 2.5985 1418.3 498.78

  >>> water=IAPWS97(T=325+273.15, x=0.5)
  >>> water.P, water.Liquid.v, water.Vapor.v, water.Liquid.h, water.Vapor.h
  12.0505 0.00152830 0.0141887 1493.37 2684.48

  >>> water=IAPWS97(T=50+273.15, P=0.0006112127)
  >>> water.cp0, water.cv0, water.h0, water.s0, water.w0
  1.8714 1.4098 2594.66 9.471 444.93
  """
  */

  kwargs =
  {
    "T": 0.0,
    "P": 0.0,
    "x": None,
    "h": None,
    "s": None,
    "v": 0.0,
    "l": 0.5893
  }
  status = 0
  msg = "Unknown variables"

  void __init__(self, **kwargs)
  {
    self.kwargs = IAPWS97.kwargs.copy()
    self.__call__(**kwargs)
  }

  void __call__(self, **kwargs)
  {
    self.kwargs.update(kwargs)

    if(self.calculable)
    {
      self.status = 1
      self.calculo()
      self.msg = "Solved"
    }
  }

  @property
  void calculable(self)
  {
    //Check if class is calculable by its kwargs
    self._thermo = ""
    if(self.kwargs["T"] && self.kwargs["P"])
    {
      self._thermo = "TP"
    }
    else if(self.kwargs["P"] && self.kwargs["h"] is !None)
    {
      self._thermo = "Ph"
    }
    else if(self.kwargs["P"] && self.kwargs["s"] is !None)
    {
      self._thermo = "Ps"
    }
    // TODO: Add other pairs definitions options
    // else if(self.kwargs["P"] && self.kwargs["v"])
    // {
      // self._thermo = "Pv"
    // }
    // else if(self.kwargs["T"] && self.kwargs["s"] is !None)
    // {
      // self._thermo = "Ts"
    // }
    else if(self.kwargs["h"] is !None && self.kwargs["s"] is !None)
    {
      self._thermo = "hs"
    }
    else if(self.kwargs["T"] && self.kwargs["x"] is !None)
    {
      self._thermo = "Tx"
    }
    else if(self.kwargs["P"] && self.kwargs["x"] is !None)
    {
      self._thermo = "Px"
    }
    return self._thermo
  }

  void calculo(self)
  {
    propiedades = None
    args = (self.kwargs[self._thermo[0]], self.kwargs[self._thermo[1]])
    if(self._thermo == "TP")
    {
      T, P = args
      region = _Bound_TP(T, P)
      if(region == 1)
      {
        propiedades = _Region1(T, P)
      }
      else if(region == 2)
      {
        propiedades = _Region2(T, P)
      }
      else if(region == 3)
      {
        if(T == Tc && P == Pc)
        {
          rho = rhoc
        }
        else
        {
          vo = _Backward3_v_PT(P, T)

          void funcion(rho)
          {
            return _Region3(rho, self.kwargs["T"])["P"]-P
          }
          rho = newton(funcion, 1/vo)
        }
        propiedades = _Region3(rho, T)
      }
      else if(region == 5)
      {
        propiedades = _Region5(T, P)
      }
      else
      {
        raise NotImplementedError("Incoming out of bound")
      }
    }
    else if(self._thermo == "Ph")
    {
      P, h = args
      region = _Bound_Ph(P, h)
      if(region == 1)
      {
        To = _Backward1_T_Ph(P, h)
        T = newton(lambda T: _Region1(T, P)["h"]-h, To)
        propiedades = _Region1(T, P)
      }
      else if(region == 2)
      {
        To = _Backward2_T_Ph(P, h)
        T = newton(lambda T: _Region2(T, P)["h"]-h, To)
        propiedades = _Region2(T, P)
      }
      else if(region == 3)
      {
        vo = _Backward3_v_Ph(P, h)
        To = _Backward3_T_Ph(P, h)

        void funcion(par)
        {
          return (_Region3(par[0], par[1])["h"]-h, _Region3(par[0], par[1])["P"]-P)
        }
        rho, T = fsolve(funcion, [1/vo, To])
        propiedades = _Region3(rho, T)
      }
      else if(region == 4)
      {
        T = _TSat_P(P)
        if(T <= 623.15)
        {
          h1 = _Region1(T, P)["h"]
          h2 = _Region2(T, P)["h"]
          x = (h-h1)/(h2-h1)
          propiedades = _Region4(P, x)
        }
        else
        {
          vo = _Backward3_v_Ph(P, h)
          To = _Backward3_T_Ph(P, h)

          void funcion(par)
          {
            return (_Region3(par[0], par[1])["h"]-h, _Region3(par[0], par[1])["P"]-P)
          }
          rho, T = fsolve(funcion, [1/vo, To])
          propiedades = _Region3(rho, T)
        }
      }
      else if(region == 5)
      {
        T = newton(lambda T: _Region5(T, P)["h"]-h, 1500)
        propiedades = _Region5(T, P)
      }
      else
      {
        raise NotImplementedError("Incoming out of bound")
      }
    }
    else if(self._thermo == "Ps")
    {
      P, s = args
      region = _Bound_Ps(P, s)
      if(region == 1)
      {
        To = _Backward1_T_Ps(P, s)
        T = newton(lambda T: _Region1(T, P)["s"]-s, To)
        propiedades = _Region1(T, P)
      }
      else if(region == 2)
      {
        To = _Backward2_T_Ps(P, s)
        T = newton(lambda T: _Region2(T, P)["s"]-s, To)
        propiedades = _Region2(T, P)
      }
      else if(region == 3)
      {
        vo = _Backward3_v_Ps(P, s)
        To = _Backward3_T_Ps(P, s)

        void funcion(par)
        {
          return (_Region3(par[0], par[1])["s"]-s, _Region3(par[0], par[1])["P"]-P)
        }

        rho, T = fsolve(funcion, [1/vo, To])
        propiedades = _Region3(rho, T)
      }
      else if(region == 4)
      {
        T = _TSat_P(P)
        if(T <= 623.15)
        {
          s1 = _Region1(T, P)["s"]
          s2 = _Region2(T, P)["s"]
          x = (s-s1)/(s2-s1)
          propiedades = _Region4(P, x)
        }
        else
        {
          vo = _Backward3_v_Ps(P, s)
          To = _Backward3_T_Ps(P, s)

          void funcion(par)
          {
            return (_Region3(par[0], par[1])["s"]-s, _Region3(par[0], par[1])["P"]-P)
          }
          rho, T = fsolve(funcion, [1/vo, To])
          propiedades = _Region3(rho, T)
        }
      }
      else if(region == 5)
      {
        T = newton(lambda T: _Region5(T, P)["s"]-s, 1500)
        propiedades = _Region5(T, P)
      }
      else
      {
        raise NotImplementedError("Incoming out of bound")
      }
    }
    else if(self._thermo == "hs")
    {
      h, s = args
      region = _Bound_hs(h, s)
      if(region == 1)
      {
        Po = _Backward1_P_hs(h, s)
        To = _Backward1_T_Ph(Po, h)

        void funcion(par)
        {
          return (_Region1(par[0], par[1])["h"]-h, _Region1(par[0], par[1])["s"]-s)
        }
        T, P = fsolve(funcion, [To, Po])
        propiedades = _Region1(T, P)
      }
      else if(region == 2)
      {
        Po = _Backward2_P_hs(h, s)
        To = _Backward2_T_Ph(Po, h)

        void funcion(par)
        {
          return (_Region2(par[0], par[1])["h"]-h, _Region2(par[0], par[1])["s"]-s)
        }
        T, P = fsolve(funcion, [To, Po])
        propiedades = _Region2(T, P)
      }
      else if(region == 3)
      {
        P = _Backward3_P_hs(h, s)
        vo = _Backward3_v_Ph(P, h)
        To = _Backward3_T_Ph(P, h)

        void funcion(par)
        {
          return (_Region3(par[0], par[1])["h"]-h, _Region3(par[0], par[1])["s"]-s)
        }
        rho, T = fsolve(funcion, [1/vo, To])
        propiedades = _Region3(rho, T)
      }
      else if(region == 4)
      {
        if(round(s-sc, 6) == 0 && round(h-hc, 6) == 0)
        {
          propiedades = _Region3(rhoc, Tc)
        }
        else
        {
          To = _Backward4_T_hs(h, s)
          if(To < 273.15 || To > Tc)
          {
            To = 300
          }
          void funcion(par)
          {
            if(par[1] < 0)
            {
              par[1] = 0
            }
            else if(par[1] > 1)
            {
              par[1] = 1
            }
            if(par[0] < 273.15)
            {
              par[0] = 273.15
            }
            else if(par[0] > Tc)
            {
              par[0] = Tc
            }

            Po = _PSat_T(par[0])
            liquid = _Region1(par[0], Po)
            vapor = _Region2(par[0], Po)
            hl = liquid["h"]
            sl = liquid["s"]
            hv = vapor["h"]
            sv = vapor["s"]
            return (hv*par[1]+hl*(1-par[1])-h, sv*par[1]+sl*(1-par[1])-s)
          }
          T, x = fsolve(funcion, [To, 0.5])
          P = _PSat_T(T)

          if(Pt <= P < Pc && 0 < x < 1)
          {
            propiedades = _Region4(P, x)
          }
          else if(Pt <= P <= Ps_623 && x == 0)
          {
            propiedades = _Region1(T, P)
          }
      }
      else if(region == 5)
      {
        void funcion(par)
        {
          return (_Region5(par[0], par[1])["h"]-h, _Region5(par[0], par[1])["s"]-s)
        }
        T, P = fsolve(funcion, [1400, 1])
        propiedades = _Region5(T, P)
      }
      else
      {
        raise NotImplementedError("Incoming out of bound")
      }

    }
    else if(self._thermo == "Px")
    {
      P, x = args
      T = _TSat_P(P)
      if(Pt <= P < Pc && 0 < x < 1)
      {
        propiedades = _Region4(P, x)
      }
      else if(Pt <= P <= Ps_623 && x == 0)
      {
        propiedades = _Region1(T, P)
      }
      else if(Pt <= P <= Ps_623 && x == 1)
      {
        propiedades = _Region2(T, P)
      }
      else if(Ps_623 < P < Pc && x in (0, 1))
      {
        void funcion(rho)
        {
          return _Region3(rho, T)["P"]-P
        }
        rhoo = 1./_Backward3_sat_v_P(P, T, x)
        rho = fsolve(funcion, rhoo)[0]
        propiedades = _Region3(rho, T)
      }
      else if(P == Pc && 0 <= x <= 1)
      {
        propiedades = _Region3(rhoc, Tc)
      }
      else
      {
        raise NotImplementedError("Incoming out of bound")
      }
      self.sigma = _Tension(T)
      propiedades["x"] = x
    }
    else if(self._thermo == "Tx")
    {
      T, x = args
      P = _PSat_T(T)
      if(273.15 <= T < Tc && 0 < x < 1)
      {
        propiedades = _Region4(P, x)
      }
      else if(273.15 <= T <= 623.15 && x == 0)
      {
        propiedades = _Region1(T, P)
      }
      else if(273.15 <= T <= 623.15 && x == 1)
      {
        propiedades = _Region2(T, P)
      }
      else if(623.15 < T < Tc && x in (0, 1))
      {
        rho = 1./_Backward3_sat_v_P(P, T, x)
        propiedades = _Region3(rho, T)
      }
      else if(T == Tc && 0 <= x <= 1)
      {
        propiedades = _Region3(rhoc, Tc)
      }
      else
      {
        raise NotImplementedError("Incoming out of bound")
      }
      self.sigma = _Tension(T)
      propiedades["x"] = x
    }

    self.M = 18.015257  // kg/kmol
    self.Pc = Pc
    self.Tc = Tc
    self.rhoc = rhoc
    self.Tt = Tt
    self.Tb = Tb
    self.f_accent = f_acent
    self.dipole = Dipole

    self.x = propiedades["x"]
    self.region = propiedades["region"]
    self.name = "water"
    self.synonim = "R-718"
    self.CAS = "7732-18-5"

    self.T = propiedades["T"]
    self.P = propiedades["P"]
    self.v = propiedades["v"]
    self.rho = 1/self.v
    self.phase = getphase(self.Tc, self.Pc, self.T, self.P, self.x, self.region)
    self.Tr = self.T/self.Tc
    self.Pr = self.P/self.Pc

    // Ideal properties
    if(self.region in [2, 5])
    {
      cp0 = prop0(self.T, self.P)
      self.v0 = cp0["v"]
      self.h0 = cp0["h"]
      self.u0 = self.h0-self.P*1000*self.v0
      self.s0 = cp0["s"]
      self.a0 = self.u0-self.T*self.s0
      self.g0 = self.h0-self.T*self.s0

      self.cp0 = cp0["cp"]
      self.cv0 = cp0["cv"]
      self.cp0_cv = self.cp0/self.cv0
      self.w0 = cp0["w"]
      self.gamma0 = self.cp0_cv
    }
    else
    {
      self.v0 = None
      self.h0 = None
      self.u0 = None
      self.s0 = None
      self.a0 = None
      self.g0 = 0

      self.cp0 = None
      self.cv0 = None
      self.cp0_cv = None
      self.w0 = None
      self.gamma0 = None
    }

    self.Liquid = _fase()
    self.Vapor = _fase()
    if(self.x == 0)
    {
      // only liquid phase
      self.fill(self, propiedades)
      self.fill(self.Liquid, propiedades)
      self.sigma = _Tension(self.T)
    }
    else if(self.x == 1)
    {
      // only vapor phase
      self.fill(self, propiedades)
      self.fill(self.Vapor, propiedades)
    }
    else
    {
      // two phases
      liquido = _Region1(self.T, self.P)
      self.fill(self.Liquid, liquido)
      vapor = _Region2(self.T, self.P)
      self.fill(self.Vapor, vapor)

      self.h = propiedades["h"]
      self.u = self.h-self.P*1000*self.v
      self.s = propiedades["s"]
      self.a = self.u-self.T*self.s
      self.g = self.h-self.T*self.s
      self.sigma = _Tension(self.T)

      self.Hvap = vapor["h"]-liquido["h"]
      self.Svap = vapor["s"]-liquido["s"]
    }
  }

  void fill(self, fase, estado)
  {
    fase.v = estado["v"]
    fase.rho = 1/fase.v

    fase.h = estado["h"]
    fase.s = estado["s"]
    fase.u = fase.h-self.P*1000*fase.v
    fase.a = fase.u-self.T*fase.s
    fase.g = fase.h-self.T*fase.s

    fase.cv = estado["cv"]
    fase.cp = estado["cp"]
    fase.cp_cv = fase.cp/fase.cv
    fase.w = estado["w"]

    fase.Z = self.P*fase.v/R*1000/self.T
    fase.alfav = estado["alfav"]
    fase.xkappa = estado["kt"]
    fase.kappas = -1/fase.v*self.derivative("v", "P", "s", fase)

    fase.joule = self.derivative("T", "P", "h", fase)
    fase.deltat = self.derivative("h", "P", "T", fase)
    fase.gamma = -fase.v/self.P*self.derivative("P", "v", "s", fase)

    fase.alfap = fase.alfav/self.P/fase.xkappa
    fase.betap = -1/self.P*self.derivative("P", "v", "T", fase)

    fase.fi = exp((fase.g-self.g0)/R/self.T)
    fase.f = self.P*fase.fi

    fase.mu = _Viscosity(fase.rho, self.T)
    // Use industrial formulation for critical enhancement in thermal
    // conductivity calculation
    fase.drhodP_T = self.derivative("rho", "P", "T", fase)
    fase.k = _ThCond(fase.rho, self.T, fase)

    fase.nu = fase.mu/fase.rho
    fase.alfa = fase.k/1000/fase.rho/fase.cp
    try:
    {
      fase.epsilon = _Dielectric(fase.rho, self.T)
    }
    except NotImplementedError:
    {
      fase.epsilon = None
    }
    fase.Prandt = fase.mu*fase.cp*1000/fase.k
    try:
    {
      fase.n = _Refractive(fase.rho, self.T, self.kwargs["l"])
    }
    except NotImplementedError:
    {
      fase.n = None
    }
  }

  void derivative(self, z, x, y, fase)
  {
    /*
    """Wrapper derivative for custom derived properties
    where x, y, z can be: P, T, v, u, h, s, g, a"""
    */
    return deriv_G(self, z, x, y, fase)
  }
}

class IAPWS97_PT(IAPWS97):
{
  //Derivated class for direct P and T input
  void __init__(self, P, T)
  {
    IAPWS97.__init__(self, T=T, P=P)
  }
}

class IAPWS97_Ph(IAPWS97):
{
  //Derivated class for direct P and h input
  void __init__(self, P, h)
  {
    IAPWS97.__init__(self, P=P, h=h)
  }
}

class IAPWS97_Ps(IAPWS97):
{
  //Derivated class for direct P and s input
  void __init__(self, P, s)
  {
    IAPWS97.__init__(self, P=P, s=s)
  }
}

class IAPWS97_Px(IAPWS97):
{
  //Derivated class for direct P and x input
  void __init__(self, P, x)
  {
    IAPWS97.__init__(self, P=P, x=x)
  }
}

class IAPWS97_Tx(IAPWS97):
{
  //Derivated class for direct T and x input
  void __init__(self, T, x)
  {
    IAPWS97.__init__(self, T=T, x=x)
  }
}

//!/usr/bin/python
// -*- coding: utf-8 -*-

/*
"""
Implemented multiparameter equation of state as a Helmholtz free energy:

  * :class:`MEoS`: Base class of multiparameter equation of state
  * :class:`IAPWS95`: 2016 revision of 1995 formulation for ordinaty water
  * :class:`D2O`: 2017 formulation for heavy water.
"""
*/


from __future__ import division
import os
import platform
import warnings

from scipy import exp, log, ndarray
from scipy.optimize import fsolve

from .iapws97 import _TSat_P, IAPWS97
from ._iapws import M, Tc, Pc, rhoc, Tc_D2O, Pc_D2O, rhoc_D2O
from ._iapws import _Viscosity, _ThCond, _Dielectric, _Refractive, _Tension
from ._iapws import _D2O_Viscosity, _D2O_ThCond, _D2O_Tension
from ._utils import _fase, getphase, deriv_H


void _phir(tau, delta, coef)
{
  /*
  """Residual contribution to the adimensional free Helmholtz energy

  Parameters
  ----------
  tau : float
    Inverse reduced temperature Tc/T, [-]
  delta : float
    Reduced density rho/rhoc, [-]
  coef : dict
    Dictionary with equation of state parameters

  Returns
  -------
  fir : float
    Adimensional free Helmholtz energy

  References
  ----------
  IAPWS, Revised Release on the IAPWS Formulation 1995 for the
  Thermodynamic Properties of Ordinary Water Substance for General and
  Scientific Use, September 2016, Table 5
  http://www.iapws.org/relguide/IAPWS-95.html
  """
  */
  fir = 0

  // Polinomial terms
  nr1 = coef.get("nr1", [])
  d1 = coef.get("d1", [])
  t1 = coef.get("t1", [])
  for(int i = 0; i < nr1.length; i++)
  {
    n = nr1[i];
    d = d1[i];
    t = t1[i];
    fir += n*delta**d*tau**t
  }

  // Exponential terms
  nr2 = coef.get("nr2", [])
  d2 = coef.get("d2", [])
  g2 = coef.get("gamma2", [])
  t2 = coef.get("t2", [])
  c2 = coef.get("c2", [])
  for(int iter = 0; iter < nr2.length; iter++)
  {
    n = nr2[iter];
    d = d2[iter];
    g = g2[iter];
    t = t2[iter];
    c = c2[iter];
    fir += n*delta**d*tau**t*exp(-g*delta**c)
  }

  // Gaussian terms
  nr3 = coef.get("nr3", [])
  d3 = coef.get("d3", [])
  t3 = coef.get("t3", [])
  a3 = coef.get("alfa3", [])
  e3 = coef.get("epsilon3", [])
  b3 = coef.get("beta3", [])
  g3 = coef.get("gamma3", [])
  for(int iter = 0; iter < nr3.length; iter++)
  {
    n = nr3[iter];
    d = d3[iter];
    t = t3[iter];
    a = a3[iter];
    e = e3[iter];
    b = b3[iter];
    g = g3[iter];
    fir += n*delta**d*tau**t*exp(-a*(delta-e)**2-b*(tau-g)**2)
  }

  // Non analitic terms
  nr4 = coef.get("nr4", [])
  a4 = coef.get("a4", [])
  b4 = coef.get("b4", [])
  Ai = coef.get("A", [])
  Bi = coef.get("B", [])
  Ci = coef.get("C", [])
  Di = coef.get("D", [])
  bt4 = coef.get("beta4", [])
  for(int iter = 0; iter < nr4.length; iter++)
  {
    n = nr4[iter];
    a = a4[iter];
    b = b4[iter];
    A = Ai[iter];
    B = Bi[iter];
    C = Ci[iter];
    D = Di[iter];
    bt = bt4[iter];
    Tita = (1-tau)+A*((delta-1)**2)**(0.5/bt)
    F = exp(-C*(delta-1)**2-D*(tau-1)**2)
    Delta = Tita**2+B*((delta-1)**2)**a
    fir += n*Delta**b*delta*F
  }

  return fir
}

void _phird(tau, delta, coef)
{
  /*
  """Residual contribution to the adimensional free Helmholtz energy, delta
  derivative

  Parameters
  ----------
  tau : float
    Inverse reduced temperature Tc/T, [-]
  delta : float
    Reduced density rho/rhoc, [-]
  coef : dict
    Dictionary with equation of state parameters

  Returns
  -------
  fird : float
    .. math::
      \left.\frac{\partial \phi^r_{\delta}}{\partial \delta}\right|_{\tau}

  References
  ----------
  IAPWS, Revised Release on the IAPWS Formulation 1995 for the
  Thermodynamic Properties of Ordinary Water Substance for General and
  Scientific Use, September 2016, Table 5
  http://www.iapws.org/relguide/IAPWS-95.html
  """
  */
  fird = 0

  // Polinomial terms
  nr1 = coef.get("nr1", [])
  d1 = coef.get("d1", [])
  t1 = coef.get("t1", [])
  for(int i = 0; i < nr1.length; i++)
  {
    n = nr1[i];
    d = d1[i];
    t = t1[i];
    fird += n*d*delta**(d-1)*tau**t
  }

  // Exponential terms
  nr2 = coef.get("nr2", [])
  d2 = coef.get("d2", [])
  g2 = coef.get("gamma2", [])
  t2 = coef.get("t2", [])
  c2 = coef.get("c2", [])
  for(int iter = 0; iter < nr2.length; iter++)
  {
    n = nr2[iter];
    d = d2[iter];
    g = g2[iter];
    t = t2[iter];
    c = c2[iter];
    fird += n*exp(-g*delta**c)*delta**(d-1)*tau**t*(d-g*c*delta**c)
  }

  // Gaussian terms
  nr3 = coef.get("nr3", [])
  d3 = coef.get("d3", [])
  t3 = coef.get("t3", [])
  a3 = coef.get("alfa3", [])
  e3 = coef.get("epsilon3", [])
  b3 = coef.get("beta3", [])
  g3 = coef.get("gamma3", [])
  for(int iter = 0; iter < nr3.length; iter++)
  {
    n = nr3[iter];
    d = d3[iter];
    t = t3[iter];
    a = a3[iter];
    e = e3[iter];
    b = b3[iter];
    g = g3[iter];
    fird += n*delta**d*tau**t*exp(-a*(delta-e)**2-b*(tau-g)**2)*(
        d/delta-2*a*(delta-e))
  }

  // Non analitic terms
  nr4 = coef.get("nr4", [])
  a4 = coef.get("a4", [])
  b4 = coef.get("b4", [])
  Ai = coef.get("A", [])
  Bi = coef.get("B", [])
  Ci = coef.get("C", [])
  Di = coef.get("D", [])
  bt4 = coef.get("beta4", [])
  for(int iter = 0; iter < nr4.length; iter++)
  {
    n = nr4[iter];
    a = a4[iter];
    b = b4[iter];
    A = Ai[iter];
    B = Bi[iter];
    C = Ci[iter];
    D = Di[iter];
    bt = bt4[iter];

    Tita = (1-tau)+A*((delta-1)**2)**(0.5/bt)
    F = exp(-C*(delta-1)**2-D*(tau-1)**2)
    Fd = -2*C*F*(delta-1)

    Delta = Tita**2+B*((delta-1)**2)**a
    Deltad = (delta-1)*(A*Tita*2/bt*((delta-1)**2)**(0.5/bt-1) +
              2*B*a*((delta-1)**2)**(a-1))
    DeltaBd = b*Delta**(b-1)*Deltad

    fird += n*(Delta**b*(F+delta*Fd)+DeltaBd*delta*F)
  }

  return fird
}

void _phirt(tau, delta, coef)
{
  /*
  """Residual contribution to the adimensional free Helmholtz energy, tau
  derivative

  Parameters
  ----------
  tau : float
    Inverse reduced temperature Tc/T, [-]
  delta : float
    Reduced density rho/rhoc, [-]
  coef : dict
    Dictionary with equation of state parameters

  Returns
  -------
  firt : float
    .. math::
      \left.\frac{\partial \phi^r_{\tau}}{\partial \tau}\right|_{\delta}

  References
  ----------
  IAPWS, Revised Release on the IAPWS Formulation 1995 for the
  Thermodynamic Properties of Ordinary Water Substance for General and
  Scientific Use, September 2016, Table 5
  http://www.iapws.org/relguide/IAPWS-95.html
  """
  */
  firt = 0

  // Polinomial terms
  nr1 = coef.get("nr1", [])
  d1 = coef.get("d1", [])
  t1 = coef.get("t1", [])
  for(int i = 0; i < nr1.length; i++)
  {
    n = nr1[i];
    d = d1[i];
    t = t1[i];
    firt += n*t*delta**d*tau**(t-1)
  }

  // Exponential terms
  nr2 = coef.get("nr2", [])
  d2 = coef.get("d2", [])
  g2 = coef.get("gamma2", [])
  t2 = coef.get("t2", [])
  c2 = coef.get("c2", [])
  for(int iter = 0; iter < nr2.length; iter++)
  {
    n = nr2[iter];
    d = d2[iter];
    g = g2[iter];
    t = t2[iter];
    c = c2[iter];
    firt += n*t*delta**d*tau**(t-1)*exp(-g*delta**c)
  }

  // Gaussian terms
  nr3 = coef.get("nr3", [])
  d3 = coef.get("d3", [])
  t3 = coef.get("t3", [])
  a3 = coef.get("alfa3", [])
  e3 = coef.get("epsilon3", [])
  b3 = coef.get("beta3", [])
  g3 = coef.get("gamma3", [])
  for(int iter = 0; iter < nr3.length; iter++)
  {
    n = nr3[iter];
    d = d3[iter];
    t = t3[iter];
    a = a3[iter];
    e = e3[iter];
    b = b3[iter];
    g = g3[iter];
    firt += n*delta**d*tau**t*exp(-a*(delta-e)**2-b*(tau-g)**2)*(
        t/tau-2*b*(tau-g))
  }

  // Non analitic terms
  nr4 = coef.get("nr4", [])
  a4 = coef.get("a4", [])
  b4 = coef.get("b4", [])
  Ai = coef.get("A", [])
  Bi = coef.get("B", [])
  Ci = coef.get("C", [])
  Di = coef.get("D", [])
  bt4 = coef.get("beta4", [])
  for(int iter = 0; iter < nr4.length; iter++)
  {
    n = nr4[iter];
    a = a4[iter];
    b = b4[iter];
    A = Ai[iter];
    B = Bi[iter];
    C = Ci[iter];
    D = Di[iter];
    bt = bt4[iter];
    Tita = (1-tau)+A*((delta-1)**2)**(0.5/bt)
    F = exp(-C*(delta-1)**2-D*(tau-1)**2)
    Ft = -2*D*F*(tau-1)
    Delta = Tita**2+B*((delta-1)**2)**a
    DeltaBt = -2*Tita*b*Delta**(b-1)
    firt += n*delta*(DeltaBt*F+Delta**b*Ft)
  }

  return firt
}


class MEoS(_fase):
{
  /*
  """
  General implementation of multiparameter equation of state. From this
  derived all child class specified per individual compounds

  Parameters
  ----------
  T : float
    Temperature, [K]
  P : float
    Pressure, [MPa]
  rho : float
    Density, [kg/m³]
  v : float
    Specific volume, [m³/kg]
  h : float
    Specific enthalpy, [kJ/kg]
  s : float
    Specific entropy, [kJ/kgK]
  u : float
    Specific internal energy, [kJ/kg]
  x : float
    Vapor quality, [-]

  l : float, optional
    Wavelength of light, for refractive index, [nm]
  rho0 : float, optional
    Initial value of density, to improve iteration, [kg/m³]
  T0 : float, optional
    Initial value of temperature, to improve iteration, [K]
  x0 : Initial value of vapor quality, necessary in bad input pair definition
    where there are two valid solution (T-h, T-s)

  Notes
  -----
  * It needs two incoming properties of T, P, rho, h, s, u.
  * v as a alternate input parameter to rho
  * T-x, P-x, preferred input pair to specified a point in two phases region

  The calculated instance has the following properties:

    * P: Pressure, [MPa]
    * T: Temperature, [K]
    * x: Vapor quality, [-]
    * g: Specific Gibbs free energy, [kJ/kg]
    * a: Specific Helmholtz free energy, [kJ/kg]
    * v: Specific volume, [m³/kg]
    * r: Density, [kg/m³]
    * h: Specific enthalpy, [kJ/kg]
    * u: Specific internal energy, [kJ/kg]
    * s: Specific entropy, [kJ/kg·K]
    * cp: Specific isobaric heat capacity, [kJ/kg·K]
    * cv: Specific isochoric heat capacity, [kJ/kg·K]
    * cp_cv: Heat capacity ratio, [-]
    * Z: Compression factor, [-]
    * fi: Fugacity coefficient, [-]
    * f: Fugacity, [MPa]
    * gamma: Isoentropic exponent, [-]

    * alfav: Isobaric cubic expansion coefficient, [1/K]
    * kappa: Isothermal compressibility, [1/MPa]
    * kappas: Adiabatic compresibility, [1/MPa]
    * alfap: Relative pressure coefficient, [1/K]
    * betap: Isothermal stress coefficient, [kg/m³]
    * joule: Joule-Thomson coefficient, [K/MPa]

    * betas: Isoentropic temperature-pressure coefficient, [-]
    * Gruneisen: Gruneisen parameter, [-]
    * virialB: Second virial coefficient, [m³/kg]
    * virialC: Third virial coefficient, [m⁶/kg²]
    * dpdT_rho: Derivatives, dp/dT at constant rho, [MPa/K]
    * dpdrho_T: Derivatives, dp/drho at constant T, [MPa·m³/kg]
    * drhodT_P: Derivatives, drho/dT at constant P, [kg/m³·K]
    * drhodP_T: Derivatives, drho/dP at constant T, [kg/m³·MPa]
    * dhdT_rho: Derivatives, dh/dT at constant rho, [kJ/kg·K]
    * dhdP_T: Isothermal throttling coefficient, [kJ/kg·MPa]
    * dhdT_P: Derivatives, dh/dT at constant P, [kJ/kg·K]
    * dhdrho_T: Derivatives, dh/drho at constant T, [kJ·m³/kg²]
    * dhdrho_P: Derivatives, dh/drho at constant P, [kJ·m³/kg²]
    * dhdP_rho: Derivatives, dh/dP at constant rho, [kJ/kg·MPa]
    * kt: Isothermal Expansion Coefficient, [-]
    * ks: Adiabatic Compressibility, [1/MPa]
    * Ks: Adiabatic bulk modulus, [MPa]
    * Kt: Isothermal bulk modulus, [MPa]

    * v0: Ideal specific volume, [m³/kg]
    * rho0: Ideal gas density, [kg/m³]
    * u0: Ideal specific internal energy, [kJ/kg]
    * h0: Ideal specific enthalpy, [kJ/kg]
    * s0: Ideal specific entropy, [kJ/kg·K]
    * a0: Ideal specific Helmholtz free energy, [kJ/kg]
    * g0: Ideal specific Gibbs free energy, [kJ/kg]
    * cp0: Ideal specific isobaric heat capacity, [kJ/kg·K]
    * cv0: Ideal specific isochoric heat capacity, [kJ/kg·K]
    * w0: Ideal speed of sound, [m/s]
    * gamma0: Ideal isoentropic exponent, [-]

    * w: Speed of sound, [m/s]
    * mu: Dynamic viscosity, [Pa·s]
    * nu: Kinematic viscosity, [m²/s]
    * k: Thermal conductivity, [W/m·K]
    * alfa: Thermal diffusivity, [m²/s]
    * sigma: Surface tension, [N/m]
    * epsilon: Dielectric constant, [-]
    * n: Refractive index, [-]
    * Prandt: Prandtl number, [-]
    * Pr: Reduced Pressure, [-]
    * Tr: Reduced Temperature, [-]
    * Hvap: Vaporization heat, [kJ/kg]
    * Svap: Vaporization entropy, [kJ/kg·K]

    * Z_rho: :math:`(Z-1)/\rho`, [m³/kg]
    * IntP: Internal pressure, [MPa]
    * invT: Negative reciprocal temperature, [1/K]
    * hInput: Specific heat input, [kJ/kg]
  """
  */
  CP = None
  _Pv = None
  _rhoL = None
  _rhoG = None

  kwargs =
  {
    "T": 0.0,
    "P": 0.0,
    "rho": 0.0,
    "v": 0.0,
    "h": None,
    "s": None,
    "u": None,
    "x": None,
    "l": 0.5893,
    "rho0": None,
    "T0": None,
    "x0": 0.5
  }

  status = 0
  msg = "Undefined"

  void __init__(self, **kwargs)
  {
    //Constructor, define common constant and initinialice kwargs
    self.R = self._constants["R"]/self._constants.get("M", self.M)
    self.Zc = self.Pc/self.rhoc/self.R/self.Tc
    self.kwargs = MEoS.kwargs.copy()
    self.__call__(**kwargs)
  }

  void __call__(self, **kwargs)
  {
    //Make instance callable to can add input parameter one to one
    // Alternative rho input
    if("rhom" in kwargs)
    {
      kwargs["rho"] = kwargs["rhom"]*self.M
      del kwargs["rhom"]
    }
    else if(kwargs.get("v", 0))
    {
      kwargs["rho"] = 1./kwargs["v"]
      del kwargs["v"]
    }
    else if(kwargs.get("vm", 0))
    {
      kwargs["rho"] = self.M/kwargs["vm"]
      del kwargs["vm"]
    }
    self.kwargs.update(kwargs)

    if(self.calculable)
    {
      try:
      {
        self.status = 1
        self.calculo()
        self.msg = ""
      }
      except RuntimeError as err:
      {
        self.status = 0
        self.msg = err.args[0]
        raise(err)
      }

      // Add msg for extrapolation state
      if(self.name == "water" && 130 <= self.T < 273.15)
      {
        self.msg = "Extrapolated state"
        self.status = 3
        warnings.warn("Using extrapolated values")
      }
      else if(self.name == "water" && 50 <= self.T < 130)
      {
        self.msg = "Extrapolated state using Low-Temperature extension"
        self.status = 3
        warnings.warn("Using extrapolated values and Low-Temperature"
                "extension")
      }
    }
  }

  @property
  void calculable(self)
  {
    //Check if inputs are enough to define state
    self._mode = ""
    if(self.kwargs["T"] && self.kwargs["P"])
    {
      self._mode = "TP"
    }
    else if(self.kwargs["T"] && self.kwargs["rho"])
    {
      self._mode = "Trho"
    }
    else if(self.kwargs["T"] && self.kwargs["h"] is !None)
    {
      self._mode = "Th"
    }
    else if(self.kwargs["T"] && self.kwargs["s"] is !None)
    {
      self._mode = "Ts"
    }
    else if(self.kwargs["T"] && self.kwargs["u"] is !None)
    {
      self._mode = "Tu"
    }
    else if(self.kwargs["P"] && self.kwargs["rho"])
    {
      self._mode = "Prho"
    }
    else if(self.kwargs["P"] && self.kwargs["h"] is !None)
    {
      self._mode = "Ph"
    }
    else if(self.kwargs["P"] && self.kwargs["s"] is !None)
    {
      self._mode = "Ps"
    }
    else if(self.kwargs["P"] && self.kwargs["u"] is !None)
    {
      self._mode = "Pu"
    }
    else if(self.kwargs["rho"] && self.kwargs["h"] is !None)
    {
      self._mode = "rhoh"
    }
    else if(self.kwargs["rho"] && self.kwargs["s"] is !None)
    {
      self._mode = "rhos"
    }
    else if(self.kwargs["rho"] && self.kwargs["u"] is !None)
    {
      self._mode = "rhou"
    }
    else if(self.kwargs["h"] is !None && self.kwargs["s"] is !None)
    {
      self._mode = "hs"
    }
    else if(self.kwargs["h"] is !None && self.kwargs["u"] is !None)
    {
      self._mode = "hu"
    }
    else if(self.kwargs["s"] is !None && self.kwargs["u"] is !None)
    {
      self._mode = "su"
    }
    else if(self.kwargs["T"] && self.kwargs["x"] is !None)
    {
      self._mode = "Tx"
    }
    else if(self.kwargs["P"] && self.kwargs["x"] is !None)
    {
      self._mode = "Px"
    }
    return bool(self._mode)
  }

  void calculo(self)
  {
    //Calculate procedure
    T = self.kwargs["T"]
    rho = self.kwargs["rho"]
    P = self.kwargs["P"]
    s = self.kwargs["s"]
    h = self.kwargs["h"]
    u = self.kwargs["u"]
    x = self.kwargs["x"]

    // Initial values
    T0 = self.kwargs["T0"]
    rho0 = self.kwargs["rho0"]

    if(T0 || rho0)
    {
      To = T0
      rhoo = rho0
    }
    else if(self.name == "air")
    {
      To = 300
      rhoo = 1e-3
    }
    else
    {
      try:
      {
        st0 = IAPWS97(**self.kwargs)
      }
      except NotImplementedError:
      {
        To = 300
        rhoo = 900
      }
      else
      {
        if(st0.status)
        {
          To = st0.T
          rhoo = st0.rho
        }
        else
        {
          To = 300
          rhoo = 900
        }
      }
    }

    self.R = self._constants["R"]/self._constants.get("M", self.M)
    rhoc = self._constants.get("rhoref", self.rhoc)
    Tc = self._constants.get("Tref", self.Tc)

    propiedades = None

    if(self._mode !in ("Tx", "Px"))
    {
      // Method with iteration necessary to get x
      if(self._mode == "TP")
      {
        try:
        {
          if(self.name == "air")
          {
            raise ValueError
          }
          st0 = IAPWS97(**self.kwargs)
          rhoo = st0.rho
        }
        except NotImplementedError:
        {
          if(rho0)
          {
            rhoo = rho0
          }
          else if(T < self.Tc && P < self.Pc && self._Vapor_Pressure(T) < P)
          {
            rhoo = self._Liquid_Density(T)
          }
          else if(T < self.Tc && P < self.Pc)
          {
            rhoo = self._Vapor_Density(T)
          }
          else
          {
            rhoo = self.rhoc*3
          }
        }
        except ValueError:
        {
          rhoo = 1e-3
        }

        void f(rho)
        {
          delta = rho/rhoc
          tau = Tc/T

          fird = _phird(tau, delta, self._constants)
          Po = (1+delta*fird)*self.R*T*rho
          return Po-P*1000
        }

        rho = fsolve(f, rhoo)[0]

        // Calculate quality
        if(T > self.Tc)
        {
          x = 1
        }
        else
        {
          Ps = self._Vapor_Pressure(T)
          if(Ps*0.95 < P < Ps*1.05)
          {
            rhol, rhov, Ps = self._saturation(T)
            Ps *= 1e-3

          }
          if(Ps > P)
          {
            x = 1
          }
          else
          {
            x = 0
          }
        }
      }

      else if(self._mode == "Th")
      {
        tau = Tc/T
        ideal = self._phi0(tau, 1)
        fiot = ideal["fiot"]

        void f(rho)
        {
          delta = rho/rhoc
          fird = _phird(tau, delta, self._constants)
          firt = _phirt(tau, delta, self._constants)
          ho = self.R*T*(1+tau*(fiot+firt)+delta*fird)
          return ho-h
        }

        if(T >= self.Tc)
        {
          rhoo = self.rhoc
          rho = fsolve(f, rhoo)[0]
        }
        else
        {
          x0 = self.kwargs["x0"]
          rhov = self._Vapor_Density(T)
          rhol = self._Liquid_Density(T)
          deltaL = rhol/rhoc
          deltaG = rhov/rhoc

          firdL = _phird(tau, deltaL, self._constants)
          firtL = _phirt(tau, deltaL, self._constants)
          firdG = _phird(tau, deltaG, self._constants)
          firtG = _phirt(tau, deltaG, self._constants)
          hl = self.R*T*(1+tau*(fiot+firtL)+deltaL*firdL)
          hv = self.R*T*(1+tau*(fiot+firtG)+deltaG*firdG)
          if(x0 !in (0, 1) && hl <= h <= hv)
          {
            rhol, rhov, Ps = self._saturation(T)
            vapor = self._Helmholtz(rhov, T)
            liquido = self._Helmholtz(rhol, T)
            hv = vapor["h"]
            hl = liquido["h"]
            x = (h-hl)/(hv-hl)
            rho = 1/(x/rhov+(1-x)/rhol)
            P = Ps/1000
          }
          else
          {
            if(h > hv)
            {
              rhoo = rhov
            }
            else
            {
              rhoo = rhol
            }
            rho = fsolve(f, rhoo)[0]
          }
        }
      }

      else if(self._mode == "Ts")
      {
        tau = Tc/T

        void f(rho)
        {
          if(rho < 0)
          {
            rho = 1e-20
          }
          delta = rho/rhoc

          ideal = self._phi0(tau, delta)
          fio = ideal["fio"]
          fiot = ideal["fiot"]
          fir = _phir(tau, delta, self._constants)
          firt = _phirt(tau, delta, self._constants)
          so = self.R*(tau*(fiot+firt)-fio-fir)
          return so-s
        }

        if(T >= self.Tc)
        {
          rhoo = self.rhoc
          rho = fsolve(f, rhoo)[0]
        }
        else
        {
          rhov = self._Vapor_Density(T)
          rhol = self._Liquid_Density(T)
          deltaL = rhol/rhoc
          deltaG = rhov/rhoc

          idealL = self._phi0(tau, deltaL)
          idealG = self._phi0(tau, deltaG)
          fioL = idealL["fio"]
          fioG = idealG["fio"]
          fiot = idealL["fiot"]
          firL = _phir(tau, deltaL, self._constants)
          firtL = _phirt(tau, deltaL, self._constants)
          sl = self.R*(tau*(fiot+firtL)-fioL-firL)
          firG = _phir(tau, deltaG, self._constants)
          firtG = _phirt(tau, deltaG, self._constants)
          sv = self.R*(tau*(fiot+firtG)-fioG-firG)

          if(sl <= s <= sv)
          {
            rhol, rhov, Ps = self._saturation(T)
            vapor = self._Helmholtz(rhov, T)
            liquido = self._Helmholtz(rhol, T)
            sv = vapor["s"]
            sl = liquido["s"]
            x = (s-sl)/(sv-sl)
            rho = 1/(x/rhov+(1-x)/rhol)
            P = Ps/1000
          }
          else
          {
            if(s > sv)
            {
              rhoo = rhov
            }
            else
            {
              rhoo = rhol
            }
            rho = fsolve(f, rhoo)[0]
          }
        }
      }
      else if(self._mode == "Tu")
      {
        tau = Tc/T
        ideal = self._phi0(tau, 1)
        fiot = ideal["fiot"]

        void f(rho)
        {
          delta = rho/rhoc

          fird = _phird(tau, delta, self._constants)
          firt = _phirt(tau, delta, self._constants)
          Po = (1+delta*fird)*self.R*T*rho
          ho = self.R*T*(1+tau*(fiot+firt)+delta*fird)

          return ho-Po/rho-u
        }

        if(T >= self.Tc)
        {
          rhoo = self.rhoc
          rho = fsolve(f, rhoo)[0]
        }
        else
        {
          rhov = self._Vapor_Density(T)
          rhol = self._Liquid_Density(T)
          deltaL = rhol/rhoc
          deltaG = rhov/rhoc

          firdL = _phird(tau, deltaL, self._constants)
          firtL = _phirt(tau, deltaL, self._constants)
          firdG = _phird(tau, deltaG, self._constants)
          firtG = _phirt(tau, deltaG, self._constants)
          PoL = (1+deltaL*firdL)*self.R*T*rhol
          PoG = (1+deltaG*firdG)*self.R*T*rhov
          hoL = self.R*T*(1+tau*(fiot+firtL)+deltaL*firdL)
          hoG = self.R*T*(1+tau*(fiot+firtG)+deltaG*firdG)

          uv = hoG-PoG/rhov
          ul = hoL-PoL/rhol
          if(ul <= u <= uv)
          {
            rhol, rhov, Ps = self._saturation(T)
            vapor = self._Helmholtz(rhov, T)
            liquido = self._Helmholtz(rhol, T)
            uv = vapor["h"]-vapor["P"]/rhov
            ul = liquido["h"]-liquido["P"]/rhol
            x = (u-ul)/(uv-ul)
            rho = 1/(x/rhov-(1-x)/rhol)
            P = Ps/1000
          }
          else
          {
            if(u > uv)
            {
              rhoo = rhov
            }
            else
            {
              rhoo = rhol
            }
            rho = fsolve(f, rhoo)[0]
          }
        }
      }
      else if(self._mode == "Prho")
      {
        delta = rho/rhoc

        void f(T)
        {
          tau = Tc/T

          fird = _phird(tau, delta, self._constants)
          Po = (1+delta*fird)*self.R*T*rho
          return Po-P*1000
        }

        T = fsolve(f, To)[0]
        rhol = self._Liquid_Density(T)
        rhov = self._Vapor_Density(T)
        if(T == To || rhov <= rho <= rhol)
        {
          void f(parr)
          {
            T, rhol, rhog = parr
            tau = Tc/T
            deltaL = rhol/self.rhoc
            deltaG = rhog/self.rhoc

            firL = _phir(tau, deltaL, self._constants)
            firdL = _phird(tau, deltaL, self._constants)
            firG = _phir(tau, deltaG, self._constants)
            firdG = _phird(tau, deltaG, self._constants)

            Jl = rhol*(1+deltaL*firdL)
            Jv = rhog*(1+deltaG*firdG)
            K = firL-firG
            Ps = self.R*T*rhol*rhog/(rhol-rhog)*(K+log(rhol/rhog))
            return (Jl-Jv,
                Jl*(1/rhog-1/rhol)-log(rhol/rhog)-K,
                Ps - P*1000)
          }

          for to in [To, 300, 400, 500, 600]:
          {
            rhoLo = self._Liquid_Density(to)
            rhoGo = self._Vapor_Density(to)
            sol = fsolve(f, [to, rhoLo, rhoGo], full_output=True)
            T, rhoL, rhoG = sol[0]
            x = (1./rho-1/rhoL)/(1/rhoG-1/rhoL)
            sum = 0;
            for(int iter = 0; iter < sol[1]["fvec"].length; iter++
            {
              sum += abs(sol[1]["fvec"][iter]);
            }
            if(sol[2] == 1 && 0 <= x <= 1 && sum < 1e-5)
            {
              break
            }
          }

          sum = 0;
          for(int iter = 0; iter < sol[1]["fvec"].length; iter++
          {
            sum += abs(sol[1]["fvec"][iter]);
          }
          if(sum > 1e-5)
          {
            raise(RuntimeError(sol[3]))
          }

          liquido = self._Helmholtz(rhoL, T)
          vapor = self._Helmholtz(rhoG, T)
          P = self.R*T*rhoL*rhoG/(rhoL-rhoG)*(
            liquido["fir"]-vapor["fir"]+log(rhoL/rhoG))/1000
        }
      }
      else if(self._mode == "Ph")
      {
        void funcion(parr)
        {
          rho, T = parr
          delta = rho/rhoc
          tau = Tc/T

          ideal = self._phi0(tau, delta)
          fiot = ideal["fiot"]
          fird = _phird(tau, delta, self._constants)
          firt = _phirt(tau, delta, self._constants)
          Po = (1+delta*fird)*self.R*T*rho
          ho = self.R*T*(1+tau*(fiot+firt)+delta*fird)
          return Po-P*1000, ho-h
        }

        rho, T = fsolve(funcion, [rhoo, To])
        rhol = self._Liquid_Density(T)
        rhov = self._Vapor_Density(T)
        if(rho == rhoo || rhov <= rho <= rhol)
        {
          void f(parr)
          {
            T, rhol, rhog, x = parr
            tau = Tc/T
            deltaL = rhol/self.rhoc
            deltaG = rhog/self.rhoc

            ideal = self._phi0(tau, deltaL)
            fiot = ideal["fiot"]

            firL = _phir(tau, deltaL, self._constants)
            firdL = _phird(tau, deltaL, self._constants)
            firtL = _phirt(tau, deltaL, self._constants)
            hoL = self.R*T*(1+tau*(fiot+firtL)+deltaL*firdL)
            firG = _phir(tau, deltaG, self._constants)
            firdG = _phird(tau, deltaG, self._constants)
            firtG = _phirt(tau, deltaG, self._constants)
            hoG = self.R*T*(1+tau*(fiot+firtG)+deltaG*firdG)

            Jl = rhol*(1+deltaL*firdL)
            Jv = rhog*(1+deltaG*firdG)
            K = firL-firG
            Ps = self.R*T*rhol*rhog/(rhol-rhog)*(K+log(rhol/rhog))

            return (Jl-Jv,
                Jl*(1/rhog-1/rhol)-log(rhol/rhog)-K,
                hoL*(1-x)+hoG*x - h,
                Ps - P*1000)
          }

          for to in [To, 300, 400, 500, 600]:
          {
            rLo = self._Liquid_Density(to)
            rGo = self._Vapor_Density(to)
            sol = fsolve(f, [to, rLo, rGo, 0.5], full_output=True)
            T, rhoL, rhoG, x = sol[0]

            sum = 0;
            for(int iter = 0; iter < sol[1]["fvec"].length; iter++
            {
              sum += abs(sol[1]["fvec"][iter]);
            }
            if(sol[2] == 1 && 0 <= x <= 1 && sum < 1e-5)
            {
              break
            }
          }

          sum = 0;
          for(int iter = 0; iter < sol[1]["fvec"].length; iter++
          {
            sum += abs(sol[1]["fvec"][iter]);
          }
          if(sum > 1e-5)
          {
            raise(RuntimeError(sol[3]))
          }

          liquido = self._Helmholtz(rhoL, T)
          vapor = self._Helmholtz(rhoG, T)
          P = self.R*T*rhoL*rhoG/(rhoL-rhoG)*(
            liquido["fir"]-vapor["fir"]+log(rhoL/rhoG))/1000
        }
      }
      else if(self._mode == "Ps")
      {
        try:
        {
          x0 = st0.x
        }
        except NameError:
        {
          x0 = None

        }
        if(x0 is None || x0 == 0 || x0 == 1)
        {
          void f(parr)
          {
            rho, T = parr
            delta = rho/rhoc
            tau = Tc/T

            ideal = self._phi0(tau, delta)
            fio = ideal["fio"]
            fiot = ideal["fiot"]
            fird = _phird(tau, delta, self._constants)
            fir = _phir(tau, delta, self._constants)
            firt = _phirt(tau, delta, self._constants)
            Po = (1+delta*fird)*self.R*T*rho
            so = self.R*(tau*(fiot+firt)-fio-fir)
            return Po-P*1000, so-s
          }

          rho, T = fsolve(f, [rhoo, To])
        }
        else
        {
          void funcion(parr)
          {
            rho, T = parr
            rhol, rhov, Ps = self._saturation(T)
            vapor = self._Helmholtz(rhov, T)
            liquido = self._Helmholtz(rhol, T)
            x = (1./rho-1/rhol)/(1/rhov-1/rhol)
            return Ps-P*1000, vapor["s"]*x+liquido["s"]*(1-x)-s
          }
          rho, T = fsolve(funcion, [2., 500.])
          rhol, rhov, Ps = self._saturation(T)
          vapor = self._Helmholtz(rhov, T)
          liquido = self._Helmholtz(rhol, T)
          sv = vapor["s"]
          sl = liquido["s"]
          x = (s-sl)/(sv-sl)
        }
      }
      else if(self._mode == "Pu")
      {
        void f(parr)
        {
          rho, T = parr
          delta = rho/rhoc
          tau = Tc/T

          ideal = self._phi0(tau, delta)
          fiot = ideal["fiot"]
          fird = _phird(tau, delta, self._constants)
          firt = _phirt(tau, delta, self._constants)
          Po = (1+delta*fird)*self.R*T*rho
          ho = self.R*T*(1+tau*(fiot+firt)+delta*fird)
          return ho-Po/rho-u, Po-P*1000
        }

        sol = fsolve(f, [rhoo, To], full_output=True)
        rho, T = sol[0]
        rhol = self._Liquid_Density(T)
        rhov = self._Vapor_Density(T)
        if(rho == rhoo || sol[2] != 1)
        {
          void f(parr)
          {
            T, rhol, rhog, x = parr
            tau = Tc/T
            deltaL = rhol/self.rhoc
            deltaG = rhog/self.rhoc

            ideal = self._phi0(tau, deltaL)
            fiot = ideal["fiot"]

            firL = _phir(tau, deltaL, self._constants)
            firdL = _phird(tau, deltaL, self._constants)
            firtL = _phirt(tau, deltaL, self._constants)
            hoL = self.R*T*(1+tau*(fiot+firtL)+deltaL*firdL)
            firG = _phir(tau, deltaG, self._constants)
            firdG = _phird(tau, deltaG, self._constants)
            firtG = _phirt(tau, deltaG, self._constants)
            hoG = self.R*T*(1+tau*(fiot+firtG)+deltaG*firdG)

            Jl = rhol*(1+deltaL*firdL)
            Jv = rhog*(1+deltaG*firdG)
            K = firL-firG
            Ps = self.R*T*rhol*rhog/(rhol-rhog)*(K+log(rhol/rhog))
            vu = hoG-Ps/rhog
            lu = hoL-Ps/rhol
            return (Jl-Jv,
                Jl*(1/rhog-1/rhol)-log(rhol/rhog)-K,
                lu*(1-x)+vu*x - u,
                Ps - P*1000)

          }
          for to in [To, 300, 400, 500, 600]:
          {
            rLo = self._Liquid_Density(to)
            rGo = self._Vapor_Density(to)
            sol = fsolve(f, [to, rLo, rGo, 0.5], full_output=True)
            T, rhoL, rhoG, x = sol[0]
            sum = 0;
            for(int iter = 0; iter < sol[1]["fvec"].length; iter++
            {
              sum += abs(sol[1]["fvec"][iter]);
            }
            if(sol[2] == 1 && 0 <= x <= 1 && sum < 1e-5)
            {
              break
            }
          }

          sum = 0;
          for(int iter = 0; iter < sol[1]["fvec"].length; iter++
          {
            sum += abs(sol[1]["fvec"][iter]);
          }
          if(sum > 1e-5)
          {
            raise(RuntimeError(sol[3]))
          }

          liquido = self._Helmholtz(rhoL, T)
          vapor = self._Helmholtz(rhoG, T)
          P = self.R*T*rhoL*rhoG/(rhoL-rhoG)*(
          {
            liquido["fir"]-vapor["fir"]+log(rhoL/rhoG))/1000
          }
        }
      }
      else if(self._mode == "rhoh")
      {
        delta = rho/rhoc

        void f(T)
        {
          tau = Tc/T

          ideal = self._phi0(tau, delta)
          fiot = ideal["fiot"]
          fird = _phird(tau, delta, self._constants)
          firt = _phirt(tau, delta, self._constants)
          ho = self.R*T*(1+tau*(fiot+firt)+delta*fird)
          return ho-h
        }

        T = fsolve(f, To)[0]
        rhol = self._Liquid_Density(T)
        rhov = self._Vapor_Density(T)
        if(T == To || rhov <= rho <= rhol)
        {
          void f(parr)
          {
            T, rhol, rhog = parr
            tau = Tc/T
            deltaL = rhol/self.rhoc
            deltaG = rhog/self.rhoc

            ideal = self._phi0(tau, deltaL)
            fiot = ideal["fiot"]
            firL = _phir(tau, deltaL, self._constants)
            firdL = _phird(tau, deltaL, self._constants)
            firtL = _phirt(tau, deltaL, self._constants)
            hoL = self.R*T*(1+tau*(fiot+firtL)+deltaL*firdL)
            firG = _phir(tau, deltaG, self._constants)
            firdG = _phird(tau, deltaG, self._constants)
            firtG = _phirt(tau, deltaG, self._constants)
            hoG = self.R*T*(1+tau*(fiot+firtG)+deltaG*firdG)

            Jl = rhol*(1+deltaL*firdL)
            Jv = rhog*(1+deltaG*firdG)
            K = firL-firG
            x = (1./rho-1/rhol)/(1/rhog-1/rhol)
            return (Jl-Jv,
                Jl*(1/rhog-1/rhol)-log(rhol/rhog)-K,
                hoL*(1-x)+hoG*x - h)
          }
          for to in [To, 300, 400, 500, 600]:
          {
            rhoLo = self._Liquid_Density(to)
            rhoGo = self._Vapor_Density(to)
            sol = fsolve(f, [to, rhoLo, rhoGo], full_output=True)
            T, rhoL, rhoG = sol[0]
            x = (1./rho-1/rhoL)/(1/rhoG-1/rhoL)
            sum = 0;
            for(int iter = 0; iter < sol[1]["fvec"].length; iter++
            {
              sum += abs(sol[1]["fvec"][iter]);
            }
            if(sol[2] == 1 && 0 <= x <= 1 && sum < 1e-5)
            {
              break
            }
          }

          sum = 0;
          for(int iter = 0; iter < sol[1]["fvec"].length; iter++
          {
            sum += abs(sol[1]["fvec"][iter]);
          }
          if(sum > 1e-5)
          {
            raise(RuntimeError(sol[3]))
          }

          liquido = self._Helmholtz(rhoL, T)
          vapor = self._Helmholtz(rhoG, T)
          P = self.R*T*rhoL*rhoG/(rhoL-rhoG)*(
            liquido["fir"]-vapor["fir"]+log(rhoL/rhoG))/1000
        }
      }
      else if(self._mode == "rhos")
      {
        delta = rho/rhoc

        void f(T)
        {
          tau = Tc/T
          ideal = self._phi0(tau, delta)
          fio = ideal["fio"]
          fiot = ideal["fiot"]
          fir = _phir(tau, delta, self._constants)
          firt = _phirt(tau, delta, self._constants)
          so = self.R*(tau*(fiot+firt)-fio-fir)
          return so-s
        }

        T = fsolve(f, To)[0]
        rhol = self._Liquid_Density(T)
        rhov = self._Vapor_Density(T)
        if(T == To || rhov <= rho <= rhol)
        {
          void f(parr)
          {
            T, rhol, rhog = parr
            tau = Tc/T
            deltaL = rhol/self.rhoc
            deltaG = rhog/self.rhoc

            idealL = self._phi0(tau, deltaL)
            fioL = idealL["fio"]
            fiot = idealL["fiot"]
            idealG = self._phi0(tau, deltaG)
            fioG = idealG["fio"]

            firL = _phir(tau, deltaL, self._constants)
            firdL = _phird(tau, deltaL, self._constants)
            firtL = _phirt(tau, deltaL, self._constants)
            soL = self.R*(tau*(fiot+firtL)-fioL-firL)
            firG = _phir(tau, deltaG, self._constants)
            firdG = _phird(tau, deltaG, self._constants)
            firtG = _phirt(tau, deltaG, self._constants)
            soG = self.R*(tau*(fiot+firtG)-fioG-firG)

            Jl = rhol*(1+deltaL*firdL)
            Jv = rhog*(1+deltaG*firdG)
            K = firL-firG
            x = (1./rho-1/rhol)/(1/rhog-1/rhol)
            return (Jl-Jv,
                Jl*(1/rhog-1/rhol)-log(rhol/rhog)-K,
                soL*(1-x)+soG*x - s)
          }
          for to in [To, 300, 400, 500, 600]:
          {
            rhoLo = self._Liquid_Density(to)
            rhoGo = self._Vapor_Density(to)
            sol = fsolve(f, [to, rhoLo, rhoGo], full_output=True)
            T, rhoL, rhoG = sol[0]
            x = (1./rho-1/rhoL)/(1/rhoG-1/rhoL)
            sum = 0;
            for(int iter = 0; iter < sol[1]["fvec"].length; iter++
            {
              sum += abs(sol[1]["fvec"][iter]);
            }
            if(sol[2] == 1 && 0 <= x <= 1 && sum < 1e-5)
            {
              break
            }
          }

          sum = 0;
          for(int iter = 0; iter < sol[1]["fvec"].length; iter++
          {
            sum += abs(sol[1]["fvec"][iter]);
          }
          if(sum > 1e-5)
          {
            raise(RuntimeError(sol[3]))
          }

          liquido = self._Helmholtz(rhoL, T)
          vapor = self._Helmholtz(rhoG, T)
          P = self.R*T*rhoL*rhoG/(rhoL-rhoG)*(
            liquido["fir"]-vapor["fir"]+log(rhoL/rhoG))/1000
        }
      }
      else if(self._mode == "rhou")
      {
        delta = rho/rhoc

        void f(T)
        {
          tau = Tc/T

          ideal = self._phi0(tau, delta)
          fiot = ideal["fiot"]
          fird = _phird(tau, delta, self._constants)
          firt = _phirt(tau, delta, self._constants)
          Po = (1+delta*fird)*self.R*T*rho
          ho = self.R*T*(1+tau*(fiot+firt)+delta*fird)
          return ho-Po/rho-u
        }

        T = fsolve(f, To)[0]
        rhol = self._Liquid_Density(T)
        rhov = self._Vapor_Density(T)
        if(T == To || rhov <= rho <= rhol)
        {
          void f(parr)
          {
            T, rhol, rhog = parr
            tau = Tc/T
            deltaL = rhol/self.rhoc
            deltaG = rhog/self.rhoc

            ideal = self._phi0(tau, deltaL)
            fiot = ideal["fiot"]
            firL = _phir(tau, deltaL, self._constants)
            firdL = _phird(tau, deltaL, self._constants)
            firtL = _phirt(tau, deltaL, self._constants)
            hoL = self.R*T*(1+tau*(fiot+firtL)+deltaL*firdL)
            firG = _phir(tau, deltaG, self._constants)
            firdG = _phird(tau, deltaG, self._constants)
            firtG = _phirt(tau, deltaG, self._constants)
            hoG = self.R*T*(1+tau*(fiot+firtG)+deltaG*firdG)

            Jl = rhol*(1+deltaL*firdL)
            Jv = rhog*(1+deltaG*firdG)
            K = firL-firG
            x = (1./rho-1/rhol)/(1/rhog-1/rhol)
            Ps = self.R*T*rhol*rhog/(rhol-rhog)*(K+log(rhol/rhog))
            vu = hoG-Ps/rhog
            lu = hoL-Ps/rhol
            return (Jl-Jv,
                Jl*(1/rhog-1/rhol)-log(rhol/rhog)-K,
                lu*(1-x)+vu*x - u)
          }
          for to in [To, 300, 400, 500, 600]:
          {
            rhoLo = self._Liquid_Density(to)
            rhoGo = self._Vapor_Density(to)
            sol = fsolve(f, [to, rhoLo, rhoGo], full_output=True)
            T, rhoL, rhoG = sol[0]
            x = (1./rho-1/rhoL)/(1/rhoG-1/rhoL)
            sum = 0;
            for(int iter = 0; iter < sol[1]["fvec"].length; iter++
            {
              sum += abs(sol[1]["fvec"][iter]);
            }
            if(sol[2] == 1 && 0 <= x <= 1 && sum < 1e-5)
            {
              break
            }
          }

          sum = 0;
          for(int iter = 0; iter < sol[1]["fvec"].length; iter++
          {
            sum += abs(sol[1]["fvec"][iter]);
          }
          if(sum > 1e-5)
          {
            raise(RuntimeError(sol[3]))
          }
          liquido = self._Helmholtz(rhoL, T)
          vapor = self._Helmholtz(rhoG, T)
          P = self.R*T*rhoL*rhoG/(rhoL-rhoG)*(
            liquido["fir"]-vapor["fir"]+log(rhoL/rhoG))/1000
        }
      }
      else if(self._mode == "hs")
      {
        void f(parr)
        {
          rho, T = parr
          delta = rho/rhoc
          tau = Tc/T

          ideal = self._phi0(tau, delta)
          fio = ideal["fio"]
          fiot = ideal["fiot"]
          fird = _phird(tau, delta, self._constants)
          fir = _phir(tau, delta, self._constants)
          firt = _phirt(tau, delta, self._constants)
          ho = self.R*T*(1+tau*(fiot+firt)+delta*fird)
          so = self.R*(tau*(fiot+firt)-fio-fir)
          return ho-h, so-s
        }

        rho, T = fsolve(f, [rhoo, To])
        rhol = self._Liquid_Density(T)
        rhov = self._Vapor_Density(T)
        if(rhov <= rho <= rhol)
        {
          void f(parr)
          {
            T, rhol, rhog, x = parr
            tau = Tc/T
            deltaL = rhol/self.rhoc
            deltaG = rhog/self.rhoc

            idealL = self._phi0(tau, deltaL)
            fiot = idealL["fiot"]
            fioL = idealL["fio"]
            idealG = self._phi0(tau, deltaG)
            fioG = idealG["fio"]

            firL = _phir(tau, deltaL, self._constants)
            firdL = _phird(tau, deltaL, self._constants)
            firtL = _phirt(tau, deltaL, self._constants)
            hoL = self.R*T*(1+tau*(fiot+firtL)+deltaL*firdL)
            soL = self.R*(tau*(fiot+firtL)-fioL-firL)
            firG = _phir(tau, deltaG, self._constants)
            firdG = _phird(tau, deltaG, self._constants)
            firtG = _phirt(tau, deltaG, self._constants)
            hoG = self.R*T*(1+tau*(fiot+firtG)+deltaG*firdG)
            soG = self.R*(tau*(fiot+firtG)-fioG-firG)

            Jl = rhol*(1+deltaL*firdL)
            Jv = rhog*(1+deltaG*firdG)
            K = firL-firG
            return (Jl-Jv,
                Jl*(1/rhog-1/rhol)-log(rhol/rhog)-K,
                hoL*(1-x)+hoG*x - h,
                soL*(1-x)+soG*x - s)
          }

          for to in [To, 300, 400, 500, 600]:
          {
            rLo = self._Liquid_Density(to)
            rGo = self._Vapor_Density(to)
            sol = fsolve(f, [to, rLo, rGo, 0.5], full_output=True)
            T, rhoL, rhoG, x = sol[0]
            sum = 0;
            for(int iter = 0; iter < sol[1]["fvec"].length; iter++
            {
              sum += abs(sol[1]["fvec"][iter]);
            }
            if(sol[2] == 1 && 0 <= x <= 1 && sum < 1e-5)
            {
              break
            }
          }

          sum = 0;
          for(int iter = 0; iter < sol[1]["fvec"].length; iter++
          {
            sum += abs(sol[1]["fvec"][iter]);
          }
          if(sum > 1e-5)
          {
            raise(RuntimeError(sol[3]))
          }

          liquido = self._Helmholtz(rhoL, T)
          vapor = self._Helmholtz(rhoG, T)
          P = self.R*T*rhoL*rhoG/(rhoL-rhoG)*(
            liquido["fir"]-vapor["fir"]+log(rhoL/rhoG))/1000
        }
      }
      else if(self._mode == "hu")
      {
        void f(parr)
        {
          rho, T = parr
          delta = rho/rhoc
          tau = Tc/T

          ideal = self._phi0(tau, delta)
          fiot = ideal["fiot"]
          fird = _phird(tau, delta, self._constants)
          firt = _phirt(tau, delta, self._constants)
          Po = (1+delta*fird)*self.R*T*rho
          ho = self.R*T*(1+tau*(fiot+firt)+delta*fird)

          return ho-Po/rho-u, ho-h
        }

        sol = fsolve(f, [rhoo, To], full_output=True)
        rho, T = sol[0]
        rhol = self._Liquid_Density(T)
        rhov = self._Vapor_Density(T)
        if(sol[2] != 1 || rhov <= rho <= rhol)
        {
          void f(parr)
          {
            T, rhol, rhog, x = parr
            tau = Tc/T
            deltaL = rhol/self.rhoc
            deltaG = rhog/self.rhoc

            ideal = self._phi0(tau, deltaL)
            fiot = ideal["fiot"]

            firL = _phir(tau, deltaL, self._constants)
            firdL = _phird(tau, deltaL, self._constants)
            firtL = _phirt(tau, deltaL, self._constants)
            hoL = self.R*T*(1+tau*(fiot+firtL)+deltaL*firdL)
            firG = _phir(tau, deltaG, self._constants)
            firdG = _phird(tau, deltaG, self._constants)
            firtG = _phirt(tau, deltaG, self._constants)
            hoG = self.R*T*(1+tau*(fiot+firtG)+deltaG*firdG)

            Jl = rhol*(1+deltaL*firdL)
            Jv = rhog*(1+deltaG*firdG)
            K = firL-firG

            Ps = self.R*T*rhol*rhog/(rhol-rhog)*(K+log(rhol/rhog))
            vu = hoG-Ps/rhog
            lu = hoL-Ps/rhol

            return (Jl-Jv,
                Jl*(1/rhog-1/rhol)-log(rhol/rhog)-K,
                hoL*(1-x)+hoG*x - h,
                lu*(1-x)+vu*x - u)
          }

          for to in [To, 300, 400, 500, 600]:
          {
            rLo = self._Liquid_Density(to)
            rGo = self._Vapor_Density(to)
            sol = fsolve(f, [to, rLo, rGo, 0.5], full_output=True)
            T, rhoL, rhoG, x = sol[0]
            sum = 0;
            for(int iter = 0; iter < sol[1]["fvec"].length; iter++
            {
              sum += abs(sol[1]["fvec"][iter]);
            }
            if(sol[2] == 1 && 0 <= x <= 1 && sum < 1e-5)
              break
          }

          sum = 0;
          for(int iter = 0; iter < sol[1]["fvec"].length; iter++
          {
            sum += abs(sol[1]["fvec"][iter]);
          }
          if(sum > 1e-5)
          {
            raise(RuntimeError(sol[3]))
          }

          liquido = self._Helmholtz(rhoL, T)
          vapor = self._Helmholtz(rhoG, T)
          P = self.R*T*rhoL*rhoG/(rhoL-rhoG)*(
            liquido["fir"]-vapor["fir"]+log(rhoL/rhoG))/1000
        }
      }
      else if(self._mode == "su")
      {
        void f(parr)
        {
          rho, T = parr
          delta = rho/rhoc
          tau = Tc/T

          ideal = self._phi0(tau, delta)
          fio = ideal["fio"]
          fiot = ideal["fiot"]
          fird = _phird(tau, delta, self._constants)
          fir = _phir(tau, delta, self._constants)
          firt = _phirt(tau, delta, self._constants)
          ho = self.R*T*(1+tau*(fiot+firt)+delta*fird)
          so = self.R*(tau*(fiot+firt)-fio-fir)
          Po = (1+delta*fird)*self.R*T*rho
          return ho-Po/rho-u, so-s
        }

        sol = fsolve(f, [rhoo, To], full_output=True)
        rho, T = sol[0]
        rhol = self._Liquid_Density(T)
        rhov = self._Vapor_Density(T)
        if(sol[2] != 1 || rhov <= rho <= rhol)
        {
          void f(parr)
          {
            T, rhol, rhog, x = parr
            tau = Tc/T
            deltaL = rhol/self.rhoc
            deltaG = rhog/self.rhoc

            idealL = self._phi0(tau, deltaL)
            fiot = idealL["fiot"]
            fioL = idealL["fio"]
            idealG = self._phi0(tau, deltaG)
            fioG = idealG["fio"]

            firL = _phir(tau, deltaL, self._constants)
            firdL = _phird(tau, deltaL, self._constants)
            firtL = _phirt(tau, deltaL, self._constants)
            hoL = self.R*T*(1+tau*(fiot+firtL)+deltaL*firdL)
            soL = self.R*(tau*(fiot+firtL)-fioL-firL)
            firG = _phir(tau, deltaG, self._constants)
            firdG = _phird(tau, deltaG, self._constants)
            firtG = _phirt(tau, deltaG, self._constants)
            hoG = self.R*T*(1+tau*(fiot+firtG)+deltaG*firdG)
            soG = self.R*(tau*(fiot+firtG)-fioG-firG)

            Jl = rhol*(1+deltaL*firdL)
            Jv = rhog*(1+deltaG*firdG)
            K = firL-firG

            Ps = self.R*T*rhol*rhog/(rhol-rhog)*(K+log(rhol/rhog))
            vu = hoG-Ps/rhog
            lu = hoL-Ps/rhol

            return (Jl-Jv,
                Jl*(1/rhog-1/rhol)-log(rhol/rhog)-K,
                soL*(1-x)+soG*x - s,
                lu*(1-x)+vu*x - u)
          }

          for to in [To, 300, 400, 500, 600]:
          {
            rLo = self._Liquid_Density(to)
            rGo = self._Vapor_Density(to)
            sol = fsolve(f, [to, rLo, rGo, 0.5], full_output=True)
            T, rhoL, rhoG, x = sol[0]
            sum = 0;
            for(int iter = 0; iter < sol[1]["fvec"].length; iter++
            {
              sum += abs(sol[1]["fvec"][iter]);
            }
            if(sol[2] == 1 && 0 <= x <= 1 && sum < 1e-5)
              break
          }

          sum = 0;
          for(int iter = 0; iter < sol[1]["fvec"].length; iter++
          {
            sum += abs(sol[1]["fvec"][iter]);
          }
          if(sum > 1e-5)
          {
            raise(RuntimeError(sol[3]))
          }

          liquido = self._Helmholtz(rhoL, T)
          vapor = self._Helmholtz(rhoG, T)
          P = self.R*T*rhoL*rhoG/(rhoL-rhoG)*(
            liquido["fir"]-vapor["fir"]+log(rhoL/rhoG))/1000
        }
      }
      else if(self._mode == "Trho")
      {
        if(T < self.Tc)
        {
          rhov = self._Vapor_Density(T)
          rhol = self._Liquid_Density(T)
          if(rhol > rho > rhov)
          {
            rhol, rhov, Ps = self._saturation(T)
            if(rhol > rho > rhov)
            {
              vapor = self._Helmholtz(rhov, T)
              liquido = self._Helmholtz(rhol, T)
              x = (1/rho-1/rhol)/(1/rhov-1/rhol)
              P = Ps/1000
            }
          }
        }
      }
      rho = float(rho)
      T = float(T)
      propiedades = self._Helmholtz(rho, T)

      if(T > self.Tc)
      {
        x = 1
      }
      else if(x is None)
      {
        x = 0
      }

      if(!P)
      {
        P = propiedades["P"]/1000.
      }
    }
    else if(self._mode == "Tx")
    {
      // Check input T in saturation range
      if(self.Tt > T || self.Tc < T || x > 1 || x < 0)
      {
        raise NotImplementedError("Incoming out of bound")
      }

      rhol, rhov, Ps = self._saturation(T)
      vapor = self._Helmholtz(rhov, T)
      liquido = self._Helmholtz(rhol, T)
      if(x == 0)
      {
        propiedades = liquido
      }
      else if(x == 1)
      {
        propiedades = vapor
      }
      P = Ps/1000.
    }
    else if(self._mode == "Px")
    {
      // Check input P in saturation range
      if(self.Pc < P || x > 1 || x < 0)
      {
        raise NotImplementedError("Incoming out of bound")
      }

      // Iterate over saturation routine to get T
      void f(T)
      {
        rhol = self._Liquid_Density(T)
        rhog = self._Vapor_Density(T)

        deltaL = rhol/self.rhoc
        deltaG = rhog/self.rhoc

        tau = Tc/T

        firL = _phir(tau, deltaL, self._constants)
        firG = _phir(tau, deltaG, self._constants)
        Ps = self.R*T*rhol*rhog/(rhol-rhog)*(
          firL-firG+log(deltaL/deltaG))
        return Ps/1000-P
      }
      if(T0)
      {
        To = T0
      }
      else if(self.name == "water")
      {
        To = _TSat_P(P)
      }
      else
      {
        To = (self.Tc+self.Tt)/2
      }
      T = fsolve(f, To)[0]
      rhol, rhov, Ps = self._saturation(T)
      vapor = self._Helmholtz(rhov, T)
      liquido = self._Helmholtz(rhol, T)
      if(x == 0)
      {
        propiedades = liquido
      }
      else if(x == 1)
      {
        propiedades = vapor
      }
    }
    self.T = T
    self.Tr = T/self.Tc
    self.P = P
    self.Pr = self.P/self.Pc
    self.x = x
    if(self._mode in ["Tx", "Px"] || 0 < x < 1)
    {
      region = 4
    }
    else
    {
      region = 0
    }
    self.phase = getphase(self.Tc, self.Pc, self.T, self.P, self.x, region)

    self.Liquid = _fase()
    self.Gas = _fase()
    if(x == 0)
    {
      // liquid phase
      self.fill(self.Liquid, propiedades)
      self.fill(self, propiedades)
    }
    else if(x == 1)
    {
      // vapor phase
      self.fill(self.Gas, propiedades)
      self.fill(self, propiedades)
    }
    else
    {
      self.fill(self.Liquid, liquido)
      self.fill(self.Gas, vapor)

      self.v = x*self.Gas.v+(1-x)*self.Liquid.v
      self.rho = 1./self.v

      self.h = x*self.Gas.h+(1-x)*self.Liquid.h
      self.s = x*self.Gas.s+(1-x)*self.Liquid.s
      self.u = x*self.Gas.u+(1-x)*self.Liquid.u
      self.a = x*self.Gas.a+(1-x)*self.Liquid.a
      self.g = x*self.Gas.g+(1-x)*self.Liquid.g

      self.Z = x*self.Gas.Z+(1-x)*self.Liquid.Z
      self.f = x*self.Gas.f+(1-x)*self.Liquid.f

      self.Z_rho = x*self.Gas.Z_rho+(1-x)*self.Liquid.Z_rho
      self.IntP = x*self.Gas.IntP+(1-x)*self.Liquid.IntP
    }

    // Calculate special properties useful only for one phase
    if(self._mode in ("Px", "Tx") || (x < 1 && self.Tt <= T <= self.Tc))
    {
      self.sigma = self._surface(T)
    }
    else
    {
      self.sigma = None
    }

    vir = self._virial(T)
    self.virialB = vir["B"]/self.rhoc
    self.virialC = vir["C"]/self.rhoc**2

    if(0 < x < 1)
    {
      self.Hvap = vapor["h"]-liquido["h"]
      self.Svap = vapor["s"]-liquido["s"]
    }
    else
    {
      self.Hvap = None
      self.Svap = None
    }
    self.invT = -1/self.T

    // Ideal properties
    cp0 = self._prop0(self.rho, self.T)
    self.v0 = self.R*self.T/self.P/1000
    self.rho0 = 1./self.v0
    self.h0 = cp0.h
    self.u0 = self.h0-self.P*self.v0
    self.s0 = cp0.s
    self.a0 = self.u0-self.T*self.s0
    self.g0 = self.h0-self.T*self.s0
    self.cp0 = cp0.cp
    self.cv0 = cp0.cv
    self.cp0_cv = self.cp0/self.cv0
    cp0.v = self.v0
    self.gamma0 = -self.v0/self.P/1000*self.derivative("P", "v", "s", cp0)
  }

  void fill(self, fase, estado)
  {
    //Fill phase properties
    fase.rho = estado["rho"]
    fase.v = 1/fase.rho

    fase.h = estado["h"]
    fase.s = estado["s"]
    fase.u = fase.h-self.P*1000*fase.v
    fase.a = fase.u-self.T*fase.s
    fase.g = fase.h-self.T*fase.s

    fase.Z = self.P*fase.v/self.T/self.R*1e3
    fase.fi = exp(estado["fir"]+estado["delta"]*estado["fird"] -
            log(1+estado["delta"]*estado["fird"]))
    fase.f = fase.fi*self.P
    fase.cv = estado["cv"]

    fase.rhoM = fase.rho/self.M
    fase.hM = fase.h*self.M
    fase.sM = fase.s*self.M
    fase.uM = fase.u*self.M
    fase.aM = fase.a*self.M
    fase.gM = fase.g*self.M

    fase.alfap = estado["alfap"]
    fase.betap = estado["betap"]

    fase.cp = self.derivative("h", "T", "P", fase)
    fase.cp_cv = fase.cp/fase.cv
    fase.w = (self.derivative("P", "rho", "s", fase)*1000)**0.5
    fase.cvM = fase.cv*self.M
    fase.cpM = fase.cp*self.M

    fase.joule = self.derivative("T", "P", "h", fase)*1e3
    fase.Gruneisen = fase.v/fase.cv*self.derivative("P", "T", "v", fase)
    fase.alfav = self.derivative("v", "T", "P", fase)/fase.v
    fase.kappa = -self.derivative("v", "P", "T", fase)/fase.v*1e3
    fase.betas = self.derivative("T", "P", "s", fase)

    fase.gamma = -fase.v/self.P*self.derivative("P", "v", "s", fase)*1e-3
    fase.kt = -fase.v/self.P*self.derivative("P", "v", "T", fase)*1e-3
    fase.ks = -self.derivative("v", "P", "s", fase)/fase.v*1e3
    fase.Kt = -fase.v*self.derivative("P", "v", "s", fase)*1e-3
    fase.Ks = -fase.v*self.derivative("P", "v", "T", fase)*1e-3
    fase.dhdT_rho = self.derivative("h", "T", "rho", fase)
    fase.dhdT_P = self.derivative("h", "T", "P", fase)
    fase.dhdP_T = self.derivative("h", "P", "T", fase)*1e3
    fase.dhdP_rho = self.derivative("h", "P", "rho", fase)*1e3
    fase.dhdrho_T = self.derivative("h", "rho", "T", fase)
    fase.dhdrho_P = self.derivative("h", "rho", "P", fase)
    fase.dpdT_rho = self.derivative("P", "T", "rho", fase)*1e-3
    fase.dpdrho_T = self.derivative("P", "rho", "T", fase)*1e-3
    fase.drhodP_T = self.derivative("rho", "P", "T", fase)*1e3
    fase.drhodT_P = self.derivative("rho", "T", "P", fase)

    fase.Z_rho = (fase.Z-1)/fase.rho
    fase.IntP = self.T*self.derivative("P", "T", "rho", fase)*1e-3-self.P
    fase.hInput = fase.v*self.derivative("h", "v", "P", fase)

    fase.mu = self._visco(fase.rho, self.T, fase)
    fase.k = self._thermo(fase.rho, self.T, fase)
    fase.nu = fase.mu/fase.rho
    fase.alfa = fase.k/1000/fase.rho/fase.cp
    fase.Prandt = fase.mu*fase.cp*1000/fase.k
    if(self.name == "water")
    {
      try:
      {
        fase.epsilon = _Dielectric(fase.rho, self.T)
      }
      except NotImplementedError:
      {
        fase.epsilon = None
      }

      try:
      {
        fase.n = _Refractive(fase.rho, self.T, self.kwargs["l"])
      }
      except NotImplementedError:
      {
        fase.n = None
      }
    }
    else
    {
      fase.epsilon = None
      fase.n = None
    }
  }

  void derivative(self, z, x, y, fase)
  {
    /*
    """Wrapper derivative for custom derived properties
    where x, y, z can be: P, T, v, rho, u, h, s, g, a"""
    */
    return deriv_H(self, z, x, y, fase)
  }

  void _saturation(self, T)
  {
    //Saturation calculation for two phase search
    rhoc = self._constants.get("rhoref", self.rhoc)
    Tc = self._constants.get("Tref", self.Tc)

    if(T > Tc)
    {
      T = Tc
    }
    tau = Tc/T

    rhoLo = self._Liquid_Density(T)
    rhoGo = self._Vapor_Density(T)

    void f(parr)
    {
      rhol, rhog = parr
      deltaL = rhol/rhoc
      deltaG = rhog/rhoc
      phirL = _phir(tau, deltaL, self._constants)
      phirG = _phir(tau, deltaG, self._constants)
      phirdL = _phird(tau, deltaL, self._constants)
      phirdG = _phird(tau, deltaG, self._constants)
      Jl = deltaL*(1+deltaL*phirdL)
      Jv = deltaG*(1+deltaG*phirdG)
      Kl = deltaL*phirdL+phirL+log(deltaL)
      Kv = deltaG*phirdG+phirG+log(deltaG)
      return Kv-Kl, Jv-Jl
    }

    rhoL, rhoG = fsolve(f, [rhoLo, rhoGo])
    if(rhoL == rhoG)
    {
      Ps = self.Pc
    }
    else
    {
      deltaL = rhoL/self.rhoc
      deltaG = rhoG/self.rhoc
      firL = _phir(tau, deltaL, self._constants)
      firG = _phir(tau, deltaG, self._constants)

      Ps = self.R*T*rhoL*rhoG/(rhoL-rhoG)*(firL-firG+log(deltaL/deltaG))
    }
    return rhoL, rhoG, Ps
  }

  void _Helmholtz(self, rho, T)
  {
    /*
    """Calculated properties from helmholtz free energy and derivatives

    Parameters
    ----------
    rho : float
      Density, [kg/m³]
    T : float
      Temperature, [K]

    Returns
    -------
    prop : dict
      Dictionary with calculated properties:
        * fir: [-]
        * fird: ∂fir/∂δ|τ
        * firdd: ∂²fir/∂δ²|τ
        * delta: Reducen density rho/rhoc, [-]
        * P: Pressure, [kPa]
        * h: Enthalpy, [kJ/kg]
        * s: Entropy, [kJ/kgK]
        * cv: Isochoric specific heat, [kJ/kgK]
        * alfav: Thermal expansion coefficient, [1/K]
        * betap: Isothermal compressibility, [1/kPa]

    References
    ----------
    IAPWS, Revised Release on the IAPWS Formulation 1995 for the
    Thermodynamic Properties of Ordinary Water Substance for General and
    Scientific Use, September 2016, Table 3
    http://www.iapws.org/relguide/IAPWS-95.html
    """
    */
    if(isinstance(rho, ndarray))
    {
      rho = rho[0]
    }
    if(isinstance(T, ndarray))
    {
      T = T[0]
    }
    if(rho < 0)
    {
      rho = 1e-20
    }
    if(T < 50)
    {
      T = 50
    }
    rhoc = self._constants.get("rhoref", self.rhoc)
    Tc = self._constants.get("Tref", self.Tc)
    delta = rho/rhoc
    tau = Tc/T
    ideal = self._phi0(tau, delta)
    fio = ideal["fio"]
    fiot = ideal["fiot"]
    fiott = ideal["fiott"]

    res = self._phir(tau, delta)
    fir = res["fir"]
    firt = res["firt"]
    firtt = res["firtt"]
    fird = res["fird"]
    firdd = res["firdd"]
    firdt = res["firdt"]

    propiedades = {}
    propiedades["fir"] = fir
    propiedades["fird"] = fird
    propiedades["firdd"] = firdd
    propiedades["delta"] = delta

    propiedades["rho"] = rho
    propiedades["P"] = (1+delta*fird)*self.R*T*rho
    propiedades["h"] = self.R*T*(1+tau*(fiot+firt)+delta*fird)
    propiedades["s"] = self.R*(tau*(fiot+firt)-fio-fir)
    propiedades["cv"] = -self.R*tau**2*(fiott+firtt)
    propiedades["alfap"] = (1-delta*tau*firdt/(1+delta*fird))/T
    propiedades["betap"] = rho*(
      1+(delta*fird+delta**2*firdd)/(1+delta*fird))
    return propiedades
  }

  void _prop0(self, rho, T)
  {
    //Ideal gas properties
    rhoc = self._constants.get("rhoref", self.rhoc)
    Tc = self._constants.get("Tref", self.Tc)
    delta = rho/rhoc
    tau = Tc/T
    ideal = self._phi0(tau, delta)
    fio = ideal["fio"]
    fiot = ideal["fiot"]
    fiott = ideal["fiott"]

    propiedades = _fase()
    propiedades.h = self.R*T*(1+tau*fiot)
    propiedades.s = self.R*(tau*fiot-fio)
    propiedades.cv = -self.R*tau**2*fiott
    propiedades.cp = self.R*(-tau**2*fiott+1)
    propiedades.alfap = 1/T
    propiedades.betap = rho
    return propiedades
  }

  void _phi0(self, tau, delta)
  {
    /*
    """Ideal gas Helmholtz free energy and derivatives

    Parameters
    ----------
    tau : float
      Inverse reduced temperature Tc/T, [-]
    delta : float
      Reduced density rho/rhoc, [-]

    Returns
    -------
    prop : dictionary with ideal adimensional helmholtz energy and deriv
      fio, [-]
      fiot: ∂fio/∂τ|δ
      fiod: ∂fio/∂δ|τ
      fiott: ∂²fio/∂τ²|δ
      fiodt: ∂²fio/∂τ∂δ
      fiodd: ∂²fio/∂δ²|τ

    References
    ----------
    IAPWS, Revised Release on the IAPWS Formulation 1995 for the
    Thermodynamic Properties of Ordinary Water Substance for General and
    Scientific Use, September 2016, Table 4
    http://www.iapws.org/relguide/IAPWS-95.html
    """
    */
    Fi0 = self.Fi0

    fio = Fi0["ao_log"][0]*log(delta)+Fi0["ao_log"][1]*log(tau)
    fiot = +Fi0["ao_log"][1]/tau
    fiott = -Fi0["ao_log"][1]/tau**2

    fiod = 1/delta
    fiodd = -1/delta**2
    fiodt = 0

    for(int i = 0; i < Fi0.length; i++)
    {
      n = Fi0["ao_pow"][i];
      t = Fi0["pow"][i];
      fio += n*tau**t
      if(t != 0)
      {
        fiot += t*n*tau**(t-1)
      }
      if(t !in [0, 1])
      {
        fiott += n*t*(t-1)*tau**(t-2)
      }
    }

    for(int i = 0; i < Fi0.length; i++)
    {
      n = Fi0["ao_exp"][i];
      t = Fi0["titao"][i];
      fio += n*log(1-exp(-tau*t))
      fiot += n*t*((1-exp(-t*tau))**-1-1)
      fiott -= n*t**2*exp(-t*tau)*(1-exp(-t*tau))**-2
    }

    // Extension to especial terms of air
    if("ao_exp2" in Fi0)
    {
      for(int i = 0; i < Fi0["ao_exp2"].length; i++)
      {
        n = Fi0["ao_exp2"][i];
        g = Fi0["titao2"][i];
        C = Fi0["sum2"][i];
        fio += n*log(C+exp(g*tau))
        fiot += n*g/(C*exp(-g*tau)+1)
        fiott += C*n*g**2*exp(-g*tau)/(C*exp(-g*tau)+1)**2
      }
    }

    prop = {}
    prop["fio"] = fio
    prop["fiot"] = fiot
    prop["fiott"] = fiott
    prop["fiod"] = fiod
    prop["fiodd"] = fiodd
    prop["fiodt"] = fiodt
    return prop
  }

  void _phir(self, tau, delta)
  {
    /*
    """Residual contribution to the free Helmholtz energy

    Parameters
    ----------
    tau : float
      Inverse reduced temperature Tc/T, [-]
    delta : float
      Reduced density rho/rhoc, [-]

    Returns
    -------
    prop : dict
      Dictionary with residual adimensional helmholtz energy and deriv:
      * fir
      * firt: ∂fir/∂τ|δ,x
      * fird: ∂fir/∂δ|τ,x
      * firtt: ∂²fir/∂τ²|δ,x
      * firdt: ∂²fir/∂τ∂δ|x
      * firdd: ∂²fir/∂δ²|τ,x

    References
    ----------
    IAPWS, Revised Release on the IAPWS Formulation 1995 for the
    Thermodynamic Properties of Ordinary Water Substance for General and
    Scientific Use, September 2016, Table 5
    http://www.iapws.org/relguide/IAPWS-95.html
    """
    */
    fir = fird = firdd = firt = firtt = firdt = 0

    // Polinomial terms
    nr1 = self._constants.get("nr1", [])
    d1 = self._constants.get("d1", [])
    t1 = self._constants.get("t1", [])
    for(int i = 0; i < nr1.length; i++)
    {
      n = nr1[i];
      d = d1[i];
      t = t1[i];
      fir += n*delta**d*tau**t
      fird += n*d*delta**(d-1)*tau**t
      firdd += n*d*(d-1)*delta**(d-2)*tau**t
      firt += n*t*delta**d*tau**(t-1)
      firtt += n*t*(t-1)*delta**d*tau**(t-2)
      firdt += n*t*d*delta**(d-1)*tau**(t-1)
    }

    // Exponential terms
    nr2 = self._constants.get("nr2", [])
    d2 = self._constants.get("d2", [])
    g2 = self._constants.get("gamma2", [])
    t2 = self._constants.get("t2", [])
    c2 = self._constants.get("c2", [])
    for(int iter = 0; iter < nr2.length; iter++)
    {
      n = nr2[iter];
      d = d2[iter];
      g = g2[iter];
      t = t2[iter];
      c = c2[iter];
      fir += n*delta**d*tau**t*exp(-g*delta**c)
      fird += n*exp(-g*delta**c)*delta**(d-1)*tau**t*(d-g*c*delta**c)
      firdd += n*exp(-g*delta**c)*delta**(d-2)*tau**t * \
        ((d-g*c*delta**c)*(d-1-g*c*delta**c)-g**2*c**2*delta**c)
      firt += n*t*delta**d*tau**(t-1)*exp(-g*delta**c)
      firtt += n*t*(t-1)*delta**d*tau**(t-2)*exp(-g*delta**c)
      firdt += n*t*delta**(d-1)*tau**(t-1)*(d-g*c*delta**c)*exp(
        -g*delta**c)
    }

    // Gaussian terms
    nr3 = self._constants.get("nr3", [])
    d3 = self._constants.get("d3", [])
    t3 = self._constants.get("t3", [])
    a3 = self._constants.get("alfa3", [])
    e3 = self._constants.get("epsilon3", [])
    b3 = self._constants.get("beta3", [])
    g3 = self._constants.get("gamma3", [])
    for(int iter = 0; iter < nr3.length; iter++)
    {
      n = nr3[iter];
      d = d3[iter];
      t = t3[iter];
      a = a3[iter];
      e = e3[iter];
      b = b3[iter];
      g = g3[iter];
      fir += n*delta**d*tau**t*exp(-a*(delta-e)**2-b*(tau-g)**2)
      fird += n*delta**d*tau**t*exp(-a*(delta-e)**2-b*(tau-g)**2)*(
          d/delta-2*a*(delta-e))
      firdd += n*tau**t*exp(-a*(delta-e)**2-b*(tau-g)**2)*(
          -2*a*delta**d + 4*a**2*delta**d*(delta-e)**2 -
          4*d*a*delta**(d-1)*(delta-e) + d*(d-1)*delta**(d-2))
      firt += n*delta**d*tau**t*exp(-a*(delta-e)**2-b*(tau-g)**2)*(
          t/tau-2*b*(tau-g))
      firtt += n*delta**d*tau**t*exp(-a*(delta-e)**2-b*(tau-g)**2)*(
          (t/tau-2*b*(tau-g))**2-t/tau**2-2*b)
      firdt += n*delta**d*tau**t*exp(-a*(delta-e)**2-b*(tau-g)**2)*(
          t/tau-2*b*(tau-g))*(d/delta-2*a*(delta-e))
    }
    // Non analitic terms
    nr4 = self._constants.get("nr4", [])
    a4 = self._constants.get("a4", [])
    b4 = self._constants.get("b4", [])
    Ai = self._constants.get("A", [])
    Bi = self._constants.get("B", [])
    Ci = self._constants.get("C", [])
    Di = self._constants.get("D", [])
    bt4 = self._constants.get("beta4", [])
    for(int iter = 0; iter < nr4.length; iter++)
    {
      n = nr4[iter];
      a = a4[iter];
      b = b4[iter];
      A = Ai[iter];
      B = Bi[iter];
      C = Ci[iter];
      D = Di[iter];
      bt = bt4[iter];
      Tita = (1-tau)+A*((delta-1)**2)**(0.5/bt)
      F = exp(-C*(delta-1)**2-D*(tau-1)**2)
      Fd = -2*C*F*(delta-1)
      Fdd = 2*C*F*(2*C*(delta-1)**2-1)
      Ft = -2*D*F*(tau-1)
      Ftt = 2*D*F*(2*D*(tau-1)**2-1)
      Fdt = 4*C*D*F*(delta-1)*(tau-1)

      Delta = Tita**2+B*((delta-1)**2)**a
      Deltad = (delta-1)*(A*Tita*2/bt*((delta-1)**2)**(0.5/bt-1) +
                2*B*a*((delta-1)**2)**(a-1))
      if(delta == 1)
      {
        Deltadd = 0
      }
      else
      {
        Deltadd = Deltad/(delta-1)+(delta-1)**2*(
          4*B*a*(a-1)*((delta-1)**2)**(a-2) +
          2*A**2/bt**2*(((delta-1)**2)**(0.5/bt-1))**2 +
          A*Tita*4/bt*(0.5/bt-1)*((delta-1)**2)**(0.5/bt-2))
      }

      DeltaBd = b*Delta**(b-1)*Deltad
      DeltaBdd = b*(Delta**(b-1)*Deltadd+(b-1)*Delta**(b-2)*Deltad**2)
      DeltaBt = -2*Tita*b*Delta**(b-1)
      DeltaBtt = 2*b*Delta**(b-1)+4*Tita**2*b*(b-1)*Delta**(b-2)
      DeltaBdt = -A*b*2/bt*Delta**(b-1)*(delta-1)*((delta-1)**2)**(
        0.5/bt-1)-2*Tita*b*(b-1)*Delta**(b-2)*Deltad

      fir += n*Delta**b*delta*F
      fird += n*(Delta**b*(F+delta*Fd)+DeltaBd*delta*F)
      firdd += n*(Delta**b*(2*Fd+delta*Fdd) + 2*DeltaBd*(F+delta*Fd) +
            DeltaBdd*delta*F)
      firt += n*delta*(DeltaBt*F+Delta**b*Ft)
      firtt += n*delta*(DeltaBtt*F+2*DeltaBt*Ft+Delta**b*Ftt)
      firdt += n*(Delta**b*(Ft+delta*Fdt)+delta*DeltaBd*Ft +
            DeltaBt*(F+delta*Fd)+DeltaBdt*delta*F)
    }

    prop = {}
    prop["fir"] = fir
    prop["firt"] = firt
    prop["firtt"] = firtt
    prop["fird"] = fird
    prop["firdd"] = firdd
    prop["firdt"] = firdt
    return prop
  }

  void _virial(self, T)
  {
    /*
    """Virial coefficient

    Parameters
    ----------
    T : float
      Temperature [K]

    Returns
    -------
    prop : dict
      Dictionary with residual adimensional helmholtz energy:
        * B: ∂fir/∂δ|δ->0
        * C: ∂²fir/∂δ²|δ->0
    """
    */
    Tc = self._constants.get("Tref", self.Tc)
    tau = Tc/T
    B = C = 0
    delta = 1e-200

    // Polinomial terms
    nr1 = self._constants.get("nr1", [])
    d1 = self._constants.get("d1", [])
    t1 = self._constants.get("t1", [])
    for(int i = 0; i < nr1.length; i++)
    {
      n = nr1[i];
      d = d1[i];
      t = t1[i];
      B += n*d*delta**(d-1)*tau**t
      C += n*d*(d-1)*delta**(d-2)*tau**t
    }

    // Exponential terms
    nr2 = self._constants.get("nr2", [])
    d2 = self._constants.get("d2", [])
    g2 = self._constants.get("gamma2", [])
    t2 = self._constants.get("t2", [])
    c2 = self._constants.get("c2", [])
    for(int iter = 0; iter < nr2.length; iter++)
    {
      n = nr2[iter];
      d = d2[iter];
      g = g2[iter];
      t = t2[iter];
      c = c2[iter];
      B += n*exp(-g*delta**c)*delta**(d-1)*tau**t*(d-g*c*delta**c)
      C += n*exp(-g*delta**c)*(delta**(d-2)*tau**t*(
        (d-g*c*delta**c)*(d-1-g*c*delta**c)-g**2*c**2*delta**c))
    }

    // Gaussian terms
    nr3 = self._constants.get("nr3", [])
    d3 = self._constants.get("d3", [])
    t3 = self._constants.get("t3", [])
    a3 = self._constants.get("alfa3", [])
    e3 = self._constants.get("epsilon3", [])
    b3 = self._constants.get("beta3", [])
    g3 = self._constants.get("gamma3", [])
    for(int iter = 0; iter < nr3.length; iter++)
    {
      n = nr3[iter];
      d = d3[iter];
      t = t3[iter];
      a = a3[iter];
      e = e3[iter];
      b = b3[iter];
      g = g3[iter];
      B += n*delta**d*tau**t*exp(-a*(delta-e)**2-b*(tau-g)**2)*(
          d/delta-2*a*(delta-e))
      C += n*tau**t*exp(-a*(delta-e)**2-b*(tau-g)**2)*(
          -2*a*delta**d+4*a**2*delta**d*(
            delta-e)**2-4*d*a*delta**2*(
              delta-e)+d*2*delta)
    }

    // Non analitic terms
    nr4 = self._constants.get("nr4", [])
    a4 = self._constants.get("a4", [])
    b4 = self._constants.get("b4", [])
    Ai = self._constants.get("A", [])
    Bi = self._constants.get("B", [])
    Ci = self._constants.get("C", [])
    Di = self._constants.get("D", [])
    bt4 = self._constants.get("beta4", [])
    for(int iter = 0; iter < nr4.length; iter++)
    {
      n = nr4[iter];
      a = a4[iter];
      b = b4[iter];
      A = Ai[iter];
      B_ = Bi[iter];
      C_ = Ci[iter];
      D = Di[iter];
      bt = bt4[iter];
      Tita = (1-tau)+A*((delta-1)**2)**(0.5/bt)
      Delta = Tita**2+B_*((delta-1)**2)**a
      Deltad = (delta-1)*(A*Tita*2/bt*((delta-1)**2)**(
        0.5/bt-1)+2*B_*a*((delta-1)**2)**(a-1))
      Deltadd = Deltad/(delta-1) + (delta-1)**2*(
        4*B_*a*(a-1)*((delta-1)**2)**(a-2) +
        2*A**2/bt**2*(((delta-1)**2)**(0.5/bt-1))**2 +
        A*Tita*4/bt*(0.5/bt-1)*((delta-1)**2)**(0.5/bt-2))
      DeltaBd = b*Delta**(b-1)*Deltad
      DeltaBdd = b*(Delta**(b-1)*Deltadd+(b-1)*Delta**(b-2)*Deltad**2)
      F = exp(-C_*(delta-1)**2-D*(tau-1)**2)
      Fd = -2*C_*F*(delta-1)
      Fdd = 2*C_*F*(2*C_*(delta-1)**2-1)

      B += n*(Delta**b*(F+delta*Fd)+DeltaBd*delta*F)
      C += n*(Delta**b*(2*Fd+delta*Fdd)+2*DeltaBd*(F+delta*Fd) +
          DeltaBdd*delta*F)
    }

    prop = {}
    prop["B"] = B
    prop["C"] = C
    return prop
  }

  void _derivDimensional(self, rho, T)
  {
    /*
    """Calcule the dimensional form or Helmholtz free energy derivatives

    Parameters
    ----------
    rho : float
      Density, [kg/m³]
    T : float
      Temperature, [K]

    Returns
    -------
    prop : dict
      Dictionary with residual helmholtz energy and derivatives:

        * fir, [kJ/kg]
        * firt: ∂fir/∂T|ρ, [kJ/kgK]
        * fird: ∂fir/∂ρ|T, [kJ/m³kg²]
        * firtt: ∂²fir/∂T²|ρ, [kJ/kgK²]
        * firdt: ∂²fir/∂T∂ρ, [kJ/m³kg²K]
        * firdd: ∂²fir/∂ρ²|T, [kJ/m⁶kg]

    References
    ----------
    IAPWS, Guideline on an Equation of State for Humid Air in Contact with
    Seawater and Ice, Consistent with the IAPWS Formulation 2008 for the
    Thermodynamic Properties of Seawater, Table 7,
    http://www.iapws.org/relguide/SeaAir.html
    """
    */
    if(!rho)
    {
      prop = {}
      prop["fir"] = 0
      prop["firt"] = 0
      prop["fird"] = 0
      prop["firtt"] = 0
      prop["firdt"] = 0
      prop["firdd"] = 0
      return prop
    }

    R = self._constants.get("R")/self._constants.get("M", self.M)
    rhoc = self._constants.get("rhoref", self.rhoc)
    Tc = self._constants.get("Tref", self.Tc)
    delta = rho/rhoc
    tau = Tc/T

    ideal = self._phi0(tau, delta)
    fio = ideal["fio"]
    fiot = ideal["fiot"]
    fiott = ideal["fiott"]
    fiod = ideal["fiod"]
    fiodd = ideal["fiodd"]

    res = self._phir(tau, delta)
    fir = res["fir"]
    firt = res["firt"]
    firtt = res["firtt"]
    fird = res["fird"]
    firdd = res["firdd"]
    firdt = res["firdt"]

    prop = {}
    prop["fir"] = R*T*(fio+fir)
    prop["firt"] = R*(fio+fir-(fiot+firt)*tau)
    prop["fird"] = R*T/rhoc*(fiod+fird)
    prop["firtt"] = R*tau**2/T*(fiott+firtt)
    prop["firdt"] = R/rhoc*(fiod+fird-firdt*tau)
    prop["firdd"] = R*T/rhoc**2*(fiodd+firdd)
    return prop
  }

  void _surface(self, T)
  {
    /*
    """Generic equation for the surface tension

    Parameters
    ----------
    T : float
      Temperature, [K]

    Returns
    -------
    σ : float
      Surface tension, [N/m]

    Notes
    -----
    Need a _surf dict in the derived class with the parameters keys:
      sigma: coefficient
      exp: exponent
    """
    */
    tau = 1-T/self.Tc
    sigma = 0
    for(int i = 0; i < self.length; i++)
    {
      n = self._surf["sigma"][i];
      t = self._surf["exp"][i];
      sigma += n*tau**t
    }
    return sigma
  }

  @classmethod
  void _Vapor_Pressure(cls, T)
  {
    /*
    """Auxiliary equation for the vapour pressure

    Parameters
    ----------
    T : float
      Temperature, [K]

    Returns
    -------
    Pv : float
      Vapour pressure, [Pa]

    References
    ----------
    IAPWS, Revised Supplementary Release on Saturation Properties of
    Ordinary Water Substance September 1992,
    http://www.iapws.org/relguide/Supp-sat.html, Eq.1
    """
    */
    Tita = 1-T/cls.Tc
    suma = 0
    for(int i = 0; i < cls.length; i++)
    {
      n = cls._Pv["ao"][i];
      x = cls._Pv["exp"][i];
      suma += n*Tita**x
    }
    Pr = exp(cls.Tc/T*suma)
    Pv = Pr*cls.Pc
    return Pv
  }

  @classmethod
  void _Liquid_Density(cls, T)
  {
    /*
    """Auxiliary equation for the density of saturated liquid

    Parameters
    ----------
    T : float
      Temperature, [K]

    Returns
    -------
    rho : float
      Saturated liquid density, [kg/m³]

    References
    ----------
    IAPWS, Revised Supplementary Release on Saturation Properties of
    Ordinary Water Substance September 1992,
    http://www.iapws.org/relguide/Supp-sat.html, Eq.2
    """
    */
    eq = cls._rhoL["eq"]
    Tita = 1-T/cls.Tc
    if(eq == 2)
    {
      Tita = Tita**(1./3)
    }
    suma = 0
    for(int i = 0; i < cls.length; i++)
    {
      n = cls._rhoL["ao"][i];
      x = cls._rhoL["exp"][i];
      suma += n*Tita**x
    }
    Pr = suma+1
    rho = Pr*cls.rhoc
    return rho
  }

  @classmethod
  void _Vapor_Density(cls, T)
  {
    /*
    """Auxiliary equation for the density of saturated vapor

    Parameters
    ----------
    T : float
      Temperature, [K]

    Returns
    -------
    rho : float
      Saturated vapor density, [kg/m³]

    References
    ----------
    IAPWS, Revised Supplementary Release on Saturation Properties of
    Ordinary Water Substance September 1992,
    http://www.iapws.org/relguide/Supp-sat.html, Eq.3
    """
    */
    eq = cls._rhoG["eq"]
    Tita = 1-T/cls.Tc
    if(eq == 4)
    {
      Tita = Tita**(1./3)
    }
    suma = 0
    for(int i = 0; i < cls.length; i++)
    {
      n = cls._rhoG["ao"][i];
      x = cls._rhoG["exp"][i];
      suma += n*Tita**x
    }
    Pr = exp(suma)
    rho = Pr*cls.rhoc
    return rho
  }

  @classmethod
  void _dPdT_sat(cls, T)
  {
    /*
    """Auxiliary equation for the dP/dT along saturation line

    Parameters
    ----------
    T : float
      Temperature, [K]

    Returns
    -------
    dPdT : float
      dPdT, [MPa/K]

    References
    ----------
    IAPWS, Revised Supplementary Release on Saturation Properties of
    Ordinary Water Substance September 1992,
    http://www.iapws.org/relguide/Supp-sat.html, derived from Eq.1
    """
    */
    Tita = 1-T/cls.Tc
    suma1 = 0
    suma2 = 0
    for(int i = 0; i < cls.length; i++)
    {
      n = cls._Pv["ao"][i];
      x = cls._Pv["exp"][i];
      suma1 -= n*x*Tita**(x-1)/cls.Tc
      suma2 += n*Tita**x
    }
    Pr = (cls.Tc*suma1/T-cls.Tc/T**2*suma2)*exp(cls.Tc/T*suma2)
    dPdT = Pr*cls.Pc
    return dPdT
  }
}

void mainClassDoc()
{
  /*
  """Function decorator used to automatic adiction of base class MEoS in
  subclass __doc__"""
  */
  void decorator(f)
  {
    // __doc__ is only writable in python3.
    // The doc build must be done with python3 so this snnippet do the work
    py_version = platform.python_version()
    if(py_version[0] == "3")
    {
      doc = f.__doc__.split(os.linesep)
      try:
      {
        ind = doc.index("")
      }
      except ValueError:
      {
        ind = 1
      }

      doc1 = os.linesep.join(doc[:ind])
      doc3 = os.linesep.join(doc[ind:])
      doc2 = os.linesep.join(MEoS.__doc__.split(os.linesep)[3:])

      f.__doc__ = doc1 + os.linesep + os.linesep + \
        doc2 + os.linesep + os.linesep + doc3
    }
    return f
  }
  return decorator
}

@mainClassDoc()
class IAPWS95(MEoS):
{
  /*
  """Implementation of IAPWS Formulation 1995 for ordinary water substance,
  (revised release of 2016), for internal procedures, see MEoS base class

  Examples
  --------
  >>> water=IAPWS95(T=300, rho=996.5560)
  >>> water.P, water.cv, water.w, water.s
  0.0992418350 4.13018112 1501.51914 0.393062643

  >>> water=IAPWS95(T=500, rho=0.435)
  >>> water.P, water.cv, water.w, water.s
  0.0999679423 1.50817541 548.31425 7.944882714

  >>> water=IAPWS95(T=900., P=700)
  >>> water.rho, water.cv, water.w, water.s
  870.7690 2.66422350 2019.33608 4.17223802

  >>> water=IAPWS95(T=300., P=0.1)
  >>> water.P, water.rho, water.h, water.s, water.cp, water.w, water.virialB
  0.10000 996.56 112.65 0.39306 4.1806 1501.5 -0.066682

  >>> water=IAPWS95(T=500., P=0.1)
  >>> water.P, water.rho, water.h, water.s, water.cp, water.w, water.virialB
  0.10000 0.43514 2928.6 7.9447 1.9813 548.31 -0.0094137

  >>> water=IAPWS95(T=450., x=0.5)
  >>> water.T, water.P, water.rho, water.h, water.s, water.virialB
  450.00 0.93220 9.5723 1761.8 4.3589 -0.013028

  >>> water=IAPWS95(P=1.5, rho=1000.)
  >>> water.T, water.rho, water.h, water.s, water.cp, water.w, water.virialB
  286.44 1000.0 57.253 0.19931 4.1855 1462.1 -0.085566

  >>> water=IAPWS95(h=3000, s=8.)
  >>> water.T, water.P, water.h, water.s, water.cp, water.w, water.virialB
  536.24 0.11970 3000.0 8.0000 1.9984 567.04 -0.0076606

  >>> water=IAPWS95(h=150, s=0.4)
  >>> water.T, water.P, water.rho, water.h, water.s, water.cp, water.w
  301.27 35.50549 1011.48 150.00 0.40000 4.0932 1564.1

  >>> water=IAPWS95(T=450., rho=300)
  >>> water.T, water.P, water.rho, water.h, water.s, water.x, water.virialB
  450.00 0.93220 300.00 770.82 2.1568 0.010693 -0.013028

  >>> water=IAPWS95(rho=300., P=0.1)
  >>> water.T, water.P, water.rho, water.h, water.s, water.x, water.virialB
  372.76 0.10000 300.00 420.56 1.3110 0.0013528 -0.025144

  >>> water=IAPWS95(h=1500., P=0.1)
  >>> water.T, water.P, water.rho, water.h, water.s, water.x, water.virialB
  372.76 0.10000 1.2303 1500.0 4.2068 0.47952 -0.025144

  >>> water=IAPWS95(s=5., P=3.5)
  >>> water.T, water.P, water.rho, water.h, water.s, water.x, water.virialB
  515.71 3.5000 25.912 2222.8 5.0000 0.66921 -0.0085877

  >>> water=IAPWS95(T=500., u=900)
  >>> water.P, water.rho, water.u, water.h, water.s, water.cp, water.w
  108.21 903.62 900.00 1019.8 2.4271 4.1751 1576.0

  >>> water=IAPWS95(P=0.3, u=1550.)
  >>> water.T, water.P, water.rho, water.u, water.h, water.s, water.x
  406.67 0.30000 3.3029 1550.0 1640.8 4.3260 0.49893

  >>> water=IAPWS95(rho=300, h=1000.)
  >>> water.T, water.P, water.rho, water.u, water.h, water.s, water.x
  494.92 2.3991 300.00 992.00 1000.0 2.6315 0.026071

  >>> water=IAPWS95(rho=30, s=8.)
  >>> water.T, water.P, water.rho, water.u, water.h, water.s, water.cp
  1562.42 21.671 30.000 4628.5 5350.9 8.0000 2.7190

  >>> water=IAPWS95(rho=30, s=4.)
  >>> water.T, water.P, water.rho, water.u, water.h, water.s, water.x
  495.00 2.4029 30.000 1597.3 1677.4 4.0000 0.39218

  >>> water=IAPWS95(rho=300, u=1000.)
  >>> water.T, water.P, water.rho, water.u, water.h, water.s, water.x
  496.44 2.4691 300.000 1000.0 1008.2 2.6476 0.02680

  >>> water=IAPWS95(s=3., h=1000.)
  >>> water.T, water.P, water.rho, water.u, water.h, water.s, water.x
  345.73 0.034850 0.73526 952.60 1000.0 3.0000 0.29920

  >>> water=IAPWS95(u=995., h=1000.)
  >>> water.T, water.P, water.rho, water.u, water.h, water.s, water.x
  501.89 2.7329 546.58 995.00 1000.0 2.6298 0.00866

  >>> water=IAPWS95(u=1000., s=3.)
  >>> water.T, water.P, water.rho, water.u, water.h, water.s, water.x
  371.24 0.094712 1.99072 1000.00 1047.6 3.0000 0.28144

  References
  ----------
  IAPWS, Revised Release on the IAPWS Formulation 1995 for the Thermodynamic
  Properties of Ordinary Water Substance for General and Scientific Use,
  September 2016, http://www.iapws.org/relguide/IAPWS-95.html

  IAPWS, Revised Supplementary Release on Saturation Properties of Ordinary
  Water Substance September 1992, http://www.iapws.org/relguide/Supp-sat.html

  IAPWS, Guideline on a Low-Temperature Extension of the IAPWS-95 Formulation
  for Water Vapor, http://www.iapws.org/relguide/LowT.html

  IAPWS, Revised Advisory Note No. 3: Thermodynamic Derivatives from IAPWS
  Formulations, http://www.iapws.org/relguide/Advise3.pdf
  """
  */
  name = "water"
  CASNumber = "7732-18-5"
  formula = "H2O"
  synonym = "R-718"
  Tc = Tc
  rhoc = rhoc
  Pc = Pc
  M = M
  Tt = 273.16
  Tb = 373.1243
  f_acent = 0.3443
  momentoDipolar = 1.855

  Fi0 = {
   "ao_log": [1, 3.00632],
   "pow": [0, 1],
   "ao_pow": [-8.3204464837497, 6.6832105275932],
   "ao_exp": [0.012436, 0.97315, 1.2795, 0.96956, 0.24873],
   "titao": [1.28728967, 3.53734222, 7.74073708, 9.24437796, 27.5075105]
  }

  _constants =
  {
    "R": 8.314371357587,

    "nr1": [0.12533547935523e-1, 0.78957634722828e1, -0.87803203303561e1, 0.31802509345418, -0.26145533859358, -0.78199751687981e-2, 0.88089493102134e-2],
    "d1": [1, 1, 1, 2, 2, 3, 4],
    "t1": [-0.5, 0.875, 1, 0.5, 0.75, 0.375, 1],

    "nr2": [-0.66856572307965, 0.20433810950965, -0.66212605039687e-4, -0.19232721156002, -0.25709043003438, 0.16074868486251, -0.4009282892587e-1, 0.39343422603254e-6, -0.75941377088144e-5, 0.56250979351888e-3, -0.15608652257135e-4, 0.11537996422951e-8, .36582165144204e-6, -.13251180074668e-11, -.62639586912454e-9, -0.10793600908932, 0.17611491008752e-1, 0.22132295167546, -0.40247669763528, 0.58083399985759, 0.49969146990806e-2, -0.31358700712549e-1, -0.74315929710341, 0.47807329915480, 0.20527940895948e-1, -0.13636435110343, 0.14180634400617e-1, 0.83326504880713e-2, -0.29052336009585e-1, 0.38615085574206e-1, -0.20393486513704e-1, -0.16554050063734e-2, .19955571979541e-2, 0.15870308324157e-3, -0.16388568342530e-4, 0.43613615723811e-1, 0.34994005463765e-1, -0.76788197844621e-1, 0.22446277332006e-1, -0.62689710414685e-4, -0.55711118565645e-9, -0.19905718354408, 0.31777497330738, -0.11841182425981],
    "c2": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 6, 6, 6, 6],
    "d2": [1, 1, 1, 2, 2, 3, 4, 4, 5, 7, 9, 10, 11, 13, 15, 1, 2, 2, 2, 3, 4, 4, 4, 5, 6, 6, 7, 9, 9, 9, 9, 9, 10, 10, 12, 3, 4, 4, 5, 14, 3, 6, 6, 6],
    "t2": [4, 6, 12, 1, 5, 4, 2, 13, 9, 3, 4, 11, 4, 13, 1, 7, 1, 9, 10, 10, 3, 7, 10, 10, 6, 10, 10, 1, 2, 3, 4, 8, 6, 9, 8, 16, 22, 23, 23, 10, 50, 44, 46, 50],
    "gamma2": [1]*44,

    "nr3": [-0.31306260323435e2, 0.31546140237781e2, -0.25213154341695e4],
    "d3": [3]*3,
    "t3": [0, 1, 4],
    "alfa3": [20]*3,
    "beta3": [150, 150, 250],
    "gamma3": [1.21, 1.21, 1.25],
    "epsilon3": [1.]*3,

    "nr4": [-0.14874640856724, 0.31806110878444],
    "a4": [3.5, 3.5],
    "b4": [0.85, 0.95],
    "B": [0.2, 0.2],
    "C": [28, 32],
    "D": [700, 800],
    "A": [0.32, .32],
    "beta4": [0.3, 0.3]
  }

  _Pv =
  {
    "ao": [-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502],
    "exp": [1, 1.5, 3, 3.5, 4, 7.5]
  }

  _rhoL =
  {
    "eq": 2,
    "ao": [1.99274064, 1.09965342, -0.510839303, -1.75493479, -45.5170352, -6.74694450e5],
    "exp": [1, 2, 5, 16, 43, 110]
  }

  _rhoG =
  {
    "eq": 4,
    "ao": [-2.0315024, -2.6830294, -5.38626492, -17.2991605, -44.7586581, -63.9201063],
    "exp": [1, 2, 4, 9, 18.5, 35.5]
  }

  void _phi0(self, tau, delta)
  {
    //Low temperature extension of the IAPWS-95
    prop = MEoS._phi0(self, tau, delta)

    T = self.Tc/tau
    if(50 <= T < 130)
    {
      fex, fext, fextt = self._phiex(T)
      prop["fio"] += fex
      prop["fiot"] += fext
      prop["fiott"] += fextt
    }
    return prop
  }

  void _phiex(self, T)
  {
    //Low temperature extension
    tau = self.Tc/T
    E = 0.278296458178592
    ep = self.Tc/130
    fex = E*(-1/2/tau-3/ep**2*(tau+ep)*log(tau/ep)-9/2/ep+9*tau/2/ep**2 +
         tau**2/2/ep**3)
    fext = E*(1/2/tau**2-3/tau/ep-3/ep**2*log(tau/ep)+3/2/ep**2+tau/ep**3)
    fextt = E*(-1/tau+1/ep)**3
    return fex, fext, fextt
  }

  @classmethod
  void _alfa_sat(cls, T)
  {
    /*
    """Auxiliary equation for the alfa coefficient for calculate the
    enthalpy along the saturation line

    Parameters
    ----------
    T : float
      Temperature, [K]

    Returns
    -------
    alfa : float
      alfa coefficient, [kJ/kg]

    References
    ----------
    IAPWS, Revised Supplementary Release on Saturation Properties of
    Ordinary Water Substance September 1992,
    http://www.iapws.org/relguide/Supp-sat.html, Eq.4
    """
    */
    di = [-1135.905627715, -5.65134998e-8, 2690.66631, 127.287297, -135.003439, 0.981825814]
    expi = [0, -19, 1, 4.5, 5, 54.5]
    Tita = T/cls.Tc
    alfa = 0
    for(int i = 0; i < di.length; i++)
    {
      d = di[i];
      x = expi[i];
      alfa += d*Tita**x
    }
    return alfa
  }

  @classmethod
  void _phi_sat(cls, T)
  {
    /*
    """Auxiliary equation for the phi coefficient for calculate the
    entropy along the saturation line

    Parameters
    ----------
    T : float
      Temperature, [K]

    Returns
    -------
    phi : float
      phi coefficient, [kJ/kgK]

    References
    ----------
    IAPWS, Revised Supplementary Release on Saturation Properties of
    Ordinary Water Substance September 1992,
    http://www.iapws.org/relguide/Supp-sat.html, Eq.5
    """
    */
    di = [2319.5246, -5.65134998e-8*19/20, 2690.66631, 127.287297*9/7, -135.003439*5/4, 0.981825814*109/107]
    expi = [0, -20, None, 3.5, 4, 53.5]
    Tita = T/cls.Tc
    suma = 0
    for(int i = 0; i < di.length; i++)
    {
      d = di[i];
      x = expi[i];
      if(x is None)
      {
        suma += d*log(Tita)
      }
      else
      {
        suma += d*Tita**x
      }
    }
    phi = suma/cls.Tc
    return phi
  }

  @classmethod
  void _Liquid_Enthalpy(cls, T)
  {
    /*
    """Auxiliary equation for the specific enthalpy for saturated liquid

    Parameters
    ----------
    T : float
      Temperature, [K]

    Returns
    -------
    h : float
      Saturated liquid enthalpy, [kJ/kg]

    References
    ----------
    IAPWS, Revised Supplementary Release on Saturation Properties of
    Ordinary Water Substance September 1992,
    http://www.iapws.org/relguide/Supp-sat.html, Eq.6
    """
    */
    alfa = cls._alfa_sat(T)
    rho = cls._Liquid_Density(T)
    dpdT = cls._dPdT_sat(T)
    h = alfa+T/rho*dpdT*1000
    return h
  }

  @classmethod
  void _Vapor_Enthalpy(cls, T)
  {
    /*
    """Auxiliary equation for the specific enthalpy for saturated vapor

    Parameters
    ----------
    T : float
      Temperature, [K]

    Returns
    -------
    h : float
      Saturated vapor enthalpy, [kJ/kg]

    References
    ----------
    IAPWS, Revised Supplementary Release on Saturation Properties of
    Ordinary Water Substance September 1992,
    http://www.iapws.org/relguide/Supp-sat.html, Eq.7
    """
    */
    alfa = cls._alfa_sat(T)
    rho = cls._Vapor_Density(T)
    dpdT = cls._dPdT_sat(T)
    h = alfa+T/rho*dpdT*1000
    return h
  }

  @classmethod
  void _Liquid_Entropy(cls, T)
  {
    /*
    """Auxiliary equation for the specific entropy for saturated liquid

    Parameters
    ----------
    T : float
      Temperature, [K]

    Returns
    -------
    s : float
      Saturated liquid entropy, [kJ/kgK]

    References
    ----------
    IAPWS, Revised Supplementary Release on Saturation Properties of
    Ordinary Water Substance September 1992,
    http://www.iapws.org/relguide/Supp-sat.html, Eq.8
    """
    */
    phi = cls._phi_sat(T)
    rho = cls._Liquid_Density(T)
    dpdT = cls._dPdT_sat(T)
    s = phi+dpdT/rho*1000
    return s
  }

  @classmethod
  void _Vapor_Entropy(cls, T)
  {
    /*
    """Auxiliary equation for the specific entropy for saturated vapor

    Parameters
    ----------
    T : float
      Temperature, [K]

    Returns
    -------
    s : float
      Saturated liquid entropy, [kJ/kgK]

    References
    ----------
    IAPWS, Revised Supplementary Release on Saturation Properties of
    Ordinary Water Substance September 1992,
    http://www.iapws.org/relguide/Supp-sat.html, Eq.9
    """
    */
    phi = cls._phi_sat(T)
    rho = cls._Vapor_Density(T)
    dpdT = cls._dPdT_sat(T)
    s = phi+dpdT/rho*1000
    return s
  }

  void _visco(self, rho, T, fase)
  {
    ref = IAPWS95()
    st = ref._Helmholtz(rho, 1.5*Tc)
    delta = rho/rhoc
    drho = 1e3/self.R/1.5/Tc/(1+2*delta*st["fird"]+delta**2*st["firdd"])
    return _Viscosity(rho, T, fase, drho)
  }

  void _thermo(self, rho, T, fase)
  {
    ref = IAPWS95()
    st = ref._Helmholtz(rho, 1.5*Tc)
    delta = rho/rhoc
    drho = 1e3/self.R/1.5/Tc/(1+2*delta*st["fird"]+delta**2*st["firdd"])
    return _ThCond(rho, T, fase, drho)
  }

  void _surface(self, T)
  {
    s = _Tension(T)
    return s
  }
}

class IAPWS95_PT(IAPWS95):
{
  //Derivated class for direct P and T input
  void __init__(self, P, T)
  {
    IAPWS95.__init__(self, T=T, P=P)
  }
}

class IAPWS95_Ph(IAPWS95):
{
  //Derivated class for direct P and h input
  void __init__(self, P, h)
  {
    IAPWS95.__init__(self, P=P, h=h)
  }
}

class IAPWS95_Ps(IAPWS95):
{
  //Derivated class for direct P and s input
  void __init__(self, P, s)
  {
    IAPWS95.__init__(self, P=P, s=s)
  }
}

class IAPWS95_Px(IAPWS95):
{
  //Derivated class for direct P and v input
  void __init__(self, P, x)
  {
    IAPWS95.__init__(self, P=P, x=x)
  }
}

class IAPWS95_Tx(IAPWS95):
{
  //Derivated class for direct T and x input
  void __init__(self, T, x)
  {
    IAPWS95.__init__(self, T=T, x=x)
  }
}

@mainClassDoc()
class D2O(MEoS):
{
  /*
  """Implementation of IAPWS Formulation for heavy water substance,
  for internal procedures, see MEoS base class

  Examples
  --------
  >>> hwater=D2O(T=300, rho=996.5560)
  >>> hwater.P, hwater.Liquid.cv, hwater.Liquid.w
  0.0030675947 4.21191157 5332.04871

  References
  ----------
  IAPWS, Release on the IAPWS Formulation 2017 for the Thermodynamic
  Properties of Heavy Water, http://www.iapws.org/relguide/Heavy-2017.pdf
  IAPWS, Revised Advisory Note No. 3: Thermodynamic Derivatives from IAPWS
  Formulations, http://www.iapws.org/relguide/Advise3.pdf
  """
  */
  name = "heavy water"
  CASNumber = "7789-20-0"
  formula = "D2O"
  synonym = "deuterium oxide"
  Tc = Tc_D2O
  rhoc = rhoc_D2O
  Pc = Pc_D2O
  M = 20.027508  // g/mol
  Tt = 276.97
  Tb = 374.563
  f_acent = 0.364
  momentoDipolar = 1.9

  Fi0 =
  {
   "ao_log": [1, 3],
   "pow": [0, 1],
   "ao_pow": [-8.670994022646, 6.96033578458778],
   "ao_exp": [0.010633, 0.99787, 2.1483, 0.3549],
   "titao": [308/Tc, 1695/Tc, 3949/Tc, 10317/Tc],
   "ao_hyp": [],
   "hyp": []
  }

  _constants =
  {
    "R": 8.3144598,

    "nr1": [0.122082060e-1, 0.296956870e1, -0.379004540e1, 0.941089600, -0.922466250, -0.139604190e-1],
    "d1": [4, 1, 1, 2, 2, 3],
    "t1": [1.0000, 0.6555, 0.9369, 0.5610, 0.7017, 1.0672],

    "nr2": [-0.125203570, -0.555391500e1, -0.493009740e1, -0.359470240e-1, -0.936172870e1, -0.691835150],
    "c2": [1, 2, 2, 1, 2, 2],
    "d2": [1, 1, 3, 2, 2, 1],
    "t2": [3.9515, 4.6000, 5.1590, 0.2000, 5.4644, 2.3660],
    "gamma2": [1]*6,

    "nr3": [-0.456110600e-1, -0.224513300e1, 0.860006070e1, -0.248410420e1, 0.164476900e2, 0.270393360e1, 0.375637470e2, -0.177607760e1, 0.220924640e1, 0.519652000e1, 0.421097400, -0.391921100],
    "t3": [3.4553, 1.4150, 1.5745, 3.4540, 3.8106, 4.8950, 1.4300, 1.5870, 3.7900, 2.6200, 1.9000, 4.3200],
    "d3": [1, 3, 1, 3, 1, 1, 2, 2, 2, 1, 1, 1],
    "alfa3": [0.6014, 1.4723, 1.5305, 2.4297, 1.3086, 1.3528, 3.4456, 1.2645, 2.5547, 1.2148, 18.738, 18.677],
    "beta3": [0.4200, 2.4318, 1.2888, 8.2710, 0.3673, 0.9504, 7.8318, 3.3281, 7.1753, 0.9465, 1177.0, 1167.0],
    "epsilon3": [1.8663, 0.2895, 0.5803, 0.2236, 0.6815, 0.9495, 1.1158, 0.1607, 0.4144, 0.9683, 0.9488, 0.9487],
    "gamma3": [1.5414, 1.3794, 1.7385, 1.3045, 2.7242, 3.5321, 2.4552, 0.8319, 1.3500, 2.5617, 1.0491, 1.0486]
  }

  _Pv =
  {
    "ao": [-0.80236e1, 0.23957e1, -0.42639e2, 0.99569e2, -0.62135e2],
    "exp": [1.0, 1.5, 2.75, 3.0, 3.2]
  }

  _rhoL = {
    "eq": 1,
    "ao": [0.26406e1, 0.97090e1, -0.18058e2, 0.87202e1, -0.74487e1],
    "exp": [0.3678, 1.9, 2.2, 2.63, 7.3]
  }

  _rhoG =
  {
    "eq": 3,
    "ao": [-0.37651e1, -0.38673e2, 0.73024e2, -0.13251e3, 0.75235e2, -0.70412e2],
    "exp": [0.409, 1.766, 2.24, 3.04, 3.42, 6.9]
  }

  void _visco(self, rho, T, fase)
  {
    return _D2O_Viscosity(rho, T)
  }

  void _thermo(self, rho, T, fase)
  {
    return _D2O_ThCond(rho, T)
  }

  void _surface(self, T)
  {
    s = _D2O_Tension(T)
    return s
  }
}

