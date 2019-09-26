using System;
using System.Collections.Generic;
using System.Linq;

public static class _iapws
{

  static _iapws()
  {
    // Miscelaneous IAPWS standards. This module include:
    //
    // * :func:`_Ice`: Ice Ih state equation
    // * :func:`_Liquid`: Properties of liquid water at 0.1 MPa
    // * :func:`_Supercooled`: Thermodynamic properties of supercooled water
    // * :func:`_Sublimation_Pressure`: Sublimation pressure correlation
    // * :func:`_Melting_Pressure`: Melting pressure correlation
    // * :func:`_Viscosity`: Viscosity correlation
    // * :func:`_ThCond`: Themal conductivity correlation
    // * :func:`_Tension`: Surface tension correlation
    // * :func:`_Dielectric`: Dielectric constant correlation
    // * :func:`_Refractive`: Refractive index correlation
    // * :func:`_Kw`: Ionization constant correlation for ordinary water
    // * :func:`_Conductivity`: Electrolytic conductivity correlation
    // heavy water
    // * :func:`_Henry`: Henry constant for liquid-gas equilibrium
    // * :func:`_Kvalue`: Vapor-liquid distribution constant
  }

  public static double M = 18.015268;
  public static double R = 0.461526;
  public static double Tc = 647.096;
  public static double Pc = 22.064;
  public static double rhoc = 322.0;
  public static double Tt = 273.16;
  public static double Pt = 0.000611657;
  public static double Tb = 373.1243;
  public static double f_acent = 0.3443;
  public static double Dipole = 1.85498;

  // IAPWS-06 for Ice
  // Basic state equation for Ice Ih
  //
  //   Parameters
  //   ----------
  //   T : float
  //     Temperature, [K]
  //   P : float
  //     Pressure, [MPa]
  //
  //   Returns
  //   -------
  //   prop : dict
  //     Dict with calculated properties of ice. The available properties are:
  //
  //       * rho: Density, [kg/m³]
  //       * h: Specific enthalpy, [kJ/kg]
  //       * u: Specific internal energy, [kJ/kg]
  //       * a: Specific Helmholtz energy, [kJ/kg]
  //       * g: Specific Gibbs energy, [kJ/kg]
  //       * s: Specific entropy, [kJ/kgK]
  //       * cp: Specific isobaric heat capacity, [kJ/kgK]
  //       * alfav: Cubic expansion coefficient, [1/K]
  //       * beta: Pressure coefficient, [MPa/K]
  //       * xkappa: Isothermal compressibility, [1/MPa]
  //       * ks: Isentropic compressibility, [1/MPa]
  //       * gt: [∂g/∂T]P
  //       * gtt: [∂²g/∂T²]P
  //       * gp: [∂g/∂P]T
  //       * gpp: [∂²g/∂P²]T
  //       * gtp: [∂²g/∂T∂P]
  //
  //   Notes
  //   ------
  //   Raise :class:`NotImplementedError` if input isn't in limit:
  //
  //     * T ≤ 273.16
  //     * P ≤ 208.566
  //     * State below the melting and sublimation lines
  //
  //   Examples
  //   --------
  //   >>> st1 = _Ice(100, 100)
  //   >>> st1["rho"], st1["h"], st1["s"]
  //   941.678203297 -483.491635676 -2.61195122589
  //
  //   >>> st2 = _Ice(273.152519,0.101325)
  //   >>> st2["a"], st2["u"], st2["cp"]
  //   -0.00918701567 -333.465403393 2.09671391024
  //
  //   >>> st3 = _Ice(273.16,611.657e-6)
  //   >>> st3["alfav"], st3["beta"], st3["xkappa"], st3["ks"]
  //   0.000159863102566 1.35714764659 1.17793449348e-04 1.14161597779e-04
  //
  //   References
  //   ----------
  //   IAPWS, Revised Release on the Equation of State 2006 for H2O Ice Ih
  //   September 2009, http://iapws.org/relguide/Ice-2009.html
  //

  public static Dictionary<string, double> _Ice(double T, double P)
  {
    // Check input in range of validity
    if(T > 273.16)
    {
      // No Ice Ih stable
      Debug.Log("Metastable ice");
    }
    else if(P > 208.566)
    {
      // Ice Ih limit upper pressure
      throw new NotImplementedException("Incoming out of bound");
    }
    else if(P < Pt)
    {
      var Psub = _Sublimation_Pressure(T);
      if(Psub > P)
      {
        // Zone Gas
        Debug.Log("Metastable ice in vapor region");
      }
    }
    else if(251.165 < T)
    {
      var Pmel = _Melting_Pressure(T);
      if(Pmel < P)
      {
        // Zone Liquid
        Debug.Log("Metastable ice in liquid region");
      }
    }
    var Tr = T / Tt;
    var Pr = P / Pt;
    var P0 = 0.101325 / Pt;
    var s0 = -3327.33756492168 * 0.001;
    var gok = new List<double> { -632020.233335886, 0.655022213658955, -1.89369929326131E-08, 3.39746123271053E-15, -5.56464869058991E-22 };
    var r2k = new List<object>
    {
      new Complex(-72.597457432922, -78.100842711287) * 0.001,
      new Complex(-5.57107698030123E-05, 4.64578634580806E-05) * 0.001,
      new Complex(2.34801409215913E-11, -2.85651142904972E-11) * 0.001
    };
    var t1 = new Complex(0.0368017112855051, 0.0510878114959572);
    var t2 = new Complex(0.337315741065416, 0.335449415919309);
    var r1 = new Complex(44.7050716285388, 65.6876847463481) * 0.001;
    double go = 0;
    foreach (var k in Enumerable.Range(0, 5))
    {
      go += gok[k] * 0.001 * Math.Pow(Pr - P0, k);
    }
    double gop = 0;
    foreach (var k in Enumerable.Range(1, 5 - 1))
    {
      gop += gok[k] * 0.001 * k / Pt * Math.Pow(Pr - P0, k - 1);
    }
    double gopp = 0;
    foreach (var k in Enumerable.Range(2, 5 - 2))
    {
      gopp += gok[k] * 0.001 * k * (k - 1) / Math.Pow(Pt, 2) * Math.Pow(Pr - P0, k - 2);
    }
    double r2 = 0;
    foreach (var k in Enumerable.Range(0, 3))
    {
      r2 += r2k[k] * Math.Pow(Pr - P0, k);
    }
    double r2p = 0;
    foreach (var k in Enumerable.Range(1, 3 - 1))
    {
      r2p += r2k[k] * k / Pt * Math.Pow(Pr - P0, k - 1);
    }
    double r2pp = r2k[2] * 2 / Math.Pow(Pt, 2);
    double c = r1 * ((t1 - Tr) * log_c(t1 - Tr) + (t1 + Tr) * log_c(t1 + Tr) - 2 * t1 * log_c(t1) - Math.Pow(Tr, 2) / t1) + r2 * ((t2 - Tr) * log_c(t2 - Tr) + (t2 + Tr) * log_c(t2 + Tr) - 2 * t2 * log_c(t2) - Math.Pow(Tr, 2) / t2);
    double ct = r1 * (-log_c(t1 - Tr) + log_c(t1 + Tr) - 2 * Tr / t1) + r2 * (-log_c(t2 - Tr) + log_c(t2 + Tr) - 2 * Tr / t2);
    double ctt = r1 * (1 / (t1 - Tr) + 1 / (t1 + Tr) - 2 / t1) + r2 * (1 / (t2 - Tr) + 1 / (t2 + Tr) - 2 / t2);
    double cp = r2p * ((t2 - Tr) * log_c(t2 - Tr) + (t2 + Tr) * log_c(t2 + Tr) - 2 * t2 * log_c(t2) - Math.Pow(Tr, 2) / t2);
    double ctp = r2p * (-log_c(t2 - Tr) + log_c(t2 + Tr) - 2 * Tr / t2);
    double cpp = r2pp * ((t2 - Tr) * log_c(t2 - Tr) + (t2 + Tr) * log_c(t2 + Tr) - 2 * t2 * log_c(t2) - Math.Pow(Tr, 2) / t2);
    double g = go - s0 * Tt * Tr + Tt * c.real;
    double gt = -s0 + ct.real;
    double gp = gop + Tt * cp.real;
    double gtt = ctt.real / Tt;
    double gtp = ctp.real;
    double gpp = gopp + Tt * cpp.real;

    var propiedades = new Dictionary<string, double> { };
    propiedades["gt"] = gt;
    propiedades["gp"] = gp;
    propiedades["gtt"] = gtt;
    propiedades["gpp"] = gpp;
    propiedades["gtp"] = gtp;
    propiedades["T"] = T;
    propiedades["P"] = P;
    propiedades["v"] = gp / 1000;
    propiedades["rho"] = 1000.0 / gp;
    propiedades["h"] = g - T * gt;
    propiedades["s"] = -gt;
    propiedades["cp"] = -T * gtt;
    propiedades["u"] = g - T * gt - P * gp;
    propiedades["g"] = g;
    propiedades["a"] = g - P * gp;
    propiedades["alfav"] = gtp / gp;
    propiedades["beta"] = -gtp / gpp;
    propiedades["xkappa"] = -gpp / gp;
    propiedades["ks"] = (Math.Pow(gtp, 2) - gtt * gpp) / gp / gtt;

    return propiedades;
  }

  // IAPWS-08 for Liquid water at 0.1 MPa
  // Supplementary release on properties of liquid water at 0.1 MPa
  //
  //   Parameters
  //   ----------
  //   T : float
  //     Temperature, [K]
  //   P : float
  //     Pressure, [MPa]
  //     Although this relation is for P=0.1MPa, can be extrapoled at pressure
  //     0.3 MPa
  //
  //   Returns
  //   -------
  //   prop : dict
  //     Dict with calculated properties of water. The available properties are:
  //
  //       * h: Specific enthalpy, [kJ/kg]
  //       * u: Specific internal energy, [kJ/kg]
  //       * a: Specific Helmholtz energy, [kJ/kg]
  //       * g: Specific Gibbs energy, [kJ/kg]
  //       * s: Specific entropy, [kJ/kgK]
  //       * cp: Specific isobaric heat capacity, [kJ/kgK]
  //       * cv: Specific isochoric heat capacity, [kJ/kgK]
  //       * w: Speed of sound, [m/s²]
  //       * rho: Density, [kg/m³]
  //       * v: Specific volume, [m³/kg]
  //       * vt: [∂v/∂T]P, [m³/kgK]
  //       * vtt: [∂²v/∂T²]P, [m³/kgK²]
  //       * vp: [∂v/∂P]T, [m³/kg/MPa]
  //       * vtp: [∂²v/∂T∂P], [m³/kg/MPa]
  //       * alfav: Cubic expansion coefficient, [1/K]
  //       * xkappa : Isothermal compressibility, [1/MPa]
  //       * ks: Isentropic compressibility, [1/MPa]
  //       * mu: Viscosity, [mPas]
  //       * k: Thermal conductivity, [W/mK]
  //       * epsilon: Dielectric constant, [-]
  //
  //   Notes
  //   ------
  //   Raise :class:`NotImplementedError` if input isn't in limit:
  //
  //     * 253.15 ≤ T ≤ 383.15
  //     * 0.1 ≤ P ≤ 0.3
  //
  //   Examples
  //   --------
  //   >>> st1 = _Liquid(260)
  //   >>> st1["rho"], st1["h"], st1["s"]
  //   997.0683602710492 -55.86223174460868 -0.20998554842619535
  //
  //   References
  //   ----------
  //   IAPWS, Revised Supplementary Release on Properties of Liquid Water at 0.1
  //   MPa, http://www.iapws.org/relguide/LiquidWater.html
  //

  public static Dictionary<string, double> _Liquid(double T, double P = 0.1)
  {
    // Check input in range of validity
    if(T <= 253.15 || T >= 383.15 || P < 0.1 || P > 0.3)
    {
      throw new NotImplementedException("Incoming out of bound");
    }
    else if(P != 0.1)
    {
      // Raise a warning if the P value is extrapolated
      Debug.Log("Using extrapolated values");
    }

    double UNUSED = 0.0;
    var R = 0.46151805;
    var Po = 0.1;
    var Tr = 10;
    var tau = T / Tr;
    var alfa = Tr / (593 - T);
    var beta = Tr / (T - 232);
    var a = new List<double>
    {
      UNUSED,
      -166147.0539,
      2708781.64,
      -155719154.4,
      UNUSED,
      0.0193763157,
      6744.58446,
      -222521.604,
      100231247.0,
      -1635521180.0,
      8322996580.0,
      -7.5245878E-06,
      -0.013767418,
      10.627293,
      -204.57795,
      1203.7414
    };
    var b = new List<double>
    {
      UNUSED,
      -0.8237426256,
      1.908956353,
      -2.017597384,
      0.8546361348,
      0.00578545292,
      -0.0153195665,
      0.0311337859,
      -0.0423546241,
      0.0338713507,
      -0.0119946761,
      -3.109147E-06,
      2.8964919E-05,
      -0.00013112763,
      0.00030410453,
      -0.00039034594,
      0.00023403117,
      -4.8510101E-05
    };
    var c = new List<double> { UNUSED, -245.2093414, 38.69269598, -8.983025854 };
    var n = new List<double>
    {
      UNUSED,
      4,
      5,
      7,
      UNUSED,
      UNUSED,
      4,
      5,
      7,
      8,
      9,
      1,
      3,
      5,
      6,
      7
    };
    var m = new List<double>
    {
      UNUSED,
      2,
      3,
      4,
      5,
      1,
      2,
      3,
      4,
      5,
      6,
      1,
      3,
      4,
      5,
      6,
      7,
      9
    };

    var suma1 = (from i in Enumerable.Range(1, 4 - 1) select (a[i] * Math.Pow(alfa, n[i]))).ToList().Sum();
    var suma2 = (from i in Enumerable.Range(1, 5 - 1) select (b[i] * Math.Pow(beta, m[i]))).ToList().Sum();
    var go = R * Tr * (c[1] + c[2] * tau + c[3] * tau * Math.Log(tau) + suma1 + suma2);
    suma1 = (from i in Enumerable.Range(6, 11 - 6) select (a[i] * Math.Pow(alfa, n[i]))).ToList().Sum();
    suma2 = (from i in Enumerable.Range(5, 11 - 5) select (b[i] * Math.Pow(beta, m[i]))).ToList().Sum();
    var vo = R * Tr / Po / 1000 * (a[5] + suma1 + suma2);
    suma1 = (from i in Enumerable.Range(11, 16 - 11) select (a[i] * Math.Pow(alfa, n[i]))).ToList().Sum();
    suma2 = (from i in Enumerable.Range(11, 18 - 11) select (b[i] * Math.Pow(beta, m[i]))).ToList().Sum();
    var vpo = R * Tr / Math.Pow(Po, 2) / 1000 * (suma1 + suma2);
    suma1 = (from i in Enumerable.Range(1, 4 - 1) select (n[i] * a[i] * Math.Pow(alfa, n[i] + 1))).ToList().Sum();
    suma2 = (from i in Enumerable.Range(1, 5 - 1) select (m[i] * b[i] * Math.Pow(beta, m[i] + 1))).ToList().Sum();
    var so = -R * (c[2] + c[3] * (1 + Math.Log(tau)) + suma1 - suma2);
    suma1 = (from i in Enumerable.Range(1, 4 - 1) select (n[i] * (n[i] + 1) * a[i] * Math.Pow(alfa, n[i] + 2))).ToList().Sum();
    suma2 = (from i in Enumerable.Range(1, 5 - 1) select (m[i] * (m[i] + 1) * b[i] * Math.Pow(beta, m[i] + 2))).ToList().Sum();
    var cpo = -R * (c[3] + tau * suma1 + tau * suma2);
    suma1 = (from i in Enumerable.Range(6, 11 - 6) select (n[i] * a[i] * Math.Pow(alfa, n[i] + 1))).ToList().Sum();
    suma2 = (from i in Enumerable.Range(5, 11 - 5) select (m[i] * b[i] * Math.Pow(beta, m[i] + 1))).ToList().Sum();
    var vto = R / Po / 1000 * (suma1 - suma2);
    // This properties are only neccessary for computing thermodynamic
    // properties at pressures different from 0.1 MPa
    suma1 = (from i in Enumerable.Range(6, 11 - 6) select (n[i] * (n[i] + 1) * a[i] * Math.Pow(alfa, n[i] + 2))).ToList().Sum();
    suma2 = (from i in Enumerable.Range(5, 11 - 5) select (m[i] * (m[i] + 1) * b[i] * Math.Pow(beta, m[i] + 2))).ToList().Sum();
    var vtto = R / Tr / Po / 1000 * (suma1 + suma2);
    suma1 = (from i in Enumerable.Range(11, 16 - 11) select (n[i] * a[i] * Math.Pow(alfa, n[i] + 1))).ToList().Sum();
    suma2 = (from i in Enumerable.Range(11, 18 - 11) select (m[i] * b[i] * Math.Pow(beta, m[i] + 1))).ToList().Sum();
    var vpto = R / Math.Pow(Po, 2) / 1000 * (suma1 - suma2);

    if(P != 0.1)
    {
      go += vo * (P - 0.1);
      so -= vto * (P - 0.1);
      cpo -= T * vtto * (P - 0.1);
      vo -= vpo * (P - 0.1);
      vto += vpto * (P - 0.1);
      var vppo = 3.24E-10 * R * Tr / Math.Pow(0.1, 3);
      vpo += vppo * (P - 0.1);
    }
    var h = go + T * so;
    var u = h - P * vo;
    a = go - P * vo;
    var cv = cpo + T * Math.Pow(vto, 2) / vpo;
    var xkappa = -vpo / vo;
    alfa = vto / vo;
    var ks = -T * Math.Pow(vto, 2) / cpo + vpo / vo;
    var w = Math.Pow(-Math.Pow(vo, 2) * 1000000000.0 / (vpo * 1000.0 + T * Math.Pow(vto, 2) * 1000000.0 / cpo), 0.5);

    var propiedades = new Dictionary<string, double> { };
    propiedades["g"] = go;
    propiedades["T"] = T;
    propiedades["P"] = P;
    propiedades["v"] = vo;
    propiedades["vt"] = vto;
    propiedades["vp"] = vpo;
    propiedades["vpt"] = vpto;
    propiedades["vtt"] = vtto;
    propiedades["rho"] = 1 / vo;
    propiedades["h"] = h;
    propiedades["s"] = so;
    propiedades["cp"] = cpo;
    propiedades["cv"] = cv;
    propiedades["u"] = u;
    propiedades["a"] = a;
    propiedades["xkappa"] = xkappa;
    propiedades["alfav"] = vto / vo;
    propiedades["ks"] = ks;
    propiedades["w"] = w;

    // Viscosity correlation, Eq 7
    a = new List<double> { UNUSED, 280.68, 511.45, 61.131, 0.45903 };
    b = new List<double> { UNUSED, -1.9, -7.7, -19.6, -40 };
    var T_ = T / 300;
    var mu = (from i in Enumerable.Range(1, 5 - 1)
      select (a[i] * Math.Pow(T_, b[i]))).ToList().Sum() / 1000000.0;
    propiedades["mu"] = mu;
    // Thermal conductivity correlation, Eq 8
    c = new List<double> { UNUSED, 1.663, -1.7781, 1.1567, -0.432115 };
    var d = new List<double> { UNUSED, -1.15, -3.4, -6.0, -7.6 };
    var k = (from i in Enumerable.Range(1, 5 - 1)
      select (c[i] * Math.Pow(T_, d[i]))).ToList().Sum();
    propiedades["k"] = k;
    // Dielectric constant correlation, Eq 9
    var e = new List<double> { UNUSED, -43.7527, 299.504, -399.364, 221.327 };
    var f = new List<double> { UNUSED, -0.05, -1.47, -2.11, -2.31 };
    var epsilon = (from i in Enumerable.Range(1, 5 - 1)
      select (e[i] * Math.Pow(T_, f[i]))).ToList().Sum();
    propiedades["epsilon"] = epsilon;
    return propiedades;
  }

  // IAPWS-15 for supercooled liquid water
  // Guideline on thermodynamic properties of supercooled water
  //
  //   Parameters
  //   ----------
  //   T : float
  //     Temperature, [K]
  //   P : float
  //     Pressure, [MPa]
  //
  //   Returns
  //   -------
  //   prop : dict
  //     Dict with calculated properties of water. The available properties are:
  //
  //       * L: Ordering field, [-]
  //       * x: Mole fraction of low-density structure, [-]
  //       * rho: Density, [kg/m³]
  //       * s: Specific entropy, [kJ/kgK]
  //       * h: Specific enthalpy, [kJ/kg]
  //       * u: Specific internal energy, [kJ/kg]
  //       * a: Specific Helmholtz energy, [kJ/kg]
  //       * g: Specific Gibbs energy, [kJ/kg]
  //       * alfap: Thermal expansion coefficient, [1/K]
  //       * xkappa : Isothermal compressibility, [1/MPa]
  //       * cp: Specific isobaric heat capacity, [kJ/kgK]
  //       * cv: Specific isochoric heat capacity, [kJ/kgK]
  //       * w: Speed of sound, [m/s²]
  //
  //   Notes
  //   ------
  //   Raise :class:`NotImplementedError` if input isn't in limit:
  //
  //     * Tm ≤ T ≤ 300
  //     * 0 < P ≤ 1000
  //
  //   The minimum temperature in range of validity is the melting temperature, it
  //   depend of pressure
  //
  //   Examples
  //   --------
  //   >>> liq = _Supercooled(235.15, 0.101325)
  //   >>> liq["rho"], liq["cp"], liq["w"]
  //   968.09999 5.997563 1134.5855
  //
  //   References
  //   ----------
  //   IAPWS, Guideline on Thermodynamic Properties of Supercooled Water,
  //   http://iapws.org/relguide/Supercooled.html
  //

  public static Dictionary<string, double> _Supercooled(double T, double P)
  {
    // Check input in range of validity
    if(P < 198.9)
    {
      var Tita = T / 235.15;
      var Ph = 0.1 + 228.27 * (1 - Math.Pow(Tita, 6.243)) + 15.724 * (1 - Math.Pow(Tita, 79.81));
      if(P < Ph || T > 300)
      {
        throw new NotImplementedException("Incoming out of bound");
      }
    }
    else
    {
      var Th = 172.82 + 0.03718 * P + 3.403E-05 * Math.Pow(P, 2) - 1.573E-08 * Math.Pow(P, 3);
      if(T < Th || T > 300 || P > 1000)
      {
        throw new NotImplementedException("Incoming out of bound");
      }
    }
    // Parameters, Table 1
    var Tll = 228.2;
    var rho0 = 1081.6482;
    var R = 0.461523087;
    var pi0 = 300000.0 / rho0 / R / Tll;
    var omega0 = 0.5212269;
    var L0 = 0.76317954;
    var k0 = 0.072158686;
    var k1 = -0.31569232;
    var k2 = 5.2992608;
    // Reducing parameters, Eq 2
    var tau = T / Tll - 1;
    var p = P * 1000 / rho0 / R / Tll;
    var tau_ = tau + 1;
    var p_ = p + pi0;
    // Eq 3
    var ci = new List<double>
    {
      -8.1570681381655,
      1.2875032,
      7.0901673598012,
      -0.032779161,
      0.73703949,
      -0.21628622,
      -5.1782479,
      0.00042293517,
      0.023592109,
      4.3773754,
      -0.002996777,
      -0.96558018,
      3.7595286,
      1.2632441,
      0.28542697,
      -0.85994947,
      -0.32916153,
      0.090019616,
      0.081149726,
      -3.2788213
    };
    var ai = new List<double>
    {
      0,
      0,
      1,
      -0.2555,
      1.5762,
      1.64,
      3.6385,
      -0.3828,
      1.6219,
      4.3287,
      3.4763,
      5.1556,
      -0.3593,
      5.0361,
      2.9786,
      6.2373,
      4.046,
      5.3558,
      9.0157,
      1.2194
    };
    var bi = new List<double>
    {
      0,
      1,
      0,
      2.1051,
      1.1422,
      0.951,
      0,
      3.6402,
      2.076,
      -0.0016,
      2.2769,
      0.0008,
      0.3706,
      -0.3975,
      2.973,
      -0.318,
      2.9805,
      2.9265,
      0.4456,
      0.1298
    };
    var di = new List<double>
    {
      0,
      0,
      0,
      -0.0016,
      0.6894,
      0.013,
      0.0002,
      0.0435,
      0.05,
      0.0004,
      0.0528,
      0.0147,
      0.8584,
      0.9924,
      1.0041,
      1.0961,
      1.0228,
      1.0303,
      1.618,
      0.5213
    };
    var phir = 0;
    foreach (var _tup_1 in zip(ci, ai, bi, di))
    {
      var c = _tup_1.Item1;
      var a = _tup_1.Item2;
      var b = _tup_1.Item3;
      var d = _tup_1.Item4;
      phir += c * Math.Pow(tau_, a) * Math.Pow(p_, b) * Math.Exp(-d * p_);
      phirt += c * a * Math.Pow(tau_, a - 1) * Math.Pow(p_, b) * Math.Exp(-d * p_);
      phirp += c * Math.Pow(tau_, a) * Math.Pow(p_, b - 1) * (b - d * p_) * Math.Exp(-d * p_);
      phirtt += c * a * (a - 1) * Math.Pow(tau_, a - 2) * Math.Pow(p_, b) * Math.Exp(-d * p_);
      phirtp += c * a * Math.Pow(tau_, a - 1) * Math.Pow(p_, b - 1) * (b - d * p_) * Math.Exp(-d * p_);
      phirpp += c * Math.Pow(tau_, a) * Math.Pow(p_, b - 2) * (Math.Pow(d * p_ - b, 2) - b) * Math.Exp(-d * p_);
    }
    // Eq 5
    var K1 = Math.Pow(Math.Pow(1 + k0 * k2 + k1 * (p - k2 * tau), 2) - 4 * k0 * k1 * k2 * (p - k2 * tau), 0.5);
    var K2 = Math.Pow(1 + Math.Pow(k2, 2), 0.5);
    // Eq 6
    var omega = 2 + omega0 * p;
    // Eq 4
    var L = L0 * K2 / 2 / k1 / k2 * (1 + k0 * k2 + k1 * (p + k2 * tau) - K1);
    // Define interval of solution, Table 4
    if(omega < 10 / 9 * (Math.Log(19) - L))
    {
      var xmin = 0.049;
      var xmax = 0.5;
    }
    else if(10 / 9 * (Math.Log(19) - L) <= omega < 50 / 49 * (Math.Log(99) - L))
    {
      xmin = 0.0099;
      xmax = 0.051;
    }
    else
    {
      xmin = 0.99 * Math.Exp(-50 / 49 * L - omega);
      xmax = min(1.1 * Math.Exp(-L - omega), 0.0101);
    }
    Func<double, double> f = x =>
    {
      return abs(L + Math.Log(x / (1 - x)) + omega * (1 - 2 * x));
    };
    var x = minimize(f, ValueTuple.Create((xmin + xmax) / 2), bounds: ValueTuple.Create((xmin, xmax)))["x"][0];
    // Eq 12
    var fi = 2 * x - 1;
    var Xi = 1 / (2 / (1 - Math.Pow(fi, 2)) - omega);
    // Derivatives, Table 3
    var Lt = L0 * K2 / 2 * (1 + (1 - k0 * k2 + k1 * (p - k2 * tau)) / K1);
    var Lp = L0 * K2 * (K1 + k0 * k2 - k1 * p + k1 * k2 * tau - 1) / 2 / k2 / K1;
    var Ltt = -2 * L0 * K2 * k0 * k1 * Math.Pow(k2, 2) / Math.Pow(K1, 3);
    var Ltp = 2 * L0 * K2 * k0 * k1 * k2 / Math.Pow(K1, 3);
    var Lpp = -2 * L0 * K2 * k0 * k1 / Math.Pow(K1, 3);

    var prop = new Dictionary<string, double> { };
    prop["L"] = L;
    prop["x"] = x;
    // Eq 13
    prop["rho"] = rho0 / ((tau + 1) / 2 * (omega0 / 2 * (1 - Math.Pow(fi, 2)) + Lp * (fi + 1)) + phirp);
    // Eq 1
    prop["g"] = phir + (tau + 1) * (x * L + x * Math.Log(x) + (1 - x) * Math.Log(1 - x) + omega * x * (1 - x));
    // Eq 14
    prop["s"] = -R * ((tau + 1) / 2 * Lt * (fi + 1) + (x * L + x * Math.Log(x) + (1 - x) * Math.Log(1 - x) + omega * x * (1 - x)) + phirt);
    // Basic derived state properties
    prop["h"] = prop["g"] + T * prop["s"];
    prop["u"] = prop["h"] + P / prop["rho"];
    prop["a"] = prop["u"] - T * prop["s"];
    // Eq 15
    prop["xkappa"] = prop["rho"] / Math.Pow(rho0, 2) / R * 1000 / Tll * ((tau + 1) / 2 * (Xi * Math.Pow(Lp - omega0 * fi, 2) - (fi + 1) * Lpp) - phirpp);
    prop["alfap"] = prop["rho"] / rho0 / Tll * (Ltp / 2 * (tau + 1) * (fi + 1) + (omega0 * (1 - Math.Pow(fi, 2)) / 2 + Lp * (fi + 1)) / 2 - (tau + 1) * Lt / 2 * Xi * (Lp - omega0 * fi) + phirtp);
    prop["cp"] = -R * (tau + 1) * (Lt * (fi + 1) + (tau + 1) / 2 * (Ltt * (fi + 1) - Math.Pow(Lt, 2) * Xi) + phirtt);
    // Eq 16
    prop["cv"] = prop["cp"] - T * Math.Pow(prop["alfap"], 2) / prop["rho"] / prop["xkappa"] * 1000.0;
    // Eq 17
    prop["w"] = Math.Pow(prop["rho"] * prop["xkappa"] * 1E-06 * prop["cv"] / prop["cp"], -0.5);

    return prop;
  }

  // Sublimation Pressure correlation
  //
  //   Parameters
  //   ----------
  //   T : float
  //     Temperature, [K]
  //
  //   Returns
  //   -------
  //   P : float
  //     Pressure at sublimation line, [MPa]
  //
  //   Notes
  //   ------
  //   Raise :class:`NotImplementedError` if input isn't in limit:
  //
  //     * 50 ≤ T ≤ 273.16
  //
  //   Examples
  //   --------
  //   >>> _Sublimation_Pressure(230)
  //   8.947352740189152e-06
  //
  //   References
  //   ----------
  //   IAPWS, Revised Release on the Pressure along the Melting and Sublimation
  //   Curves of Ordinary Water Substance, http://iapws.org/relguide/MeltSub.html.
  //

  public static double _Sublimation_Pressure(double T)
  {
    if(50 <= T <= 273.16)
    {
      var Tita = T / Tt;
      var suma = 0;
      var a = new List<double> { -21.2144006, 27.3203819, -6.1059813 };
      var expo = new List<double> { 0.00333333333, 1.20666667, 1.70333333 };
      foreach (var _tup_1 in zip(a, expo))
      {
        var ai = _tup_1.Item1;
        var expi = _tup_1.Item2;
        suma += ai * Math.Pow(Tita, expi);
      }
      return Math.Exp(suma / Tita) * Pt;
    }
    else
    {
      throw new NotImplementedException("Incoming out of bound");
    }
  }

  // Melting Pressure correlation
  //
  //   Parameters
  //   ----------
  //   T : float
  //     Temperature, [K]
  //   ice: string
  //     Type of ice: Ih, III, V, VI, VII.
  //     Below 273.15 is a mandatory input, the ice Ih is the default value.
  //     Above 273.15, the ice type is unnecesary.
  //
  //   Returns
  //   -------
  //   P : float
  //     Pressure at sublimation line, [MPa]
  //
  //   Notes
  //   ------
  //   Raise :class:`NotImplementedError` if input isn't in limit:
  //
  //     * 251.165 ≤ T ≤ 715
  //
  //   Examples
  //   --------
  //   >>> _Melting_Pressure(260)
  //   8.947352740189152e-06
  //   >>> _Melting_Pressure(254, "III")
  //   268.6846466336108
  //
  //   References
  //   ----------
  //   IAPWS, Revised Release on the Pressure along the Melting and Sublimation
  //   Curves of Ordinary Water Substance, http://iapws.org/relguide/MeltSub.html.
  //

  public static double _Melting_Pressure(double T, string ice = "Ih")
  {
    if(ice == "Ih" && 251.165 <= T <= 273.16)
    {
      // Ice Ih
      var Tref = Tt;
      var Pref = Pt;
      var Tita = T / Tref;
      var a = new List<double>
      {
        1195393.37,
        80818.3159,
        3338.2686
      };
      var expo = new List<double>
      {
        3.0,
        25.75,
        103.75
      };
      var suma = 1;
      foreach (var _tup_1 in zip(a, expo))
      {
        var ai = _tup_1.Item1;
        var expi = _tup_1.Item2;
        suma += ai * (1 - Math.Pow(Tita, expi));
      }
      var P = suma * Pref;
    }
    else if(ice == "III" && 251.165 < T <= 256.164)
    {
      // Ice III
      Tref = 251.165;
      Pref = 208.566;
      Tita = T / Tref;
      P = Pref * (1 - 0.299948 * (1 - Math.Pow(Tita, 60.0)));
    }
    else if(ice == "V" && 256.164 < T <= 273.15 || 273.15 < T <= 273.31)
    {
      // Ice V
      Tref = 256.164;
      Pref = 350.1;
      Tita = T / Tref;
      P = Pref * (1 - 1.18721 * (1 - Math.Pow(Tita, 8.0)));
    }
    else if(273.31 < T <= 355)
    {
      // Ice VI
      Tref = 273.31;
      Pref = 632.4;
      Tita = T / Tref;
      P = Pref * (1 - 1.07476 * (1 - Math.Pow(Tita, 4.6)));
    }
    else if(355.0 < T <= 715)
    {
      // Ice VII
      Tref = 355;
      Pref = 2216.0;
      Tita = T / Tref;
      P = Pref * Math.Exp(1.73683 * (1 - 1.0 / Tita) - 0.0544606 * (1 - Math.Pow(Tita, 5)) + 8.06106E-08 * (1 - Math.Pow(Tita, 22)));
    }
    else
    {
      throw new NotImplementedException("Incoming out of bound");
    }
    return P;
  }

  // Transport properties
  // Equation for the Viscosity
  //
  //   Parameters
  //   ----------
  //   rho : float
  //     Density, [kg/m³]
  //   T : float
  //     Temperature, [K]
  //   fase: dict, optional for calculate critical enhancement
  //     phase properties
  //   drho: float, optional for calculate critical enhancement
  //     [∂ρ/∂P]T at reference state,
  //
  //   Returns
  //   -------
  //   μ : float
  //     Viscosity, [Pa·s]
  //
  //   Examples
  //   --------
  //   >>> _Viscosity(998, 298.15)
  //   0.0008897351001498108
  //   >>> _Viscosity(600, 873.15)
  //   7.743019522728247e-05
  //
  //   References
  //   ----------
  //   IAPWS, Release on the IAPWS Formulation 2008 for the Viscosity of Ordinary
  //   Water Substance, http://www.iapws.org/relguide/viscosity.html
  //

  public static double _Viscosity(double rho, double T, _utils._fase fase = null, double drho = 0.0)
  {
    object mu2;
    object Lw;
    object Y;
    var Tr = T / Tc;
    var Dr = rho / rhoc;
    // Eq 11
    var H = new List<double> { 1.67752, 2.20462, 0.6366564, -0.241605 };
    var mu0 = 100 * Math.Pow(Tr, 0.5) / (from _tup_1 in H.Select((_p_1,_p_2) => Tuple.Create(_p_2, _p_1)).Chop((i,Hi) => (i, Hi))
      let i = _tup_1.Item1
      let Hi = _tup_1.Item2
      select (Hi / Math.Pow(Tr, i))).ToList().Sum();
    // Eq 12
    var I = new List<int>
    {
      0,
      1,
      2,
      3,
      0,
      1,
      2,
      3,
      5,
      0,
      1,
      2,
      3,
      4,
      0,
      1,
      0,
      3,
      4,
      3,
      5
    };
    var J = new List<int>
    {
      0,
      0,
      0,
      0,
      1,
      1,
      1,
      1,
      1,
      2,
      2,
      2,
      2,
      2,
      3,
      3,
      4,
      4,
      5,
      6,
      6
    };
    var Hij = new List<double>
    {
      0.520094,
      0.0850895,
      -1.08374,
      -0.289555,
      0.222531,
      0.999115,
      1.88797,
      1.26613,
      0.120573,
      -0.281378,
      -0.906851,
      -0.772479,
      -0.489837,
      -0.25704,
      0.161913,
      0.257399,
      -0.0325372,
      0.0698452,
      0.00872102,
      -0.00435673,
      -0.000593264
    };
    var mu1 = Math.Exp(Dr * (from _tup_2 in zip(I, J, Hij).Chop((i,j,h) => (i, j, h))
      let i = _tup_2.Item1
      let j = _tup_2.Item2
      let h = _tup_2.Item3
      select (Math.Pow(1 / Tr - 1, i) * h * Math.Pow(Dr - 1, j))).ToList().Sum());
    // Critical enhancement
    if(fase && drho)
    {
      var qc = 1 / 1.9;
      var qd = 1 / 1.1;
      // Eq 21
      var DeltaX = Pc * Math.Pow(Dr, 2) * (fase.drhodP_T / rho - drho / rho * 1.5 / Tr);
      if(DeltaX < 0)
      {
        DeltaX = 0;
      }
      // Eq 20
      var X = 0.13 * Math.Pow(DeltaX / 0.06, 0.63 / 1.239);
      if(X <= 0.3817016416)
      {
        // Eq 15
        Y = qc / 5 * X * Math.Pow(qd * X, 5) * (1 - qc * X + Math.Pow(qc * X, 2) - 765.0 / 504 * Math.Pow(qd * X, 2));
      }
      else
      {
        var Fid = Math.Acos(Math.Pow(1 + Math.Pow(qd, 2) * Math.Pow(X, 2), -0.5));
        var w = Math.Pow(abs((qc * X - 1) / (qc * X + 1)), 0.5) * Matn.Tan(Fid / 2);
        // Eq 18
        if(qc * X > 1)
        {
          Lw = Math.Log((1 + w) / (1 - w));
        }
        else
        {
          Lw = 2 * Math.Atan(abs(w));
        }
        // Eq 16
        Y = Math.Sin(3 * Fid) / 12 - Math.Sin(2 * Fid) / 4 / qc / X + (1 - 5 / 4 * Math.Pow(qc * X, 2)) / Math.Pow(qc * X, 2) * Math.Sin(Fid) - ((1 - 3 / 2 * Math.Pow(qc * X, 2)) * Fid - Math.Pow(abs(Math.Pow(qc * X, 2) - 1), 1.5) * Lw) / Math.Pow(qc * X, 3);
      }
      // Eq 14
      mu2 = Math.Exp(0.068 * Y);
    }
    else
    {
      mu2 = 1;
    }
    // Eq 10
    var mu = mu0 * mu1 * mu2;
    return mu * 1E-06;
  }

  // Equation for the thermal conductivity
  //
  //   Parameters
  //   ----------
  //   rho : float
  //     Density, [kg/m³]
  //   T : float
  //     Temperature, [K]
  //   fase: dict, optional for calculate critical enhancement
  //     phase properties
  //   drho: float, optional for calculate critical enhancement
  //     [∂ρ/∂P]T at reference state,
  //
  //   Returns
  //   -------
  //   k : float
  //     Thermal conductivity, [W/mK]
  //
  //   Examples
  //   --------
  //   >>> _ThCond(998, 298.15)
  //   0.6077128675880629
  //   >>> _ThCond(0, 873.15)
  //   0.07910346589648833
  //
  //   References
  //   ----------
  //   IAPWS, Release on the IAPWS Formulation 2011 for the Thermal Conductivity
  //   of Ordinary Water Substance, http://www.iapws.org/relguide/ThCond.html
  //

  public static double _ThCond(double rho, double T, _utils._fase fase = null, double drho = 0.0)
  {
    object k2;
    object Z;
    var d = rho / rhoc;
    var Tr = T / Tc;
    // Eq 16
    var no = new List<double> { 0.002443221, 0.01323095, 0.006770357, -0.003454586, 0.0004096266 };
    var k0 = Math.Pow(Tr, 0.5) / (from _tup_1 in no.Select((_p_1,_p_2) => Tuple.Create(_p_2, _p_1)).Chop((i,n) => (i, n))
      let i = _tup_1.Item1
      let n = _tup_1.Item2
      select (n / Math.Pow(Tr, i))).ToList().Sum();
    // Eq 17
    var I = new List<int>
    {
      0,
      0,
      0,
      0,
      0,
      0,
      1,
      1,
      1,
      1,
      1,
      1,
      2,
      2,
      2,
      2,
      2,
      2,
      3,
      3,
      3,
      3,
      4,
      4,
      4,
      4,
      4,
      4
    };
    var J = new List<int>
    {
      0,
      1,
      2,
      3,
      4,
      5,
      0,
      1,
      2,
      3,
      4,
      5,
      0,
      1,
      2,
      3,
      4,
      5,
      0,
      1,
      2,
      3,
      0,
      1,
      2,
      3,
      4,
      5
    };
    var nij = new List<double>
    {
      1.60397357,
      -0.646013523,
      0.111443906,
      0.102997357,
      -0.0504123634,
      0.00609859258,
      2.33771842,
      -2.78843778,
      1.53616167,
      -0.463045512,
      0.0832827019,
      -0.00719201245,
      2.19650529,
      -4.54580785,
      3.55777244,
      -1.40944978,
      0.275418278,
      -0.0205938816,
      -1.21051378,
      1.60812989,
      -0.621178141,
      0.0716373224,
      -2.720337,
      4.57586331,
      -3.18369245,
      1.1168348,
      -0.19268305,
      0.012913842
    };
    var k1 = Math.Exp(d * (from _tup_2 in zip(I, J, nij).Chop((i,j,n) => (i, j, n))
      let i = _tup_2.Item1
      let j = _tup_2.Item2
      let n = _tup_2.Item3
      select (Math.Pow(1 / Tr - 1, i) * n * Math.Pow(d - 1, j))).ToList().Sum());
    // Critical enhancement
    if(fase)
    {
      var R = 0.46151805;
      if(drho == 0.0)
      {
        // Industrial formulation
        // Eq 25
        if(d <= 0.310559006)
        {
          var ai = new List<double> { 6.53786807199516, -5.61149954923348, 3.39624167361325, -2.27492629730878, 10.2631854662709, 1.97815050331519 };
        }
        else if(d <= 0.776397516)
        {
          ai = new List<double> { 6.52717759281799, -6.30816983387575, 8.08379285492595, -9.82240510197603, 12.1358413791395, -5.54349664571295 };
        }
        else if(d <= 1.242236025)
        {
          ai = new List<double> { 5.35500529896124, -3.96415689925446, 8.91990208918795, -12.033872950579, 9.19494865194302, -2.16866274479712 };
        }
        else if(d <= 1.863354037)
        {
          ai = new List<double> { 1.55225959906681, 0.464621290821181, 8.93237374861479, -11.0321960061126, 6.1678099993336, -0.965458722086812 };
        }
        else
        {
          ai = new List<double> { 1.11999926419994, 0.595748562571649, 9.8895256507892, -10.325505114704, 4.66861294457414, -0.503243546373828 };
        }
        drho = 1 / (from _tup_3 in ai.Select((_p_3,_p_4) => Tuple.Create(_p_4, _p_3)).Chop((i,a) => (i, a))
          let i = _tup_3.Item1
          let a = _tup_3.Item2
          select (a * Math.Pow(d, i))).ToList().Sum() * rhoc / Pc;
      }
      var DeltaX = d * (Pc / rhoc * fase.drhodP_T - Pc / rhoc * drho * 1.5 / Tr);
      if(DeltaX < 0)
      {
        DeltaX = 0;
      }
      var X = 0.13 * Math.Pow(DeltaX / 0.06, 0.63 / 1.239);
      var y = X / 0.4;
      // Eq 19
      if(y < 1.2E-07)
      {
        Z = 0;
      }
      else
      {
        Z = 2 / Math.PI / y * ((1 - 1 / fase.cp_cv) * Math.Atan(y) + y / fase.cp_cv - (1 - Math.Exp(-1 / (1 / y + Math.Pow(y, 2) / 3 / Math.Pow(d, 2)))));
      }
      // Eq 18
      k2 = 177.8514 * d * fase.cp / R * Tr / fase.mu * 1E-06 * Z;
    }
    else
    {
      // No critical enhancement
      k2 = 0;
    }
    // Eq 10
    var k = k0 * k1 + k2;
    return 0.001 * k;
  }

  // Equation for the surface tension
  //
  //   Parameters
  //   ----------
  //   T : float
  //     Temperature, [K]
  //
  //   Returns
  //   -------
  //   σ : float
  //     Surface tension, [N/m]
  //
  //   Notes
  //   ------
  //   Raise :class:`NotImplementedError` if input isn't in limit:
  //
  //     * 248.15 ≤ T ≤ 647
  //     * Estrapolate to -25ºC in supercooled liquid metastable state
  //
  //   Examples
  //   --------
  //   >>> _Tension(300)
  //   0.0716859625
  //   >>> _Tension(450)
  //   0.0428914992
  //
  //   References
  //   ----------
  //   IAPWS, Revised Release on Surface Tension of Ordinary Water Substance
  //   June 2014, http://www.iapws.org/relguide/Surf-H2O.html
  //

  public static double _Tension(double T)
  {
    if(248.15 <= T <= Tc)
    {
      var Tr = T / Tc;
      return 0.001 * (235.8 * Math.Pow(1 - Tr, 1.256) * (1 - 0.625 * (1 - Tr)));
    }
    else
    {
      throw new NotImplementedException("Incoming out of bound");
    }
  }

  // Equation for the Dielectric constant
  //
  //   Parameters
  //   ----------
  //   rho : float
  //     Density, [kg/m³]
  //   T : float
  //     Temperature, [K]
  //
  //   Returns
  //   -------
  //   epsilon : float
  //     Dielectric constant, [-]
  //
  //   Notes
  //   ------
  //   Raise :class:`NotImplementedError` if input isn't in limit:
  //
  //     * 238 ≤ T ≤ 1200
  //
  //   Examples
  //   --------
  //   >>> _Dielectric(999.242866, 298.15)
  //   78.5907250
  //   >>> _Dielectric(26.0569558, 873.15)
  //   1.12620970
  //
  //   References
  //   ----------
  //   IAPWS, Release on the Static Dielectric Constant of Ordinary Water
  //   Substance for Temperatures from 238 K to 873 K and Pressures up to 1000
  //   MPa, http://www.iapws.org/relguide/Dielec.html
  //

  public static double _Dielectric(double rho, double T)
  {
    // Check input parameters
    if(T < 238 || T > 1200)
    {
      throw new NotImplementedException("Incoming out of bound");
    }
    var k = 1.380658E-23;
    var Na = 6.0221367E+23;
    var alfa = 1.636E-40;
    var epsilon0 = 8.854187817E-12;
    var mu = 6.138E-30;
    var d = rho / rhoc;
    var Tr = Tc / T;
    var I = new List<int>
    {
      1,
      1,
      1,
      2,
      3,
      3,
      4,
      5,
      6,
      7,
      10,
      null
    };
    var J = new List<double>
    {
      0.25,
      1,
      2.5,
      1.5,
      1.5,
      2.5,
      2,
      2,
      5,
      0.5,
      10,
      null
    };
    var n = new List<double>
    {
      0.978224486826,
      -0.957771379375,
      0.237511794148,
      0.714692244396,
      -0.298217036956,
      -0.108863472196,
      0.0949327488264,
      -0.00980469816509,
      1.6516763497E-05,
      9.37359795772E-05,
      -1.2317921872E-10,
      0.00196096504426
    };
    var g = 1 + n[11] * d / Math.Pow(Tc / 228 / Tr - 1, 1.2);
    foreach (var i in Enumerable.Range(0, 11))
    {
      g += n[i] * Math.Pow(d, I[i]) * Math.Pow(Tr, J[i]);
    }
    var A = Na * Math.Pow(mu, 2) * rho * g / M * 1000 / epsilon0 / k / T;
    var B = Na * alfa * rho / 3 / M * 1000 / epsilon0;
    var e = (1 + A + 5 * B + Math.Pow(9 + 2 * A + 18 * B + Math.Pow(A, 2) + 10 * A * B + 9 * Math.Pow(B, 2), 0.5)) / 4 / (1 - B);
    return e;
  }

  // Equation for the refractive index
  //
  //   Parameters
  //   ----------
  //   rho : float
  //     Density, [kg/m³]
  //   T : float
  //     Temperature, [K]
  //   l : float, optional
  //     Light Wavelength, [μm]
  //
  //   Returns
  //   -------
  //   n : float
  //     Refractive index, [-]
  //
  //   Notes
  //   ------
  //   Raise :class:`NotImplementedError` if input isn't in limit:
  //
  //     * 0 ≤ ρ ≤ 1060
  //     * 261.15 ≤ T ≤ 773.15
  //     * 0.2 ≤ λ ≤ 1.1
  //
  //   Examples
  //   --------
  //   >>> _Refractive(997.047435, 298.15, 0.2265)
  //   1.39277824
  //   >>> _Refractive(30.4758534, 773.15, 0.5893)
  //   1.00949307
  //
  //   References
  //   ----------
  //   IAPWS, Release on the Refractive Index of Ordinary Water Substance as a
  //   Function of Wavelength, Temperature and Pressure,
  //   http://www.iapws.org/relguide/rindex.pdf
  //

  public static double _Refractive(double rho, double T, double l = 0.5893)
  {
    // Check input parameters
    if(rho < 0 || rho > 1060 || T < 261.15 || T > 773.15 || l < 0.2 || l > 1.1)
    {
      throw new NotImplementedException("Incoming out of bound");
    }
    var Lir = 5.432937;
    var Luv = 0.229202;
    var d = rho / 1000.0;
    var Tr = T / 273.15;
    var L = l / 0.589;
    var a = new List<double> { 0.244257733, 0.00974634476, -0.00373234996, 0.000268678472, 0.0015892057, 0.00245934259, 0.90070492, -0.0166626219 };
    var A = d * (a[0] + a[1] * d + a[2] * Tr + a[3] * Math.Pow(L, 2) * Tr + a[4] / Math.Pow(L, 2) + a[5] / (Math.Pow(L, 2) - Math.Pow(Luv, 2)) + a[6] / (Math.Pow(L, 2) - Math.Pow(Lir, 2)) + a[7] * Math.Pow(d, 2));
    return Math.Pow((2 * A + 1) / (1 - A), 0.5);
  }

  // Equation for the ionization constant of ordinary water
  //
  //   Parameters
  //   ----------
  //   rho : float
  //     Density, [kg/m³]
  //   T : float
  //     Temperature, [K]
  //
  //   Returns
  //   -------
  //   pKw : float
  //     Ionization constant in -log10(kw), [-]
  //
  //   Notes
  //   ------
  //   Raise :class:`NotImplementedError` if input isn't in limit:
  //
  //     * 0 ≤ ρ ≤ 1250
  //     * 273.15 ≤ T ≤ 1073.15
  //
  //   Examples
  //   --------
  //   >>> _Kw(1000, 300)
  //   13.906565
  //
  //   References
  //   ----------
  //   IAPWS, Release on the Ionization Constant of H2O,
  //   http://www.iapws.org/relguide/Ionization.pdf
  //

  public static double _Kw(double rho, double T)
  {
    // Check input parameters
    if(rho < 0 || rho > 1250 || T < 273.15 || T > 1073.15)
    {
      throw new NotImplementedException("Incoming out of bound");
    }
    // The internal method of calculation use rho in g/cm³
    var d = rho / 1000.0;
    // Water molecular weight different
    var Mw = 18.015268;
    var gamma = new List<double> { 0.61415, 48251.33, -67707.93, 10102100.0 };
    var pKg = 0;
    foreach (var _tup_1 in gamma.Select((_p_1,_p_2) => Tuple.Create(_p_2, _p_1)))
    {
      var i = _tup_1.Item1;
      var g = _tup_1.Item2;
      pKg += g / Math.Pow(T, i);
    }
    var Q = d * Math.Exp(-0.864671 + 8659.19 / T - 22786.2 / Math.Pow(T, 2) * Math.Pow(d, 2.0 / 3));
    var pKw = -12 * (Math.Log(1 + Q,10) - Q / (Q + 1) * d * (0.642044 - 56.8534 / T - 0.375754 * d)) + pKg + 2 * Math.Log(Mw / 1000,10);
    return pKw;
  }

  // Equation for the electrolytic conductivity of liquid and dense
  //   supercrítical water
  //
  //   Parameters
  //   ----------
  //   rho : float
  //     Density, [kg/m³]
  //   T : float
  //     Temperature, [K]
  //
  //   Returns
  //   -------
  //   K : float
  //     Electrolytic conductivity, [S/m]
  //
  //   Notes
  //   ------
  //   Raise :class:`NotImplementedError` if input isn't in limit:
  //
  //     * 600 ≤ ρ ≤ 1200
  //     * 273.15 ≤ T ≤ 1073.15
  //
  //   Examples
  //   --------
  //   >>> _Conductivity(1000, 373.15)
  //   1.13
  //
  //   References
  //   ----------
  //   IAPWS, Electrolytic Conductivity (Specific Conductance) of Liquid and Dense
  //   Supercritical Water from 0°C to 800°C and Pressures up to 1000 MPa,
  //   http://www.iapws.org/relguide/conduct.pdf
  //

  public static double _Conductivity(double rho, double T)
  {
    // FIXME: Dont work
    var rho_ = rho / 1000;
    var kw = Math.Pow(10, -_Kw(rho, T));
    var A = new List<double> { 1850.0, 1410.0, 2.16417E-06, 1.81609E-07, -1.75297E-09, 7.20708E-12 };
    var B = new List<double> { 16.0, 11.6, 0.000326, -2.3E-06, 1.1E-08 };
    var t = T - 273.15;
    var Loo = A[0] - 1 / (1 / A[1] + (from i in Enumerable.Range(0, 4)
      select (A[i + 2] * Math.Pow(t, i + 1))).ToList().Sum());
    var rho_h = B[0] - 1 / (1 / B[1] + (from i in Enumerable.Range(0, 3)
      select (B[i + 2] * Math.Pow(t, i + 1))).ToList().Sum());
    // Eq 4
    var L_o = (rho_h - rho_) * Loo / rho_h;
    // Eq 1
    var k = 100 * 0.001 * L_o * Math.Pow(kw, 0.5) * rho_;
    return k;
  }

  // Equation for the calculation of Henry's constant
  //
  //   Parameters
  //   ----------
  //   T : float
  //     Temperature, [K]
  //   gas : string
  //     Name of gas to calculate solubility
  //   liquid : string
  //     Name of liquid solvent, can be H20 (default)
  //
  //   Returns
  //   -------
  //   kw : float
  //     Henry's constant, [MPa]
  //
  //   Notes
  //   -----
  //   The gas availables for H2O solvent are He, Ne, Ar, Kr, Xe, H2, N2, O2, CO,
  //   CO2, H2S, CH4, C2H6, SF6
  //
  //   Raise :class:`NotImplementedError` if input gas or liquid are unsupported
  //
  //   Examples
  //   --------
  //   >>> _Henry(500, "He")
  //   1.1973
  //
  //   References
  //   ----------
  //   IAPWS, Guideline on the Henry's Constant and Vapor-Liquid Distribution
  //   Constant for Gases in H2O at High Temperatures,
  //   http://www.iapws.org/relguide/HenGuide.html
  //

  public static double _Henry(double T, string gas, string liquid = "H2O")
  {
    object bi;
    object ai;
    object Pc;
    object Tc;
    var limit = new Dictionary<string, object>
    {
      { "He", (273.21, 553.18)},
      { "Ne", (273.2, 543.36)},
      { "Ar", (273.19, 568.36)},
      { "Kr", (273.19, 525.56)},
      { "Xe", (273.22, 574.85)},
      { "H2", (273.15, 636.09)},
      { "N2", (278.12, 636.46)},
      { "O2", (274.15, 616.52)},
      { "CO", (278.15, 588.67)},
      { "CO2", (274.19, 642.66)},
      { "H2S", (273.15, 533.09)},
      { "CH4", (275.46, 633.11)},
      { "C2H6", (275.44, 473.46)},
      { "SF6", (283.14, 505.55)},
    };
    // Check input parameters
    if(liquid != "H2O")
    {
      throw new NotImplementedException("Solvent liquid unsupported");
    }
    if(!limit.Contains(gas))
    {
      throw new NotImplementedException("Gas unsupported");
    }
    var _tup_1 = limit[gas];
    var Tmin = _tup_1.Item1;
    var Tmax = _tup_1.Item2;
    if(T < Tmin || T > Tmax)
    {
      Debug.Log("Temperature out of data of correlation");
    }

    Tc = 647.096;
    Pc = 22.064;
    var Tr = T / Tc;
    var tau = 1 - Tr;
    // Eq 4
    if(liquid == "H2O")
    {
      ai = new List<double> { -7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502 };
      bi = new List<double> { 1, 1.5, 3, 3.5, 4, 7.5 };
    }
    else
    {
      ai = new List<double> { -7.896657, 24.73308, -27.81128, 9.355913, -9.220083 };
      bi = new List<double> { 1, 1.89, 2, 3, 3.6 };
    }
    var ps = Pc * Math.Exp(1 / Tr * (from _tup_2 in zip(ai, bi).Chop((a,b) => (a, b))
      let a = _tup_2.Item1
      let b = _tup_2.Item2
      select (a * Math.Pow(tau, b))).ToList().Sum());
    // Select values from Table 2
    var par = new Dictionary<string, object>
    {
      { "He", (-3.52839, 7.12983, 4.4777)},
      { "Ne", (-3.18301, 5.31448, 5.43774)},
      { "Ar", (-8.40954, 4.29587, 10.52779)},
      { "Kr", (-8.97358, 3.61508, 11.29963)},
      { "Xe", (-14.21635, 4.00041, 15.60999)},
      { "H2", (-4.73284, 6.08954, 6.06066)},
      { "N2", (-9.67578, 4.72162, 11.70585)},
      { "O2", (-9.44833, 4.43822, 11.42005)},
      { "CO", (-10.52862, 5.13259, 12.01421)},
      { "CO2", (-8.55445, 4.01195, 9.52345)},
      { "H2S", (-4.51499, 5.23538, 4.42126)},
      { "CH4", (-10.44708, 4.66491, 12.12986)},
      { "C2H6", (-19.67563, 4.51222, 20.62567)},
      { "SF6", (-16.56118, 2.15289, 20.3544)},
    };
    var _tup_3 = par[gas];
    var A = _tup_3.Item1;
    var B = _tup_3.Item2;
    var C = _tup_3.Item3;
    // Eq 3
    var kh = ps * Math.Exp(A / Tr + B * Math.Pow(tau, 0.355) / Tr + C * Math.Pow(Tr, -0.41) * Math.Exp(tau));
    return kh;
  }

  // Equation for the vapor-liquid distribution constant
  //
  //   Parameters
  //   ----------
  //   T : float
  //     Temperature, [K]
  //   gas : string
  //     Name of gas to calculate solubility
  //   liquid : string
  //     Name of liquid solvent, can be H20 (default)
  //
  //   Returns
  //   -------
  //   kd : float
  //     Vapor-liquid distribution constant, [-]
  //
  //   Notes
  //   -----
  //   The gas availables for H2O solvent are He, Ne, Ar, Kr, Xe, H2, N2, O2, CO,
  //   CO2, H2S, CH4, C2H6, SF6
  //
  //   Raise :class:`NotImplementedError` if input gas or liquid are unsupported
  //
  //   Examples
  //   --------
  //   >>> _Kvalue(600, "He")
  //   3.8019
  //
  //   References
  //   ----------
  //   IAPWS, Guideline on the Henry's Constant and Vapor-Liquid Distribution
  //   Constant for Gases in H2O at High Temperatures,
  //   http://www.iapws.org/relguide/HenGuide.html
  //

  public static double _Kvalue(double T, string gas, string liquid = "H2O")
  {
    double q;
    List<int> di;
    List<double> ci;
    double Tc;
    var limit = new Dictionary<string, object>
    {
      { "He", (273.21, 553.18)},
      { "Ne", (273.2, 543.36)},
      { "Ar", (273.19, 568.36)},
      { "Kr", (273.19, 525.56)},
      { "Xe", (273.22, 574.85)},
      { "H2", (273.15, 636.09)},
      { "N2", (278.12, 636.46)},
      { "O2", (274.15, 616.52)},
      { "CO", (278.15, 588.67)},
      { "CO2", (274.19, 642.66)},
      { "H2S", (273.15, 533.09)},
      { "CH4", (275.46, 633.11)},
      { "C2H6", (275.44, 473.46)},
      { "SF6", (283.14, 505.55)},
    };
    // Check input parameters
    if(liquid != "H2O")
    {
      throw new NotImplementedException("Solvent liquid unsupported");
    }
    if(!limit.Contains(gas))
    {
      throw new NotImplementedException("Gas unsupported");
    }
    var _tup_1 = limit[gas];
    var Tmin = _tup_1.Item1;
    var Tmax = _tup_1.Item2;
    if(T < Tmin || T > Tmax)
    {
      Debug.Log("Temperature out of data of correlation");
    }

    Tc = 647.096;
    var Tr = T / Tc;
    var tau = 1 - Tr;
    // Eq 6
    if(liquid == "H2O")
    {
      ci = new List<double> { 1.99274064, 1.09965342, -0.510839303, -1.75493479, -45.5170352, -674694.45 };
      di = new List<int> { 1 / 3, 2 / 3, 5 / 3, 16 / 3, 43 / 3, 110 / 3 };
      q = -0.023767;
    }
    else
    {
      ci = new List<double> { 2.7072, 0.58662, -1.3069, -45.663 };
      di = new List<double> { 0.374, 1.45, 2.6, 12.3 };
      q = -0.024552;
    }
    var f = (from _tup_2 in zip(ci, di).Chop((c,d) => (c, d))
      let c = _tup_2.Item1
      let d = _tup_2.Item2
      select (c * Math.Pow(tau, d))).ToList().Sum();
    // Select values from Table 2
    var par = new Dictionary<string, object>
    {
      { "He", (2267.4082, -2.9616, -3.2604, 7.8819)},
      { "Ne", (2507.3022, -38.6955, 110.3992, -71.9096)},
      { "Ar", (2310.5463, -46.7034, 160.4066, -118.3043)},
      { "Kr", (2276.9722, -61.1494, 214.0117, -159.0407)},
      { "Xe", (2022.8375, 16.7913, -61.2401, 41.9236)},
      { "H2", (2286.4159, 11.3397, -70.7279, 63.0631)},
      { "N2", (2388.8777, -14.9593, 42.0179, -29.4396)},
      { "O2", (2305.0674, -11.324, 25.3224, -15.6449)},
      { "CO", (2346.2291, -57.6317, 204.5324, -152.6377)},
      { "CO2", (1672.9376, 28.1751, -112.4619, 85.3807)},
      { "H2S", (1319.1205, 14.1571, -46.8361, 33.2266)},
      { "CH4", (2215.6977, -0.1089, -6.624, 4.6789)},
      { "C2H6", (2143.8121, 6.8859, -12.6084, 0)},
      { "SF6", (2871.7265, -66.7556, 229.7191, -172.74)},
    };
    var _tup_3 = par[gas];
    var E = _tup_3.Item1;
    var F = _tup_3.Item2;
    var G = _tup_3.Item3;
    var H = _tup_3.Item4;
    // Eq 5
    var kd = Math.Exp(q * F + E / T * f + (F + G * Math.Pow(tau, 2.0 / 3) + H * tau) * Math.Exp((273.15 - T) / 100));
    return kd;
  }
}

