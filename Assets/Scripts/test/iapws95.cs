/*
using M = _iapws.M;
using Tc = _iapws.Tc;
using Pc = _iapws.Pc;
using rhoc = _iapws.rhoc;
*/
using System.Collections.Generic;
using System;
using System.Collections;
using System.Linq;

public static class iapws95
{

  static iapws95()
  {
  //Implemented multiparameter equation of state as a Helmholtz free energy:
  //
  //* :class:`MEoS`: Base class of multiparameter equation of state
  //* :class:`IAPWS95`: 2016 revision of 1995 formulation for ordinaty water
  }

  //!/usr/bin/python
  // -*- coding: utf-8 -*-
  // Residual contribution to the adimensional free Helmholtz energy
  //
  //   Parameters
  //   ----------
  //   tau : float
  //     Inverse reduced temperature Tc/T, [-]
  //   delta : float
  //     Reduced density rho/rhoc, [-]
  //   coef : dict
  //     Dictionary with equation of state parameters
  //
  //   Returns
  //   -------
  //   fir : float
  //     Adimensional free Helmholtz energy
  //
  //   References
  //   ----------
  //   IAPWS, Revised Release on the IAPWS Formulation 1995 for the
  //   Thermodynamic Properties of Ordinary Water Substance for General and
  //   Scientific Use, September 2016, Table 5
  //   http://www.iapws.org/relguide/IAPWS-95.html
  //

  public static object _phir(object tau, object delta, object coef)
  {
    object b;
    object a;
    object g;
    object t;
    object d;
    object n;
    var fir = 0;
    // Polinomial terms
    var nr1 = coef.get("nr1", new List<object>());
    var d1 = coef.get("d1", new List<object>());
    var t1 = coef.get("t1", new List<object>());
    foreach (var _tup_1 in zip(nr1, d1, t1))
    {
      n = _tup_1.Item1;
      d = _tup_1.Item2;
      t = _tup_1.Item3;
      fir += n * Math.Pow(delta, d) * Math.Pow(tau, t);
    }
    // Exponential terms
    var nr2 = coef.get("nr2", new List<object>());
    var d2 = coef.get("d2", new List<object>());
    var g2 = coef.get("gamma2", new List<object>());
    var t2 = coef.get("t2", new List<object>());
    var c2 = coef.get("c2", new List<object>());
    foreach (var _tup_2 in zip(nr2, d2, g2, t2, c2))
    {
      n = _tup_2.Item1;
      d = _tup_2.Item2;
      g = _tup_2.Item3;
      t = _tup_2.Item4;
      var c = _tup_2.Item5;
      fir += n * Math.Pow(delta, d) * Math.Pow(tau, t) * Math.Exp(-g * Math.Pow(delta, c));
    }
    // Gaussian terms
    var nr3 = coef.get("nr3", new List<object>());
    var d3 = coef.get("d3", new List<object>());
    var t3 = coef.get("t3", new List<object>());
    var a3 = coef.get("alfa3", new List<object>());
    var e3 = coef.get("epsilon3", new List<object>());
    var b3 = coef.get("beta3", new List<object>());
    var g3 = coef.get("gamma3", new List<object>());
    foreach (var _tup_3 in zip(nr3, d3, t3, a3, e3, b3, g3))
    {
      n = _tup_3.Item1;
      d = _tup_3.Item2;
      t = _tup_3.Item3;
      a = _tup_3.Item4;
      var e = _tup_3.Item5;
      b = _tup_3.Item6;
      g = _tup_3.Item7;
      fir += n * Math.Pow(delta, d) * Math.Pow(tau, t) * Math.Exp(-a * Math.Pow(delta - e, 2) - b * Math.Pow(tau - g, 2));
    }
    // Non analitic terms
    var nr4 = coef.get("nr4", new List<object>());
    var a4 = coef.get("a4", new List<object>());
    var b4 = coef.get("b4", new List<object>());
    var Ai = coef.get("A", new List<object>());
    var Bi = coef.get("B", new List<object>());
    var Ci = coef.get("C", new List<object>());
    var Di = coef.get("D", new List<object>());
    var bt4 = coef.get("beta4", new List<object>());
    foreach (var _tup_4 in zip(nr4, a4, b4, Ai, Bi, Ci, Di, bt4))
    {
      n = _tup_4.Item1;
      a = _tup_4.Item2;
      b = _tup_4.Item3;
      var A = _tup_4.Item4;
      var B = _tup_4.Item5;
      var C = _tup_4.Item6;
      var D = _tup_4.Item7;
      var bt = _tup_4.Item8;
      var Tita = 1 - tau + A * Math.Pow(Math.Pow(delta - 1, 2), 0.5 / bt);
      var F = Math.Exp(-C * Math.Pow(delta - 1, 2) - D * Math.Pow(tau - 1, 2));
      var Delta = Math.Pow(Tita, 2) + B * Math.Pow(Math.Pow(delta - 1, 2), a);
      fir += n * Math.Pow(Delta, b) * delta * F;
    }
    return fir;
  }

  // Residual contribution to the adimensional free Helmholtz energy, delta
  //   derivative
  //
  //   Parameters
  //   ----------
  //   tau : float
  //     Inverse reduced temperature Tc/T, [-]
  //   delta : float
  //     Reduced density rho/rhoc, [-]
  //   coef : dict
  //     Dictionary with equation of state parameters
  //
  //   Returns
  //   -------
  //   fird : float
  //     .. math::
  //       \left.\frac{\partial \phi^r_{\delta}}{\partial \delta}\right|_{\tau}
  //
  //   References
  //   ----------
  //   IAPWS, Revised Release on the IAPWS Formulation 1995 for the
  //   Thermodynamic Properties of Ordinary Water Substance for General and
  //   Scientific Use, September 2016, Table 5
  //   http://www.iapws.org/relguide/IAPWS-95.html
  //

  public static object _phird(object tau, object delta, object coef)
  {
    object b;
    object a;
    object g;
    object t;
    object d;
    object n;
    var fird = 0;
    // Polinomial terms
    var nr1 = coef.get("nr1", new List<object>());
    var d1 = coef.get("d1", new List<object>());
    var t1 = coef.get("t1", new List<object>());
    foreach (var _tup_1 in zip(nr1, d1, t1))
    {
      n = _tup_1.Item1;
      d = _tup_1.Item2;
      t = _tup_1.Item3;
      fird += n * d * Math.Pow(delta, d - 1) * Math.Pow(tau, t);
    }
    // Exponential terms
    var nr2 = coef.get("nr2", new List<object>());
    var d2 = coef.get("d2", new List<object>());
    var g2 = coef.get("gamma2", new List<object>());
    var t2 = coef.get("t2", new List<object>());
    var c2 = coef.get("c2", new List<object>());
    foreach (var _tup_2 in zip(nr2, d2, g2, t2, c2))
    {
      n = _tup_2.Item1;
      d = _tup_2.Item2;
      g = _tup_2.Item3;
      t = _tup_2.Item4;
      var c = _tup_2.Item5;
      fird += n * Math.Exp(-g * Math.Pow(delta, c)) * Math.Pow(delta, d - 1) * Math.Pow(tau, t) * (d - g * c * Math.Pow(delta, c));
    }
    // Gaussian terms
    var nr3 = coef.get("nr3", new List<object>());
    var d3 = coef.get("d3", new List<object>());
    var t3 = coef.get("t3", new List<object>());
    var a3 = coef.get("alfa3", new List<object>());
    var e3 = coef.get("epsilon3", new List<object>());
    var b3 = coef.get("beta3", new List<object>());
    var g3 = coef.get("gamma3", new List<object>());
    foreach (var _tup_3 in zip(nr3, d3, t3, a3, e3, b3, g3))
    {
      n = _tup_3.Item1;
      d = _tup_3.Item2;
      t = _tup_3.Item3;
      a = _tup_3.Item4;
      var e = _tup_3.Item5;
      b = _tup_3.Item6;
      g = _tup_3.Item7;
      fird += n * Math.Pow(delta, d) * Math.Pow(tau, t) * Math.Exp(-a * Math.Pow(delta - e, 2) - b * Math.Pow(tau - g, 2)) * (d / delta - 2 * a * (delta - e));
    }
    // Non analitic terms
    var nr4 = coef.get("nr4", new List<object>());
    var a4 = coef.get("a4", new List<object>());
    var b4 = coef.get("b4", new List<object>());
    var Ai = coef.get("A", new List<object>());
    var Bi = coef.get("B", new List<object>());
    var Ci = coef.get("C", new List<object>());
    var Di = coef.get("D", new List<object>());
    var bt4 = coef.get("beta4", new List<object>());
    foreach (var _tup_4 in zip(nr4, a4, b4, Ai, Bi, Ci, Di, bt4))
    {
      n = _tup_4.Item1;
      a = _tup_4.Item2;
      b = _tup_4.Item3;
      var A = _tup_4.Item4;
      var B = _tup_4.Item5;
      var C = _tup_4.Item6;
      var D = _tup_4.Item7;
      var bt = _tup_4.Item8;
      var Tita = 1 - tau + A * Math.Pow(Math.Pow(delta - 1, 2), 0.5 / bt);
      var F = Math.Exp(-C * Math.Pow(delta - 1, 2) - D * Math.Pow(tau - 1, 2));
      var Fd = -2 * C * F * (delta - 1);
      var Delta = Math.Pow(Tita, 2) + B * Math.Pow(Math.Pow(delta - 1, 2), a);
      var Deltad = (delta - 1) * (A * Tita * 2 / bt * Math.Pow(Math.Pow(delta - 1, 2), 0.5 / bt - 1) + 2 * B * a * Math.Pow(Math.Pow(delta - 1, 2), a - 1));
      var DeltaBd = b * Math.Pow(Delta, b - 1) * Deltad;
      fird += n * (Math.Pow(Delta, b) * (F + delta * Fd) + DeltaBd * delta * F);
    }
    return fird;
  }

  // Residual contribution to the adimensional free Helmholtz energy, tau
  //   derivative
  //
  //   Parameters
  //   ----------
  //   tau : float
  //     Inverse reduced temperature Tc/T, [-]
  //   delta : float
  //     Reduced density rho/rhoc, [-]
  //   coef : dict
  //     Dictionary with equation of state parameters
  //
  //   Returns
  //   -------
  //   firt : float
  //     .. math::
  //       \left.\frac{\partial \phi^r_{\tau}}{\partial \tau}\right|_{\delta}
  //
  //   References
  //   ----------
  //   IAPWS, Revised Release on the IAPWS Formulation 1995 for the
  //   Thermodynamic Properties of Ordinary Water Substance for General and
  //   Scientific Use, September 2016, Table 5
  //   http://www.iapws.org/relguide/IAPWS-95.html
  //

  public static object _phirt(object tau, object delta, object coef)
  {
    object b;
    object a;
    object g;
    object t;
    object d;
    object n;
    var firt = 0;
    // Polinomial terms
    var nr1 = coef.get("nr1", new List<object>());
    var d1 = coef.get("d1", new List<object>());
    var t1 = coef.get("t1", new List<object>());
    foreach (var _tup_1 in zip(nr1, d1, t1))
    {
      n = _tup_1.Item1;
      d = _tup_1.Item2;
      t = _tup_1.Item3;
      firt += n * t * Math.Pow(delta, d) * Math.Pow(tau, t - 1);
    }
    // Exponential terms
    var nr2 = coef.get("nr2", new List<object>());
    var d2 = coef.get("d2", new List<object>());
    var g2 = coef.get("gamma2", new List<object>());
    var t2 = coef.get("t2", new List<object>());
    var c2 = coef.get("c2", new List<object>());
    foreach (var _tup_2 in zip(nr2, d2, g2, t2, c2))
    {
      n = _tup_2.Item1;
      d = _tup_2.Item2;
      g = _tup_2.Item3;
      t = _tup_2.Item4;
      var c = _tup_2.Item5;
      firt += n * t * Math.Pow(delta, d) * Math.Pow(tau, t - 1) * Math.Exp(-g * Math.Pow(delta, c));
    }
    // Gaussian terms
    var nr3 = coef.get("nr3", new List<object>());
    var d3 = coef.get("d3", new List<object>());
    var t3 = coef.get("t3", new List<object>());
    var a3 = coef.get("alfa3", new List<object>());
    var e3 = coef.get("epsilon3", new List<object>());
    var b3 = coef.get("beta3", new List<object>());
    var g3 = coef.get("gamma3", new List<object>());
    foreach (var _tup_3 in zip(nr3, d3, t3, a3, e3, b3, g3))
    {
      n = _tup_3.Item1;
      d = _tup_3.Item2;
      t = _tup_3.Item3;
      a = _tup_3.Item4;
      var e = _tup_3.Item5;
      b = _tup_3.Item6;
      g = _tup_3.Item7;
      firt += n * Math.Pow(delta, d) * Math.Pow(tau, t) * Math.Exp(-a * Math.Pow(delta - e, 2) - b * Math.Pow(tau - g, 2)) * (t / tau - 2 * b * (tau - g));
    }
    // Non analitic terms
    var nr4 = coef.get("nr4", new List<object>());
    var a4 = coef.get("a4", new List<object>());
    var b4 = coef.get("b4", new List<object>());
    var Ai = coef.get("A", new List<object>());
    var Bi = coef.get("B", new List<object>());
    var Ci = coef.get("C", new List<object>());
    var Di = coef.get("D", new List<object>());
    var bt4 = coef.get("beta4", new List<object>());
    foreach (var _tup_4 in zip(nr4, a4, b4, Ai, Bi, Ci, Di, bt4))
    {
      n = _tup_4.Item1;
      a = _tup_4.Item2;
      b = _tup_4.Item3;
      var A = _tup_4.Item4;
      var B = _tup_4.Item5;
      var C = _tup_4.Item6;
      var D = _tup_4.Item7;
      var bt = _tup_4.Item8;
      var Tita = 1 - tau + A * Math.Pow(Math.Pow(delta - 1, 2), 0.5 / bt);
      var F = Math.Exp(-C * Math.Pow(delta - 1, 2) - D * Math.Pow(tau - 1, 2));
      var Ft = -2 * D * F * (tau - 1);
      var Delta = Math.Pow(Tita, 2) + B * Math.Pow(Math.Pow(delta - 1, 2), a);
      var DeltaBt = -2 * Tita * b * Math.Pow(Delta, b - 1);
      firt += n * delta * (DeltaBt * F + Math.Pow(Delta, b) * Ft);
    }
    return firt;
  }

  //
  //   General implementation of multiparameter equation of state. From this
  //   derived all child class specified per individual compounds
  //
  //   Parameters
  //   ----------
  //   T : float
  //     Temperature, [K]
  //   P : float
  //     Pressure, [MPa]
  //   rho : float
  //     Density, [kg/m³]
  //   v : float
  //     Specific volume, [m³/kg]
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //   s : float
  //     Specific entropy, [kJ/kgK]
  //   u : float
  //     Specific internal energy, [kJ/kg]
  //   x : float
  //     Vapor quality, [-]
  //
  //   l : float, optional
  //     Wavelength of light, for refractive index, [nm]
  //   rho0 : float, optional
  //     Initial value of density, to improve iteration, [kg/m³]
  //   T0 : float, optional
  //     Initial value of temperature, to improve iteration, [K]
  //   x0 : Initial value of vapor quality, necessary in bad input pair definition
  //     where there are two valid solution (T-h, T-s)
  //
  //   Notes
  //   -----
  //   * It needs two incoming properties of T, P, rho, h, s, u.
  //   * v as a alternate input parameter to rho
  //   * T-x, P-x, preferred input pair to specified a point in two phases region
  //
  //   The calculated instance has the following properties:
  //
  //     * P: Pressure, [MPa]
  //     * T: Temperature, [K]
  //     * x: Vapor quality, [-]
  //     * g: Specific Gibbs free energy, [kJ/kg]
  //     * a: Specific Helmholtz free energy, [kJ/kg]
  //     * v: Specific volume, [m³/kg]
  //     * r: Density, [kg/m³]
  //     * h: Specific enthalpy, [kJ/kg]
  //     * u: Specific internal energy, [kJ/kg]
  //     * s: Specific entropy, [kJ/kg·K]
  //     * cp: Specific isobaric heat capacity, [kJ/kg·K]
  //     * cv: Specific isochoric heat capacity, [kJ/kg·K]
  //     * cp_cv: Heat capacity ratio, [-]
  //     * Z: Compression factor, [-]
  //     * fi: Fugacity coefficient, [-]
  //     * f: Fugacity, [MPa]
  //     * gamma: Isoentropic exponent, [-]
  //
  //     * alfav: Isobaric cubic expansion coefficient, [1/K]
  //     * kappa: Isothermal compressibility, [1/MPa]
  //     * kappas: Adiabatic compresibility, [1/MPa]
  //     * alfap: Relative pressure coefficient, [1/K]
  //     * betap: Isothermal stress coefficient, [kg/m³]
  //     * joule: Joule-Thomson coefficient, [K/MPa]
  //
  //     * betas: Isoentropic temperature-pressure coefficient, [-]
  //     * Gruneisen: Gruneisen parameter, [-]
  //     * virialB: Second virial coefficient, [m³/kg]
  //     * virialC: Third virial coefficient, [m⁶/kg²]
  //     * dpdT_rho: Derivatives, dp/dT at constant rho, [MPa/K]
  //     * dpdrho_T: Derivatives, dp/drho at constant T, [MPa·m³/kg]
  //     * drhodT_P: Derivatives, drho/dT at constant P, [kg/m³·K]
  //     * drhodP_T: Derivatives, drho/dP at constant T, [kg/m³·MPa]
  //     * dhdT_rho: Derivatives, dh/dT at constant rho, [kJ/kg·K]
  //     * dhdP_T: Isothermal throttling coefficient, [kJ/kg·MPa]
  //     * dhdT_P: Derivatives, dh/dT at constant P, [kJ/kg·K]
  //     * dhdrho_T: Derivatives, dh/drho at constant T, [kJ·m³/kg²]
  //     * dhdrho_P: Derivatives, dh/drho at constant P, [kJ·m³/kg²]
  //     * dhdP_rho: Derivatives, dh/dP at constant rho, [kJ/kg·MPa]
  //     * kt: Isothermal Expansion Coefficient, [-]
  //     * ks: Adiabatic Compressibility, [1/MPa]
  //     * Ks: Adiabatic bulk modulus, [MPa]
  //     * Kt: Isothermal bulk modulus, [MPa]
  //
  //     * v0: Ideal specific volume, [m³/kg]
  //     * rho0: Ideal gas density, [kg/m³]
  //     * u0: Ideal specific internal energy, [kJ/kg]
  //     * h0: Ideal specific enthalpy, [kJ/kg]
  //     * s0: Ideal specific entropy, [kJ/kg·K]
  //     * a0: Ideal specific Helmholtz free energy, [kJ/kg]
  //     * g0: Ideal specific Gibbs free energy, [kJ/kg]
  //     * cp0: Ideal specific isobaric heat capacity, [kJ/kg·K]
  //     * cv0: Ideal specific isochoric heat capacity, [kJ/kg·K]
  //     * w0: Ideal speed of sound, [m/s]
  //     * gamma0: Ideal isoentropic exponent, [-]
  //
  //     * w: Speed of sound, [m/s]
  //     * mu: Dynamic viscosity, [Pa·s]
  //     * nu: Kinematic viscosity, [m²/s]
  //     * k: Thermal conductivity, [W/m·K]
  //     * alfa: Thermal diffusivity, [m²/s]
  //     * sigma: Surface tension, [N/m]
  //     * epsilon: Dielectric constant, [-]
  //     * n: Refractive index, [-]
  //     * Prandt: Prandtl number, [-]
  //     * Pr: Reduced Pressure, [-]
  //     * Tr: Reduced Temperature, [-]
  //     * Hvap: Vaporization heat, [kJ/kg]
  //     * Svap: Vaporization entropy, [kJ/kg·K]
  //
  //     * Z_rho: :math:`(Z-1)/\rho`, [m³/kg]
  //     * IntP: Internal pressure, [MPa]
  //     * invT: Negative reciprocal temperature, [1/K]
  //     * hInput: Specific heat input, [kJ/kg]
  //

  public class MEoS : _utils._fase
  {
    public string _mode;
    public int a;
    public int a0;
    public object cp0;
    public object cp0_cv;
    public object cv0;
    public int f;
    public int g;
    public object g0;
    public int gamma0;
    public object Gas;
    public int h;
    public object h0;
    public double Hvap;
    public int IntP;
    public int invT;
    public object Liquid;
    public object P;
    public object phase;
    public object Pr;
    public object R;
    public double rho;
    public double rho0;
    public int s;
    public object s0;
    public double sigma;
    public double Svap;
    public object T;
    public object Tr;
    public int u;
    public int u0;
    public int v;
    public int v0;
    public object virialB;
    public int virialC;
    public object x;
    public int Z;
    public int Z_rho;
    public object Zc;
    public double CP = null;
    public double _Pv = null;
    public double _rhoL = null;
    public double _rhoG = null;

    public Dictionary<string, double> kwargs = new Dictionary<object, object>
    {
      { "T", 0.0},
      { "P", 0.0},
      { "rho", 0.0},
      { "v", 0.0},
      { "h", null},
      { "s", null},
      { "u", null},
      { "x", null},
      { "l", 0.5893},
      { "rho0", null},
      { "T0", null},
      { "x0", 0.5}
    };

    public int status = 0;
    public string msg = "Undefined";

    public MEoS(Hashtable kwargs)
    {
      this.R = this._constants["R"] / this._constants.get("M", this.M);
      this.Zc = this.Pc / this.rhoc / this.R / this.Tc;
      this.kwargs = MEoS.kwargs.copy();
      this.@__call__(kwargs);
    }

    // Make instance callable to can add input parameter one to one
    public virtual object @__call__(Hashtable kwargs)
    {
      // Alternative rho input
      if(kwargs.Contains("rhom"))
      {
        kwargs["rho"] = kwargs["rhom"] * this.M;
        kwargs.Remove("rhom");
      }
      else if(kwargs.get("v", 0))
      {
        kwargs["rho"] = 1.0 / kwargs["v"];
        kwargs.Remove("v");
      }
      else if(kwargs.get("vm", 0))
      {
        kwargs["rho"] = this.M / kwargs["vm"];
        kwargs.Remove("vm");
      }
      this.kwargs.update(kwargs);
      if(this.calculable)
      {
        try
        {
          this.status = 1;
          this.calculo();
          this.msg = "";
        }
        catch(RuntimeError)
        {
          this.status = 0;
          this.msg = err.args[0];
          throw new RuntimeError();
        }
        if(this.name == "water" && 130 <= this.T < 273.15)
        {
          this.msg = "Extrapolated state";
          this.status = 3;
          Debug.Log("Using extrapolated values");
        }
        else if(this.name == "water" && 50 <= this.T < 130)
        {
          this.msg = "Extrapolated state using Low-Temperature extension";
          this.status = 3;
          Debug.Log("Using extrapolated values and Low-Temperatureextension");
        }
      }
    }

    // Check if inputs are enough to define state
    public object calculable
    {
      get
      {
        this._mode = "";
        if(this.kwargs["T"] && this.kwargs["P"]) this._mode = "TP";
        else if(this.kwargs["T"] && this.kwargs["rho"]) this._mode = "Trho";
        else if(this.kwargs["T"] && this.kwargs["h"] != null) this._mode = "Th";
        else if(this.kwargs["T"] && this.kwargs["s"] != null) this._mode = "Ts";
        else if(this.kwargs["T"] && this.kwargs["u"] != null) this._mode = "Tu";
        else if(this.kwargs["P"] && this.kwargs["rho"]) this._mode = "Prho";
        else if(this.kwargs["P"] && this.kwargs["h"] != null) this._mode = "Ph";
        else if(this.kwargs["P"] && this.kwargs["s"] != null) this._mode = "Ps";
        else if(this.kwargs["P"] && this.kwargs["u"] != null) this._mode = "Pu";
        else if(this.kwargs["rho"] && this.kwargs["h"] != null) this._mode = "rhoh";
        else if(this.kwargs["rho"] && this.kwargs["s"] != null) this._mode = "rhos";
        else if(this.kwargs["rho"] && this.kwargs["u"] != null) this._mode = "rhou";
        else if(this.kwargs["h"] != null && this.kwargs["s"] != null) this._mode = "hs";
        else if(this.kwargs["h"] != null && this.kwargs["u"] != null) this._mode = "hu";
        else if(this.kwargs["s"] != null && this.kwargs["u"] != null) this._mode = "su";
        else if(this.kwargs["T"] && this.kwargs["x"] != null) this._mode = "Tx";
        else if(this.kwargs["P"] && this.kwargs["x"] != null) this._mode = "Px";
        return this._mode != "";
      }
    }

    // Calculate procedure
    public virtual object calculo()
    {
      object region;
      object rGo;
      object rLo;
      object rhoG;
      object rhoL;
      object sol;
      object rhoGo;
      object rhoLo;
      object delta;
      object sv;
      object sl;
      object liquido;
      object vapor;
      object firtG;
      object firdG;
      object firtL;
      object firdL;
      object deltaG;
      object deltaL;
      object x0;
      object ideal;
      object rhov;
      object rhol;
      object Ps;
      object st0;
      object rhoo;
      object To;
      var T = this.kwargs["T"];
      var rho = this.kwargs["rho"];
      var P = this.kwargs["P"];
      var s = this.kwargs["s"];
      var h = this.kwargs["h"];
      var u = this.kwargs["u"];
      var x = this.kwargs["x"];
      // Initial values
      var T0 = this.kwargs["T0"];
      var rho0 = this.kwargs["rho0"];
      if(T0 || rho0)
      {
        To = T0;
        rhoo = rho0;
      }
      else if(this.name == "air")
      {
        To = 300;
        rhoo = 0.001;
      }
      else
      {
        try
        {
          st0 = iapws97.IAPWS97(this.kwargs);
        } catch (NotImplementedError)
        {
          To = 300;
          rhoo = 900;
        }
      }
      this.R = this._constants["R"] / this._constants.get("M", this.M);
      var rhoc = this._constants.get("rhoref", this.rhoc);
      var Tc = this._constants.get("Tref", this.Tc);
      object propiedades = null;
      if(!("Tx", "Px").Contains(this._mode))
      {
        // Method with iteration necessary to get x
        if(this._mode == "TP")
        {
          try
          {
            if(this.name == "air")
            {
              throw new ValueError();
            }
            st0 = iapws97.IAPWS97(this.kwargs);
            rhoo = st0.rho;
          } catch (NotImplementedError)
          {
            if(rho0)
            {
              rhoo = rho0;
            }
            else if(T < this.Tc && P < this.Pc && this._Vapor_Pressure(T) < P)
            {
              rhoo = this._Liquid_Density(T);
            }
            else if(T < this.Tc && P < this.Pc)
            {
              rhoo = this._Vapor_Density(T);
            }
            else
            {
              rhoo = this.rhoc * 3;
            }
          } catch (ValueError)
          {
            rhoo = 0.001;
          }

          Func<double, double> f = rho =>
          {
            var delta = rho / rhoc;
            var tau = Tc / T;
            var fird = _phird(tau, delta, this._constants);
            var Po = (1 + delta * fird) * this.R * T * rho;
            return Po - P * 1000;
          };

          rho = fsolve(f, rhoo)[0];

          // Calculate quality
          if(T > this.Tc)
          {
            x = 1;
          }
          else
          {
            Ps = this._Vapor_Pressure(T);
            if(Ps * 0.95 < P < Ps * 1.05)
            {
              var _tup_1 = this._saturation(T);
              rhol = _tup_1.Item1;
              rhov = _tup_1.Item2;
              Ps = _tup_1.Item3;
              Ps *= 0.001;
            }
            if(Ps > P) { x = 1; }
            else { x = 0; }
          }
        }
        else if(this._mode == "Th")
        {
          var tau = Tc / T;
          ideal = this._phi0(tau, 1);
          var fiot = ideal["fiot"];

          Func<double, double> f = rho =>
          {
            var delta = rho / rhoc;
            var fird = _phird(tau, delta, this._constants);
            var firt = _phirt(tau, delta, this._constants);
            var ho = this.R * T * (1 + tau * (fiot + firt) + delta * fird);
            return ho - h;
          };

          if(T >= this.Tc)
          {
            rhoo = this.rhoc;
            rho = fsolve(f, rhoo)[0];
          }
          else
          {
            x0 = this.kwargs["x0"];
            rhov = this._Vapor_Density(T);
            rhol = this._Liquid_Density(T);
            deltaL = rhol / rhoc;
            deltaG = rhov / rhoc;
            firdL = _phird(tau, deltaL, this._constants);
            firtL = _phirt(tau, deltaL, this._constants);
            firdG = _phird(tau, deltaG, this._constants);
            firtG = _phirt(tau, deltaG, this._constants);
            var hl = this.R * T * (1 + tau * (fiot + firtL) + deltaL * firdL);
            var hv = this.R * T * (1 + tau * (fiot + firtG) + deltaG * firdG);
            if(!(0, 1).Contains(x0) && hl <= h <= hv)
            {
              var _tup_2 = this._saturation(T);
              rhol = _tup_2.Item1;
              rhov = _tup_2.Item2;
              Ps = _tup_2.Item3;
              vapor = this._Helmholtz(rhov, T);
              liquido = this._Helmholtz(rhol, T);
              hv = vapor["h"];
              hl = liquido["h"];
              x = (h - hl) / (hv - hl);
              rho = 1 / (x / rhov + (1 - x) / rhol);
              P = Ps / 1000;
            }
            else
            {
              if(h > hv) { rhoo = rhov; }
              else { rhoo = rhol; }
              rho = fsolve(f, rhoo)[0];
            }
          }
        }
        else if(this._mode == "Ts")
        {
          tau = Tc / T;

          Func<double, double> f = rho =>
          {
            if(rho < 0) { rho = 1E-20; }
            var delta = rho / rhoc;
            var ideal = this._phi0(tau, delta);
            var fio = ideal["fio"];
            var fiot = ideal["fiot"];
            var fir = _phir(tau, delta, this._constants);
            var firt = _phirt(tau, delta, this._constants);
            var so = this.R * (tau * (fiot + firt) - fio - fir);
            return so - s;
          };

          if(T >= this.Tc)
          {
            rhoo = this.rhoc;
            rho = fsolve(f, rhoo)[0];
          }
          else
          {
            rhov = this._Vapor_Density(T);
            rhol = this._Liquid_Density(T);
            deltaL = rhol / rhoc;
            deltaG = rhov / rhoc;
            var idealL = this._phi0(tau, deltaL);
            var idealG = this._phi0(tau, deltaG);
            var fioL = idealL["fio"];
            var fioG = idealG["fio"];
            fiot = idealL["fiot"];
            var firL = _phir(tau, deltaL, this._constants);
            firtL = _phirt(tau, deltaL, this._constants);
            sl = this.R * (tau * (fiot + firtL) - fioL - firL);
            var firG = _phir(tau, deltaG, this._constants);
            firtG = _phirt(tau, deltaG, this._constants);
            sv = this.R * (tau * (fiot + firtG) - fioG - firG);
            if(sl <= s <= sv)
            {
              var _tup_3 = this._saturation(T);
              rhol = _tup_3.Item1;
              rhov = _tup_3.Item2;
              Ps = _tup_3.Item3;
              vapor = this._Helmholtz(rhov, T);
              liquido = this._Helmholtz(rhol, T);
              sv = vapor["s"];
              sl = liquido["s"];
              x = (s - sl) / (sv - sl);
              rho = 1 / (x / rhov + (1 - x) / rhol);
              P = Ps / 1000;
            }
            else
            {
              if(s > sv)
              {
                rhoo = rhov;
              }
              else
              {
                rhoo = rhol;
              }
              rho = fsolve(f, rhoo)[0];
            }
          }
        }
        else if(this._mode == "Tu")
        {
          tau = Tc / T;
          ideal = this._phi0(tau, 1);
          fiot = ideal["fiot"];

          Func<double, double> f = rho =>
          {
            var delta = rho / rhoc;
            var fird = _phird(tau, delta, this._constants);
            var firt = _phirt(tau, delta, this._constants);
            var Po = (1 + delta * fird) * this.R * T * rho;
            var ho = this.R * T * (1 + tau * (fiot + firt) + delta * fird);
            return ho - Po / rho - u;
          };

          if(T >= this.Tc)
          {
            rhoo = this.rhoc;
            rho = fsolve(f, rhoo)[0];
          }
          else
          {
            rhov = this._Vapor_Density(T);
            rhol = this._Liquid_Density(T);
            deltaL = rhol / rhoc;
            deltaG = rhov / rhoc;
            firdL = _phird(tau, deltaL, this._constants);
            firtL = _phirt(tau, deltaL, this._constants);
            firdG = _phird(tau, deltaG, this._constants);
            firtG = _phirt(tau, deltaG, this._constants);
            var PoL = (1 + deltaL * firdL) * this.R * T * rhol;
            var PoG = (1 + deltaG * firdG) * this.R * T * rhov;
            var hoL = this.R * T * (1 + tau * (fiot + firtL) + deltaL * firdL);
            var hoG = this.R * T * (1 + tau * (fiot + firtG) + deltaG * firdG);
            var uv = hoG - PoG / rhov;
            var ul = hoL - PoL / rhol;
            if(ul <= u <= uv)
            {
              var _tup_4 = this._saturation(T);
              rhol = _tup_4.Item1;
              rhov = _tup_4.Item2;
              Ps = _tup_4.Item3;
              vapor = this._Helmholtz(rhov, T);
              liquido = this._Helmholtz(rhol, T);
              uv = vapor["h"] - vapor["P"] / rhov;
              ul = liquido["h"] - liquido["P"] / rhol;
              x = (u - ul) / (uv - ul);
              rho = 1 / (x / rhov - (1 - x) / rhol);
              P = Ps / 1000;
            }
            else
            {
              if(u > uv) { rhoo = rhov; }
              else { rhoo = rhol; }
              rho = fsolve(f, rhoo)[0];
            }
          }
        }
        else if(this._mode == "Prho")
        {
          delta = rho / rhoc;

          Func<double, double> f = T =>
          {
            var tau = Tc / T;
            var fird = _phird(tau, delta, this._constants);
            var Po = (1 + delta * fird) * this.R * T * rho;
            return Po - P * 1000;
          };

          T = fsolve(f, To)[0];
          rhol = this._Liquid_Density(T);
          rhov = this._Vapor_Density(T);
          if(T == To || rhov <= rho <= rhol)
          {

            Func<double, double> f = parr =>
            {
              var _tup_1 = parr;
              var T = _tup_1.Item1;
              var rhol = _tup_1.Item2;
              var rhog = _tup_1.Item3;
              var tau = Tc / T;
              var deltaL = rhol / this.rhoc;
              var deltaG = rhog / this.rhoc;
              var firL = _phir(tau, deltaL, this._constants);
              var firdL = _phird(tau, deltaL, this._constants);
              var firG = _phir(tau, deltaG, this._constants);
              var firdG = _phird(tau, deltaG, this._constants);
              var Jl = rhol * (1 + deltaL * firdL);
              var Jv = rhog * (1 + deltaG * firdG);
              var K = firL - firG;
              var Ps = this.R * T * rhol * rhog / (rhol - rhog) * (K + Math.Log(rhol / rhog));
              return (Jl - Jv, Jl * (1 / rhog - 1 / rhol) - Math.Log(rhol / rhog) - K, Ps - P * 1000);
            };

            foreach (var to in new List<int> { To, 300, 400, 500, 600 })
            {
              rhoLo = this._Liquid_Density(to);
              rhoGo = this._Vapor_Density(to);
              sol = fsolve(f, new List<object> { to, rhoLo, rhoGo }, full_output: true);
              var _tup_5 = sol[0];
              T = _tup_5.Item1;
              rhoL = _tup_5.Item2;
              rhoG = _tup_5.Item3;
              x = (1.0 / rho - 1 / rhoL) / (1 / rhoG - 1 / rhoL);
              if(sol[2] == 1 && 0 <= x <= 1 && abs(sol[1]["fvec"]).Sum() < 1E-05)
              {
                break;
              }
            }
            if(abs(sol[1]["fvec"]).Sum() > 1E-05)
            {
              throw new RuntimeError(sol[3]);
            }
            liquido = this._Helmholtz(rhoL, T);
            vapor = this._Helmholtz(rhoG, T);
            P = this.R * T * rhoL * rhoG / (rhoL - rhoG) * (liquido["fir"] - vapor["fir"] + Math.Log(rhoL / rhoG)) / 1000;
          }
        }
        else if(this._mode == "Ph")
        {

          Func<double, double> funcion = parr =>
          {
            var _tup_1 = parr;
            var rho = _tup_1.Item1;
            var T = _tup_1.Item2;
            var delta = rho / rhoc;
            var tau = Tc / T;
            var ideal = this._phi0(tau, delta);
            var fiot = ideal["fiot"];
            var fird = _phird(tau, delta, this._constants);
            var firt = _phirt(tau, delta, this._constants);
            var Po = (1 + delta * fird) * this.R * T * rho;
            var ho = this.R * T * (1 + tau * (fiot + firt) + delta * fird);
            return Tuple.Create(Po - P * 1000, ho - h);
          };

          var _tup_6 = fsolve(funcion, new List<object> { rhoo, To });
          rho = _tup_6.Item1;
          T = _tup_6.Item2;
          rhol = this._Liquid_Density(T);
          rhov = this._Vapor_Density(T);
          if(rho == rhoo || rhov <= rho <= rhol)
          {

            Func<double, double> f = parr =>
            {
              var _tup_1 = parr;
              var T = _tup_1.Item1;
              var rhol = _tup_1.Item2;
              var rhog = _tup_1.Item3;
              var x = _tup_1.Item4;
              var tau = Tc / T;
              var deltaL = rhol / this.rhoc;
              var deltaG = rhog / this.rhoc;
              var ideal = this._phi0(tau, deltaL);
              var fiot = ideal["fiot"];
              var firL = _phir(tau, deltaL, this._constants);
              var firdL = _phird(tau, deltaL, this._constants);
              var firtL = _phirt(tau, deltaL, this._constants);
              var hoL = this.R * T * (1 + tau * (fiot + firtL) + deltaL * firdL);
              var firG = _phir(tau, deltaG, this._constants);
              var firdG = _phird(tau, deltaG, this._constants);
              var firtG = _phirt(tau, deltaG, this._constants);
              var hoG = this.R * T * (1 + tau * (fiot + firtG) + deltaG * firdG);
              var Jl = rhol * (1 + deltaL * firdL);
              var Jv = rhog * (1 + deltaG * firdG);
              var K = firL - firG;
              var Ps = this.R * T * rhol * rhog / (rhol - rhog) * (K + Math.Log(rhol / rhog));
              return (Jl - Jv, Jl * (1 / rhog - 1 / rhol) - Math.Log(rhol / rhog) - K, hoL * (1 - x) + hoG * x - h, Ps - P * 1000);
            };

            foreach(var to in new List<int> { To, 300, 400, 500, 600 })
            {
              rLo = this._Liquid_Density(to);
              rGo = this._Vapor_Density(to);
              sol = fsolve(f, new List<object>
                { to, rLo, rGo, 0.5 },
                full_output: true);
              var _tup_7 = sol[0];
              T = _tup_7.Item1;
              rhoL = _tup_7.Item2;
              rhoG = _tup_7.Item3;
              x = _tup_7.Item4;
              if(sol[2] == 1 && 0 <= x <= 1 && abs(sol[1]["fvec"]).Sum() < 1E-05)
              {
                break;
              }
            }
            if(abs(sol[1]["fvec"]).Sum() > 1E-05)
            {
              throw new RuntimeError(sol[3]);
            }
            liquido = this._Helmholtz(rhoL, T);
            vapor = this._Helmholtz(rhoG, T);
            P = this.R * T * rhoL * rhoG / (rhoL - rhoG) * (liquido["fir"] - vapor["fir"] + Math.Log(rhoL / rhoG)) / 1000;
          }
        }
        else if(this._mode == "Ps")
        {
          try
          {
            x0 = st0.x;
          } catch (NameError)
          {
            x0 = null;
          }
          if(x0 == null || x0 == 0 || x0 == 1)
          {

            Func<double, double> f = parr =>
            {
              var _tup_1 = parr;
              var rho = _tup_1.Item1;
              var T = _tup_1.Item2;
              var delta = rho / rhoc;
              var tau = Tc / T;
              var ideal = this._phi0(tau, delta);
              var fio = ideal["fio"];
              var fiot = ideal["fiot"];
              var fird = _phird(tau, delta, this._constants);
              var fir = _phir(tau, delta, this._constants);
              var firt = _phirt(tau, delta, this._constants);
              var Po = (1 + delta * fird) * this.R * T * rho;
              var so = this.R * (tau * (fiot + firt) - fio - fir);
              return Tuple.Create(Po - P * 1000, so - s);
            };

            var _tup_8 = fsolve(f, new List<object> { rhoo, To });
            rho = _tup_8.Item1;
            T = _tup_8.Item2;
          }
          else
          {

            Func<double, double> funcion = parr =>
            {
              var _tup_1 = parr;
              var rho = _tup_1.Item1;
              var T = _tup_1.Item2;
              var _tup_2 = this._saturation(T);
              var rhol = _tup_2.Item1;
              var rhov = _tup_2.Item2;
              var Ps = _tup_2.Item3;
              var vapor = this._Helmholtz(rhov, T);
              var liquido = this._Helmholtz(rhol, T);
              var x = (1.0 / rho - 1 / rhol) / (1 / rhov - 1 / rhol);
              return Tuple.Create(Ps - P * 1000, vapor["s"] * x + liquido["s"] * (1 - x) - s);
            };

            var _tup_9 = fsolve(funcion, new List<double> { 2.0, 500.0 });
            rho = _tup_9.Item1;
            T = _tup_9.Item2;
            var _tup_10 = this._saturation(T);
            rhol = _tup_10.Item1;
            rhov = _tup_10.Item2;
            Ps = _tup_10.Item3;
            vapor = this._Helmholtz(rhov, T);
            liquido = this._Helmholtz(rhol, T);
            sv = vapor["s"];
            sl = liquido["s"];
            x = (s - sl) / (sv - sl);
          }
        }
        else if(this._mode == "Pu")
        {

          Func<double, double> f = parr =>
          {
            var _tup_1 = parr;
            var rho = _tup_1.Item1;
            var T = _tup_1.Item2;
            var delta = rho / rhoc;
            var tau = Tc / T;
            var ideal = this._phi0(tau, delta);
            var fiot = ideal["fiot"];
            var fird = _phird(tau, delta, this._constants);
            var firt = _phirt(tau, delta, this._constants);
            var Po = (1 + delta * fird) * this.R * T * rho;
            var ho = this.R * T * (1 + tau * (fiot + firt) + delta * fird);
            return Tuple.Create(ho - Po / rho - u, Po - P * 1000);
          };

          sol = fsolve(f, new List<object> { rhoo, To }, full_output: true);
          var _tup_11 = sol[0];
          rho = _tup_11.Item1;
          T = _tup_11.Item2;
          rhol = this._Liquid_Density(T);
          rhov = this._Vapor_Density(T);
          if(rho == rhoo || sol[2] != 1)
          {

            Func<double, double> f = parr =>
            {
              var _tup_1 = parr;
              var T = _tup_1.Item1;
              var rhol = _tup_1.Item2;
              var rhog = _tup_1.Item3;
              var x = _tup_1.Item4;
              var tau = Tc / T;
              var deltaL = rhol / this.rhoc;
              var deltaG = rhog / this.rhoc;
              var ideal = this._phi0(tau, deltaL);
              var fiot = ideal["fiot"];
              var firL = _phir(tau, deltaL, this._constants);
              var firdL = _phird(tau, deltaL, this._constants);
              var firtL = _phirt(tau, deltaL, this._constants);
              var hoL = this.R * T * (1 + tau * (fiot + firtL) + deltaL * firdL);
              var firG = _phir(tau, deltaG, this._constants);
              var firdG = _phird(tau, deltaG, this._constants);
              var firtG = _phirt(tau, deltaG, this._constants);
              var hoG = this.R * T * (1 + tau * (fiot + firtG) + deltaG * firdG);
              var Jl = rhol * (1 + deltaL * firdL);
              var Jv = rhog * (1 + deltaG * firdG);
              var K = firL - firG;
              var Ps = this.R * T * rhol * rhog / (rhol - rhog) * (K + Math.Log(rhol / rhog));
              var vu = hoG - Ps / rhog;
              var lu = hoL - Ps / rhol;
              return (Jl - Jv, Jl * (1 / rhog - 1 / rhol) - Math.Log(rhol / rhog) - K, lu * (1 - x) + vu * x - u, Ps - P * 1000);
            };

            foreach (var to in new List<int> { To, 300, 400, 500, 600 })
            {
              rLo = this._Liquid_Density(to);
              rGo = this._Vapor_Density(to);
              sol = fsolve(f, new List<object> { to, rLo, rGo, 0.5 }, full_output: true);
              var _tup_12 = sol[0];
              T = _tup_12.Item1;
              rhoL = _tup_12.Item2;
              rhoG = _tup_12.Item3;
              x = _tup_12.Item4;
              if(sol[2] == 1 && 0 <= x <= 1 && abs(sol[1]["fvec"]).Sum() < 1E-05)
              {
                break;
              }
            }
            if(abs(sol[1]["fvec"]).Sum() > 1E-05)
            {
              throw new RuntimeError(sol[3]);
            }
            liquido = this._Helmholtz(rhoL, T);
            vapor = this._Helmholtz(rhoG, T);
            P = this.R * T * rhoL * rhoG / (rhoL - rhoG) * (liquido["fir"] - vapor["fir"] + Math.Log(rhoL / rhoG)) / 1000;
          }
        }
        else if(this._mode == "rhoh")
        {
          delta = rho / rhoc;

          Func<double, double> f = T =>
          {
            var tau = Tc / T;
            var ideal = this._phi0(tau, delta);
            var fiot = ideal["fiot"];
            var fird = _phird(tau, delta, this._constants);
            var firt = _phirt(tau, delta, this._constants);
            var ho = this.R * T * (1 + tau * (fiot + firt) + delta * fird);
            return ho - h;
          };

          T = fsolve(f, To)[0];
          rhol = this._Liquid_Density(T);
          rhov = this._Vapor_Density(T);
          if(T == To || rhov <= rho <= rhol)
          {

            Func<double, double> f = parr =>
            {
              var _tup_1 = parr;
              var T = _tup_1.Item1;
              var rhol = _tup_1.Item2;
              var rhog = _tup_1.Item3;
              var tau = Tc / T;
              var deltaL = rhol / this.rhoc;
              var deltaG = rhog / this.rhoc;
              var ideal = this._phi0(tau, deltaL);
              var fiot = ideal["fiot"];
              var firL = _phir(tau, deltaL, this._constants);
              var firdL = _phird(tau, deltaL, this._constants);
              var firtL = _phirt(tau, deltaL, this._constants);
              var hoL = this.R * T * (1 + tau * (fiot + firtL) + deltaL * firdL);
              var firG = _phir(tau, deltaG, this._constants);
              var firdG = _phird(tau, deltaG, this._constants);
              var firtG = _phirt(tau, deltaG, this._constants);
              var hoG = this.R * T * (1 + tau * (fiot + firtG) + deltaG * firdG);
              var Jl = rhol * (1 + deltaL * firdL);
              var Jv = rhog * (1 + deltaG * firdG);
              var K = firL - firG;
              var x = (1.0 / rho - 1 / rhol) / (1 / rhog - 1 / rhol);
              return (Jl - Jv, Jl * (1 / rhog - 1 / rhol) - Math.Log(rhol / rhog) - K, hoL * (1 - x) + hoG * x - h);
            };

            foreach (var to in new List<int> { To, 300, 400, 500, 600 })
            {
              rhoLo = this._Liquid_Density(to);
              rhoGo = this._Vapor_Density(to);
              sol = fsolve(f, new List<object> { to, rhoLo, rhoGo }, full_output: true);
              var _tup_13 = sol[0];
              T = _tup_13.Item1;
              rhoL = _tup_13.Item2;
              rhoG = _tup_13.Item3;
              x = (1.0 / rho - 1 / rhoL) / (1 / rhoG - 1 / rhoL);
              if(sol[2] == 1 && 0 <= x <= 1 && abs(sol[1]["fvec"]).Sum() < 1E-05)
              {
                break;
              }
            }
            if(abs(sol[1]["fvec"]).Sum() > 1E-05)
            {
              throw new RuntimeError(sol[3]);
            }
            liquido = this._Helmholtz(rhoL, T);
            vapor = this._Helmholtz(rhoG, T);
            P = this.R * T * rhoL * rhoG / (rhoL - rhoG) * (liquido["fir"] - vapor["fir"] + Math.Log(rhoL / rhoG)) / 1000;
          }
        }
        else if(this._mode == "rhos")
        {
          delta = rho / rhoc;

          Func<double, double> f = T =>
          {
            var tau = Tc / T;
            var ideal = this._phi0(tau, delta);
            var fio = ideal["fio"];
            var fiot = ideal["fiot"];
            var fir = _phir(tau, delta, this._constants);
            var firt = _phirt(tau, delta, this._constants);
            var so = this.R * (tau * (fiot + firt) - fio - fir);
            return so - s;
          };

          T = fsolve(f, To)[0];
          rhol = this._Liquid_Density(T);
          rhov = this._Vapor_Density(T);
          if(T == To || rhov <= rho <= rhol)
          {

            Func<double, double> f = parr =>
            {
              var _tup_1 = parr;
              var T = _tup_1.Item1;
              var rhol = _tup_1.Item2;
              var rhog = _tup_1.Item3;
              var tau = Tc / T;
              var deltaL = rhol / this.rhoc;
              var deltaG = rhog / this.rhoc;
              var idealL = this._phi0(tau, deltaL);
              var fioL = idealL["fio"];
              var fiot = idealL["fiot"];
              var idealG = this._phi0(tau, deltaG);
              var fioG = idealG["fio"];
              var firL = _phir(tau, deltaL, this._constants);
              var firdL = _phird(tau, deltaL, this._constants);
              var firtL = _phirt(tau, deltaL, this._constants);
              var soL = this.R * (tau * (fiot + firtL) - fioL - firL);
              var firG = _phir(tau, deltaG, this._constants);
              var firdG = _phird(tau, deltaG, this._constants);
              var firtG = _phirt(tau, deltaG, this._constants);
              var soG = this.R * (tau * (fiot + firtG) - fioG - firG);
              var Jl = rhol * (1 + deltaL * firdL);
              var Jv = rhog * (1 + deltaG * firdG);
              var K = firL - firG;
              var x = (1.0 / rho - 1 / rhol) / (1 / rhog - 1 / rhol);
              return (Jl - Jv, Jl * (1 / rhog - 1 / rhol) - Math.Log(rhol / rhog) - K, soL * (1 - x) + soG * x - s);
            };

            foreach (var to in new List<int> { To, 300, 400, 500, 600 })
            {
              rhoLo = this._Liquid_Density(to);
              rhoGo = this._Vapor_Density(to);
              sol = fsolve(f, new List<object>
              {
                to,
                rhoLo,
                rhoGo
              }, full_output: true);
              var _tup_14 = sol[0];
              T = _tup_14.Item1;
              rhoL = _tup_14.Item2;
              rhoG = _tup_14.Item3;
              x = (1.0 / rho - 1 / rhoL) / (1 / rhoG - 1 / rhoL);
              if(sol[2] == 1 && 0 <= x <= 1 && abs(sol[1]["fvec"]).Sum() < 1E-05)
              {
                break;
              }
            }
            if(abs(sol[1]["fvec"]).Sum() > 1E-05)
            {
              throw new RuntimeError(sol[3]);
            }
            liquido = this._Helmholtz(rhoL, T);
            vapor = this._Helmholtz(rhoG, T);
            P = this.R * T * rhoL * rhoG / (rhoL - rhoG) * (liquido["fir"] - vapor["fir"] + Math.Log(rhoL / rhoG)) / 1000;
          }
        }
        else if(this._mode == "rhou")
        {
          delta = rho / rhoc;

          Func<double, double> f = T =>
          {
            var tau = Tc / T;
            var ideal = this._phi0(tau, delta);
            var fiot = ideal["fiot"];
            var fird = _phird(tau, delta, this._constants);
            var firt = _phirt(tau, delta, this._constants);
            var Po = (1 + delta * fird) * this.R * T * rho;
            var ho = this.R * T * (1 + tau * (fiot + firt) + delta * fird);
            return ho - Po / rho - u;
          };

          T = fsolve(f, To)[0];
          rhol = this._Liquid_Density(T);
          rhov = this._Vapor_Density(T);
          if(T == To || rhov <= rho <= rhol)
          {

            Func<double, double> f = parr =>
            {
              var _tup_1 = parr;
              var T = _tup_1.Item1;
              var rhol = _tup_1.Item2;
              var rhog = _tup_1.Item3;
              var tau = Tc / T;
              var deltaL = rhol / this.rhoc;
              var deltaG = rhog / this.rhoc;
              var ideal = this._phi0(tau, deltaL);
              var fiot = ideal["fiot"];
              var firL = _phir(tau, deltaL, this._constants);
              var firdL = _phird(tau, deltaL, this._constants);
              var firtL = _phirt(tau, deltaL, this._constants);
              var hoL = this.R * T * (1 + tau * (fiot + firtL) + deltaL * firdL);
              var firG = _phir(tau, deltaG, this._constants);
              var firdG = _phird(tau, deltaG, this._constants);
              var firtG = _phirt(tau, deltaG, this._constants);
              var hoG = this.R * T * (1 + tau * (fiot + firtG) + deltaG * firdG);
              var Jl = rhol * (1 + deltaL * firdL);
              var Jv = rhog * (1 + deltaG * firdG);
              var K = firL - firG;
              var x = (1.0 / rho - 1 / rhol) / (1 / rhog - 1 / rhol);
              var Ps = this.R * T * rhol * rhog / (rhol - rhog) * (K + Math.Log(rhol / rhog));
              var vu = hoG - Ps / rhog;
              var lu = hoL - Ps / rhol;
              return (Jl - Jv, Jl * (1 / rhog - 1 / rhol) - Math.Log(rhol / rhog) - K, lu * (1 - x) + vu * x - u);
            };

            foreach (var to in new List<int> { To, 300, 400, 500, 600 })
            {
              rhoLo = this._Liquid_Density(to);
              rhoGo = this._Vapor_Density(to);
              sol = fsolve(f, new List<object> { to, rhoLo, rhoGo }, full_output: true);
              var _tup_15 = sol[0];
              T = _tup_15.Item1;
              rhoL = _tup_15.Item2;
              rhoG = _tup_15.Item3;
              x = (1.0 / rho - 1 / rhoL) / (1 / rhoG - 1 / rhoL);
              if(sol[2] == 1 && 0 <= x <= 1 && abs(sol[1]["fvec"]).Sum() < 1E-05)
              {
                break;
              }
            }
            if(abs(sol[1]["fvec"]).Sum() > 1E-05)
            {
              throw new RuntimeError(sol[3]);
            }
            liquido = this._Helmholtz(rhoL, T);
            vapor = this._Helmholtz(rhoG, T);
            P = this.R * T * rhoL * rhoG / (rhoL - rhoG) * (liquido["fir"] - vapor["fir"] + Math.Log(rhoL / rhoG)) / 1000;
          }
        }
        else if(this._mode == "hs")
        {

          Func<double, Tuple> f = parr =>
          {
            var _tup_1 = parr;
            var rho = _tup_1.Item1;
            var T = _tup_1.Item2;
            var delta = rho / rhoc;
            var tau = Tc / T;
            var ideal = this._phi0(tau, delta);
            var fio = ideal["fio"];
            var fiot = ideal["fiot"];
            var fird = _phird(tau, delta, this._constants);
            var fir = _phir(tau, delta, this._constants);
            var firt = _phirt(tau, delta, this._constants);
            var ho = this.R * T * (1 + tau * (fiot + firt) + delta * fird);
            var so = this.R * (tau * (fiot + firt) - fio - fir);
            return Tuple.Create(ho - h, so - s);
          };

          var _tup_16 = fsolve(f, new List<object> { rhoo, To });
          rho = _tup_16.Item1;
          T = _tup_16.Item2;
          rhol = this._Liquid_Density(T);
          rhov = this._Vapor_Density(T);
          if(rhov <= rho <= rhol)
          {

            Func<double, double> f = parr =>
            {
              var _tup_1 = parr;
              var T = _tup_1.Item1;
              var rhol = _tup_1.Item2;
              var rhog = _tup_1.Item3;
              var x = _tup_1.Item4;
              var tau = Tc / T;
              var deltaL = rhol / this.rhoc;
              var deltaG = rhog / this.rhoc;
              var idealL = this._phi0(tau, deltaL);
              var fiot = idealL["fiot"];
              var fioL = idealL["fio"];
              var idealG = this._phi0(tau, deltaG);
              var fioG = idealG["fio"];
              var firL = _phir(tau, deltaL, this._constants);
              var firdL = _phird(tau, deltaL, this._constants);
              var firtL = _phirt(tau, deltaL, this._constants);
              var hoL = this.R * T * (1 + tau * (fiot + firtL) + deltaL * firdL);
              var soL = this.R * (tau * (fiot + firtL) - fioL - firL);
              var firG = _phir(tau, deltaG, this._constants);
              var firdG = _phird(tau, deltaG, this._constants);
              var firtG = _phirt(tau, deltaG, this._constants);
              var hoG = this.R * T * (1 + tau * (fiot + firtG) + deltaG * firdG);
              var soG = this.R * (tau * (fiot + firtG) - fioG - firG);
              var Jl = rhol * (1 + deltaL * firdL);
              var Jv = rhog * (1 + deltaG * firdG);
              var K = firL - firG;
              return (Jl - Jv, Jl * (1 / rhog - 1 / rhol) - Math.Log(rhol / rhog) - K, hoL * (1 - x) + hoG * x - h, soL * (1 - x) + soG * x - s);
            };

            foreach (var to in new List<int> { To, 300, 400, 500, 600 })
            {
              rLo = this._Liquid_Density(to);
              rGo = this._Vapor_Density(to);
              sol = fsolve(f, new List<object> { to, rLo, rGo, 0.5 }, full_output: true);
              var _tup_17 = sol[0];
              T = _tup_17.Item1;
              rhoL = _tup_17.Item2;
              rhoG = _tup_17.Item3;
              x = _tup_17.Item4;
              if(sol[2] == 1 && 0 <= x <= 1 && abs(sol[1]["fvec"]).Sum() < 1E-05)
              {
                break;
              }
            }
            if(abs(sol[1]["fvec"]).Sum() > 1E-05)
            {
              throw new RuntimeError(sol[3]);
            }
            liquido = this._Helmholtz(rhoL, T);
            vapor = this._Helmholtz(rhoG, T);
            P = this.R * T * rhoL * rhoG / (rhoL - rhoG) * (liquido["fir"] - vapor["fir"] + Math.Log(rhoL / rhoG)) / 1000;
          }
        }
        else if(this._mode == "hu")
        {

          Func<double, Tuple> f = parr =>
          {
            var _tup_1 = parr;
            var rho = _tup_1.Item1;
            var T = _tup_1.Item2;
            var delta = rho / rhoc;
            var tau = Tc / T;
            var ideal = this._phi0(tau, delta);
            var fiot = ideal["fiot"];
            var fird = _phird(tau, delta, this._constants);
            var firt = _phirt(tau, delta, this._constants);
            var Po = (1 + delta * fird) * this.R * T * rho;
            var ho = this.R * T * (1 + tau * (fiot + firt) + delta * fird);
            return Tuple.Create(ho - Po / rho - u, ho - h);
          };

          sol = fsolve(f, new List<object> { rhoo, To }, full_output: true);
          var _tup_18 = sol[0];
          rho = _tup_18.Item1;
          T = _tup_18.Item2;
          rhol = this._Liquid_Density(T);
          rhov = this._Vapor_Density(T);
          if(sol[2] != 1 || rhov <= rho <= rhol)
          {

            Func<double, double> f = parr =>
            {
              var _tup_1 = parr;
              var T = _tup_1.Item1;
              var rhol = _tup_1.Item2;
              var rhog = _tup_1.Item3;
              var x = _tup_1.Item4;
              var tau = Tc / T;
              var deltaL = rhol / this.rhoc;
              var deltaG = rhog / this.rhoc;
              var ideal = this._phi0(tau, deltaL);
              var fiot = ideal["fiot"];
              var firL = _phir(tau, deltaL, this._constants);
              var firdL = _phird(tau, deltaL, this._constants);
              var firtL = _phirt(tau, deltaL, this._constants);
              var hoL = this.R * T * (1 + tau * (fiot + firtL) + deltaL * firdL);
              var firG = _phir(tau, deltaG, this._constants);
              var firdG = _phird(tau, deltaG, this._constants);
              var firtG = _phirt(tau, deltaG, this._constants);
              var hoG = this.R * T * (1 + tau * (fiot + firtG) + deltaG * firdG);
              var Jl = rhol * (1 + deltaL * firdL);
              var Jv = rhog * (1 + deltaG * firdG);
              var K = firL - firG;
              var Ps = this.R * T * rhol * rhog / (rhol - rhog) * (K + Math.Log(rhol / rhog));
              var vu = hoG - Ps / rhog;
              var lu = hoL - Ps / rhol;
              return (Jl - Jv, Jl * (1 / rhog - 1 / rhol) - Math.Log(rhol / rhog) - K, hoL * (1 - x) + hoG * x - h, lu * (1 - x) + vu * x - u);
            };

            foreach (var to in new List<int> { To, 300, 400, 500, 600 })
            {
              rLo = this._Liquid_Density(to);
              rGo = this._Vapor_Density(to);
              sol = fsolve(f, new List<object> { to, rLo, rGo, 0.5 }, full_output: true);
              var _tup_19 = sol[0];
              T = _tup_19.Item1;
              rhoL = _tup_19.Item2;
              rhoG = _tup_19.Item3;
              x = _tup_19.Item4;
              if(sol[2] == 1 && 0 <= x <= 1 && abs(sol[1]["fvec"]).Sum() < 1E-05)
              {
                break;
              }
            }
            if(abs(sol[1]["fvec"]).Sum() > 1E-05)
            {
              throw new RuntimeError(sol[3]);
            }
            liquido = this._Helmholtz(rhoL, T);
            vapor = this._Helmholtz(rhoG, T);
            P = this.R * T * rhoL * rhoG / (rhoL - rhoG) * (liquido["fir"] - vapor["fir"] + Math.Log(rhoL / rhoG)) / 1000;
          }
        }
        else if(this._mode == "su")
        {

          Func<double, Tuple> f = parr =>
          {
            var _tup_1 = parr;
            var rho = _tup_1.Item1;
            var T = _tup_1.Item2;
            var delta = rho / rhoc;
            var tau = Tc / T;
            var ideal = this._phi0(tau, delta);
            var fio = ideal["fio"];
            var fiot = ideal["fiot"];
            var fird = _phird(tau, delta, this._constants);
            var fir = _phir(tau, delta, this._constants);
            var firt = _phirt(tau, delta, this._constants);
            var ho = this.R * T * (1 + tau * (fiot + firt) + delta * fird);
            var so = this.R * (tau * (fiot + firt) - fio - fir);
            var Po = (1 + delta * fird) * this.R * T * rho;
            return Tuple.Create(ho - Po / rho - u, so - s);
          };

          sol = fsolve(f, new List<object> { rhoo, To }, full_output: true);
          var _tup_20 = sol[0];
          rho = _tup_20.Item1;
          T = _tup_20.Item2;
          rhol = this._Liquid_Density(T);
          rhov = this._Vapor_Density(T);
          if(sol[2] != 1 || rhov <= rho <= rhol)
          {

            Func<double, double> f = parr =>
            {
              var _tup_1 = parr;
              var T = _tup_1.Item1;
              var rhol = _tup_1.Item2;
              var rhog = _tup_1.Item3;
              var x = _tup_1.Item4;
              var tau = Tc / T;
              var deltaL = rhol / this.rhoc;
              var deltaG = rhog / this.rhoc;
              var idealL = this._phi0(tau, deltaL);
              var fiot = idealL["fiot"];
              var fioL = idealL["fio"];
              var idealG = this._phi0(tau, deltaG);
              var fioG = idealG["fio"];
              var firL = _phir(tau, deltaL, this._constants);
              var firdL = _phird(tau, deltaL, this._constants);
              var firtL = _phirt(tau, deltaL, this._constants);
              var hoL = this.R * T * (1 + tau * (fiot + firtL) + deltaL * firdL);
              var soL = this.R * (tau * (fiot + firtL) - fioL - firL);
              var firG = _phir(tau, deltaG, this._constants);
              var firdG = _phird(tau, deltaG, this._constants);
              var firtG = _phirt(tau, deltaG, this._constants);
              var hoG = this.R * T * (1 + tau * (fiot + firtG) + deltaG * firdG);
              var soG = this.R * (tau * (fiot + firtG) - fioG - firG);
              var Jl = rhol * (1 + deltaL * firdL);
              var Jv = rhog * (1 + deltaG * firdG);
              var K = firL - firG;
              var Ps = this.R * T * rhol * rhog / (rhol - rhog) * (K + Math.Log(rhol / rhog));
              var vu = hoG - Ps / rhog;
              var lu = hoL - Ps / rhol;
              return (Jl - Jv, Jl * (1 / rhog - 1 / rhol) - Math.Log(rhol / rhog) - K, soL * (1 - x) + soG * x - s, lu * (1 - x) + vu * x - u);
            };

            foreach (var to in new List<int> { To, 300, 400, 500, 600 })
            {
              rLo = this._Liquid_Density(to);
              rGo = this._Vapor_Density(to);
              sol = fsolve(f, new List<object> { to, rLo, rGo, 0.5 }, full_output: true);
              var _tup_21 = sol[0];
              T = _tup_21.Item1;
              rhoL = _tup_21.Item2;
              rhoG = _tup_21.Item3;
              x = _tup_21.Item4;
              if(sol[2] == 1 && 0 <= x <= 1 && abs(sol[1]["fvec"]).Sum() < 1E-05)
              {
                break;
              }
            }
            if(abs(sol[1]["fvec"]).Sum() > 1E-05)
            {
              throw new RuntimeError(sol[3]);
            }
            liquido = this._Helmholtz(rhoL, T);
            vapor = this._Helmholtz(rhoG, T);
            P = this.R * T * rhoL * rhoG / (rhoL - rhoG) * (liquido["fir"] - vapor["fir"] + Math.Log(rhoL / rhoG)) / 1000;
          }
        }
        else if(this._mode == "Trho")
        {
          if(T < this.Tc)
          {
            rhov = this._Vapor_Density(T);
            rhol = this._Liquid_Density(T);
            if(rhol > rho > rhov)
            {
              var _tup_22 = this._saturation(T);
              rhol = _tup_22.Item1;
              rhov = _tup_22.Item2;
              Ps = _tup_22.Item3;
              if(rhol > rho > rhov)
              {
                vapor = this._Helmholtz(rhov, T);
                liquido = this._Helmholtz(rhol, T);
                x = (1 / rho - 1 / rhol) / (1 / rhov - 1 / rhol);
                P = Ps / 1000;
              }
            }
          }
        }
        rho = (float)rho;
        T = (float)T;
        propiedades = this._Helmholtz(rho, T);
        if(T > this.Tc) { x = 1; }
        else if(x == null) { x = 0; }
        if(!P) { P = propiedades["P"] / 1000.0; }
      }
      else if(this._mode == "Tx")
      {
        // Check input T in saturation range
        if(this.Tt > T || this.Tc < T || x > 1 || x < 0)
        {
          throw new NotImplementedException("Incoming out of bound");
        }
        var _tup_23 = this._saturation(T);
        rhol = _tup_23.Item1;
        rhov = _tup_23.Item2;
        Ps = _tup_23.Item3;
        vapor = this._Helmholtz(rhov, T);
        liquido = this._Helmholtz(rhol, T);
        if(x == 0) { propiedades = liquido; }
        else if(x == 1) { propiedades = vapor; }
        P = Ps / 1000.0;
      }
      else if(this._mode == "Px")
      {
        // Check input P in saturation range
        if(this.Pc < P || x > 1 || x < 0)
        {
          throw new NotImplementedException("Incoming out of bound");
        }

        // Iterate over saturation routine to get T
        Func<double, double> f = T =>
        {
          var rhol = this._Liquid_Density(T);
          var rhog = this._Vapor_Density(T);
          var deltaL = rhol / this.rhoc;
          var deltaG = rhog / this.rhoc;
          var tau = Tc / T;
          var firL = _phir(tau, deltaL, this._constants);
          var firG = _phir(tau, deltaG, this._constants);
          var Ps = this.R * T * rhol * rhog / (rhol - rhog) * (firL - firG + Math.Log(deltaL / deltaG));
          return Ps / 1000 - P;
        };

        if(T0) { To = T0; }
        else if(this.name == "water") { To = iapws97._TSat_P(P); }
        else { To = (this.Tc + this.Tt) / 2; }
        T = fsolve(f, To)[0];
        var _tup_24 = this._saturation(T);
        rhol = _tup_24.Item1;
        rhov = _tup_24.Item2;
        Ps = _tup_24.Item3;
        vapor = this._Helmholtz(rhov, T);
        liquido = this._Helmholtz(rhol, T);
        if(x == 0) { propiedades = liquido; }
        else if(x == 1) { propiedades = vapor; }
      }

      Func<double, double> f = rho =>
      {
        var delta = rho / rhoc;
        var tau = Tc / T;
        var fird = _phird(tau, delta, this._constants);
        var Po = (1 + delta * fird) * this.R * T * rho;
        return Po - P * 1000;
      };

      Func<double, double> f = parr =>
      {
        var _tup_1 = parr;
        var T = _tup_1.Item1;
        var rhol = _tup_1.Item2;
        var rhog = _tup_1.Item3;
        var tau = Tc / T;
        var deltaL = rhol / this.rhoc;
        var deltaG = rhog / this.rhoc;
        var idealL = this._phi0(tau, deltaL);
        var fioL = idealL["fio"];
        var fiot = idealL["fiot"];
        var idealG = this._phi0(tau, deltaG);
        var fioG = idealG["fio"];
        var firL = _phir(tau, deltaL, this._constants);
        var firdL = _phird(tau, deltaL, this._constants);
        var firtL = _phirt(tau, deltaL, this._constants);
        var soL = this.R * (tau * (fiot + firtL) - fioL - firL);
        var firG = _phir(tau, deltaG, this._constants);
        var firdG = _phird(tau, deltaG, this._constants);
        var firtG = _phirt(tau, deltaG, this._constants);
        var soG = this.R * (tau * (fiot + firtG) - fioG - firG);
        var Jl = rhol * (1 + deltaL * firdL);
        var Jv = rhog * (1 + deltaG * firdG);
        var K = firL - firG;
        var x = (1.0 / rho - 1 / rhol) / (1 / rhog - 1 / rhol);
        return (Jl - Jv, Jl * (1 / rhog - 1 / rhol) - Math.Log(rhol / rhog) - K, soL * (1 - x) + soG * x - s);
      };

      this.T = T;
      this.Tr = T / this.Tc;
      this.P = P;
      this.Pr = this.P / this.Pc;
      this.x = x;
      if(new List<string> { "Tx", "Px" }.Contains(this._mode) || 0 < x < 1) { region = 4; }
      else { region = 0; }
      this.phase = _utils.getphase(this.Tc, this.Pc, this.T, this.P, this.x, region);
      this.Liquid = _utils._fase();
      this.Gas = _utils._fase();
      if(x == 0)
      {
        // liquid phase
        this.fill(this.Liquid, propiedades);
        this.fill(this, propiedades);
      }
      else if(x == 1)
      {
        // vapor phase
        this.fill(this.Gas, propiedades);
        this.fill(this, propiedades);
      }
      else
      {
        this.fill(this.Liquid, liquido);
        this.fill(this.Gas, vapor);
        this.v = x * this.Gas.v + (1 - x) * this.Liquid.v;
        this.rho = 1.0 / this.v;
        this.h = x * this.Gas.h + (1 - x) * this.Liquid.h;
        this.s = x * this.Gas.s + (1 - x) * this.Liquid.s;
        this.u = x * this.Gas.u + (1 - x) * this.Liquid.u;
        this.a = x * this.Gas.a + (1 - x) * this.Liquid.a;
        this.g = x * this.Gas.g + (1 - x) * this.Liquid.g;
        this.Z = x * this.Gas.Z + (1 - x) * this.Liquid.Z;
        this.f = x * this.Gas.f + (1 - x) * this.Liquid.f;
        this.Z_rho = x * this.Gas.Z_rho + (1 - x) * this.Liquid.Z_rho;
        this.IntP = x * this.Gas.IntP + (1 - x) * this.Liquid.IntP;
      }
      // Calculate special properties useful only for one phase
      if(("Px", "Tx").Contains(this._mode) || x < 1 && this.Tt <= T <= this.Tc)
      {
        this.sigma = this._surface(T);
      }
      else
      {
        this.sigma = null;
      }
      var vir = this._virial(T);
      this.virialB = vir["B"] / this.rhoc;
      this.virialC = vir["C"] / Math.Pow(this.rhoc, 2);
      if(0 < x < 1)
      {
        this.Hvap = vapor["h"] - liquido["h"];
        this.Svap = vapor["s"] - liquido["s"];
      }
      else
      {
        this.Hvap = null;
        this.Svap = null;
      }
      this.invT = -1 / this.T;
      // Ideal properties
      var cp0 = this._prop0(this.rho, this.T);
      this.v0 = this.R * this.T / this.P / 1000;
      this.rho0 = 1.0 / this.v0;
      this.h0 = cp0.h;
      this.u0 = this.h0 - this.P * this.v0;
      this.s0 = cp0.s;
      this.a0 = this.u0 - this.T * this.s0;
      this.g0 = this.h0 - this.T * this.s0;
      this.cp0 = cp0.cp;
      this.cv0 = cp0.cv;
      this.cp0_cv = this.cp0 / this.cv0;
      cp0.v = this.v0;
      this.gamma0 = -this.v0 / this.P / 1000 * this.derivative("P", "v", "s", cp0);
    }

    // Fill phase properties
    public virtual object fill(object fase, object estado)
    {
      fase.rho = estado["rho"];
      fase.v = 1 / fase.rho;
      fase.h = estado["h"];
      fase.s = estado["s"];
      fase.u = fase.h - this.P * 1000 * fase.v;
      fase.a = fase.u - this.T * fase.s;
      fase.g = fase.h - this.T * fase.s;
      fase.Z = this.P * fase.v / this.T / this.R * 1000.0;
      fase.fi = Math.Exp(estado["fir"] + estado["delta"] * estado["fird"] - Math.Log(1 + estado["delta"] * estado["fird"]));
      fase.f = fase.fi * this.P;
      fase.cv = estado["cv"];
      fase.rhoM = fase.rho / this.M;
      fase.hM = fase.h * this.M;
      fase.sM = fase.s * this.M;
      fase.uM = fase.u * this.M;
      fase.aM = fase.a * this.M;
      fase.gM = fase.g * this.M;
      fase.alfap = estado["alfap"];
      fase.betap = estado["betap"];
      fase.cp = this.derivative("h", "T", "P", fase);
      fase.cp_cv = fase.cp / fase.cv;
      fase.w = Math.Pow(this.derivative("P", "rho", "s", fase) * 1000, 0.5);
      fase.cvM = fase.cv * this.M;
      fase.cpM = fase.cp * this.M;
      fase.joule = this.derivative("T", "P", "h", fase) * 1000.0;
      fase.Gruneisen = fase.v / fase.cv * this.derivative("P", "T", "v", fase);
      fase.alfav = this.derivative("v", "T", "P", fase) / fase.v;
      fase.kappa = -this.derivative("v", "P", "T", fase) / fase.v * 1000.0;
      fase.betas = this.derivative("T", "P", "s", fase);
      fase.gamma = -fase.v / this.P * this.derivative("P", "v", "s", fase) * 0.001;
      fase.kt = -fase.v / this.P * this.derivative("P", "v", "T", fase) * 0.001;
      fase.ks = -this.derivative("v", "P", "s", fase) / fase.v * 1000.0;
      fase.Kt = -fase.v * this.derivative("P", "v", "s", fase) * 0.001;
      fase.Ks = -fase.v * this.derivative("P", "v", "T", fase) * 0.001;
      fase.dhdT_rho = this.derivative("h", "T", "rho", fase);
      fase.dhdT_P = this.derivative("h", "T", "P", fase);
      fase.dhdP_T = this.derivative("h", "P", "T", fase) * 1000.0;
      fase.dhdP_rho = this.derivative("h", "P", "rho", fase) * 1000.0;
      fase.dhdrho_T = this.derivative("h", "rho", "T", fase);
      fase.dhdrho_P = this.derivative("h", "rho", "P", fase);
      fase.dpdT_rho = this.derivative("P", "T", "rho", fase) * 0.001;
      fase.dpdrho_T = this.derivative("P", "rho", "T", fase) * 0.001;
      fase.drhodP_T = this.derivative("rho", "P", "T", fase) * 1000.0;
      fase.drhodT_P = this.derivative("rho", "T", "P", fase);
      fase.Z_rho = (fase.Z - 1) / fase.rho;
      fase.IntP = this.T * this.derivative("P", "T", "rho", fase) * 0.001 - this.P;
      fase.hInput = fase.v * this.derivative("h", "v", "P", fase);
      fase.mu = this._visco(fase.rho, this.T, fase);
      fase.k = this._thermo(fase.rho, this.T, fase);
      fase.nu = fase.mu / fase.rho;
      fase.alfa = fase.k / 1000 / fase.rho / fase.cp;
      fase.Prandt = fase.mu * fase.cp * 1000 / fase.k;
      if(this.name == "water")
      {
        try
        {
          fase.epsilon = _iapws._Dielectric(fase.rho, this.T);
        } catch (NotImplementedError)
        {
          fase.epsilon = null;
        }
        try
        {
          fase.n = _iapws._Refractive(fase.rho, this.T, this.kwargs["l"]);
        } catch (NotImplementedError)
        {
          fase.n = null;
        }
      }
      else
      {
        fase.epsilon = null;
        fase.n = null;
      }
    }

    // Wrapper derivative for custom derived properties
    //     where x, y, z can be: P, T, v, rho, u, h, s, g, a
    public virtual object derivative(object z, object x, object y, object fase)
    {
      return _utils.deriv_H(this, z, x, y, fase);
    }

    // Saturation calculation for two phase search
    public virtual object _saturation(object T)
    {
      object Ps;
      var rhoc = this._constants.get("rhoref", this.rhoc);
      var Tc = this._constants.get("Tref", this.Tc);
      if(T > Tc) { T = Tc; }
      var tau = Tc / T;
      var rhoLo = this._Liquid_Density(T);
      var rhoGo = this._Vapor_Density(T);

      Func<double, Tuple> f = parr =>
      {
        var _tup_1 = parr;
        var rhol = _tup_1.Item1;
        var rhog = _tup_1.Item2;
        var deltaL = rhol / rhoc;
        var deltaG = rhog / rhoc;
        var phirL = _phir(tau, deltaL, this._constants);
        var phirG = _phir(tau, deltaG, this._constants);
        var phirdL = _phird(tau, deltaL, this._constants);
        var phirdG = _phird(tau, deltaG, this._constants);
        var Jl = deltaL * (1 + deltaL * phirdL);
        var Jv = deltaG * (1 + deltaG * phirdG);
        var Kl = deltaL * phirdL + phirL + Math.Log(deltaL);
        var Kv = deltaG * phirdG + phirG + Math.Log(deltaG);
        return Tuple.Create(Kv - Kl, Jv - Jl);
      };

      var _tup_1 = fsolve(f, new List<object> { rhoLo, rhoGo });

      var rhoL = _tup_1.Item1;
      var rhoG = _tup_1.Item2;
      if(rhoL == rhoG)
      {
        Ps = this.Pc;
      }
      else
      {
        var deltaL = rhoL / this.rhoc;
        var deltaG = rhoG / this.rhoc;
        var firL = _phir(tau, deltaL, this._constants);
        var firG = _phir(tau, deltaG, this._constants);
        Ps = this.R * T * rhoL * rhoG / (rhoL - rhoG) * (firL - firG + Math.Log(deltaL / deltaG));
      }
      return Tuple.Create(rhoL, rhoG, Ps);
    }

    // Calculated properties from helmholtz free energy and derivatives
    //
    //     Parameters
    //     ----------
    //     rho : float
    //       Density, [kg/m³]
    //     T : float
    //       Temperature, [K]
    //
    //     Returns
    //     -------
    //     prop : dict
    //       Dictionary with calculated properties:
    //         * fir: [-]
    //         * fird: ∂fir/∂δ|τ
    //         * firdd: ∂²fir/∂δ²|τ
    //         * delta: Reducen density rho/rhoc, [-]
    //         * P: Pressure, [kPa]
    //         * h: Enthalpy, [kJ/kg]
    //         * s: Entropy, [kJ/kgK]
    //         * cv: Isochoric specific heat, [kJ/kgK]
    //         * alfav: Thermal expansion coefficient, [1/K]
    //         * betap: Isothermal compressibility, [1/kPa]
    //
    //     References
    //     ----------
    //     IAPWS, Revised Release on the IAPWS Formulation 1995 for the
    //     Thermodynamic Properties of Ordinary Water Substance for General and
    //     Scientific Use, September 2016, Table 3
    //     http://www.iapws.org/relguide/IAPWS-95.html
    //

    public virtual Dictionary<string, double> _Helmholtz(object rho, object T)
    {
      if(rho is ndarray) { rho = rho[0]; }
      if(T is ndarray) { T = T[0]; }
      if(rho < 0) { rho = 1E-20; }
      if(T < 50) { T = 50; }

      var rhoc = this._constants.get("rhoref", this.rhoc);
      var Tc = this._constants.get("Tref", this.Tc);
      var delta = rho / rhoc;
      var tau = Tc / T;
      var ideal = this._phi0(tau, delta);
      var fio = ideal["fio"];
      var fiot = ideal["fiot"];
      var fiott = ideal["fiott"];
      var res = this._phir(tau, delta);
      var fir = res["fir"];
      var firt = res["firt"];
      var firtt = res["firtt"];
      var fird = res["fird"];
      var firdd = res["firdd"];
      var firdt = res["firdt"];
      var propiedades = new Dictionary<string, double> { };

      propiedades["fir"] = fir;
      propiedades["fird"] = fird;
      propiedades["firdd"] = firdd;
      propiedades["delta"] = delta;
      propiedades["rho"] = rho;
      propiedades["P"] = (1 + delta * fird) * this.R * T * rho;
      propiedades["h"] = this.R * T * (1 + tau * (fiot + firt) + delta * fird);
      propiedades["s"] = this.R * (tau * (fiot + firt) - fio - fir);
      propiedades["cv"] = -this.R * Math.Pow(tau, 2) * (fiott + firtt);
      propiedades["alfap"] = (1 - delta * tau * firdt / (1 + delta * fird)) / T;
      propiedades["betap"] = rho * (1 + (delta * fird + Math.Pow(delta, 2) * firdd) / (1 + delta * fird));

      return propiedades;
    }

    // Ideal gas properties
    public virtual _utils._fase _prop0(object rho, object T)
    {
      var rhoc = this._constants.get("rhoref", this.rhoc);
      var Tc = this._constants.get("Tref", this.Tc);
      var delta = rho / rhoc;
      var tau = Tc / T;
      var ideal = this._phi0(tau, delta);
      var fio = ideal["fio"];
      var fiot = ideal["fiot"];
      var fiott = ideal["fiott"];

      var propiedades = _utils._fase();

      propiedades.h = this.R * T * (1 + tau * fiot);
      propiedades.s = this.R * (tau * fiot - fio);
      propiedades.cv = -this.R * Math.Pow(tau, 2) * fiott;
      propiedades.cp = this.R * (-Math.Pow(tau, 2) * fiott + 1);
      propiedades.alfap = 1 / T;
      propiedades.betap = rho;

      return propiedades;
    }

    // Ideal gas Helmholtz free energy and derivatives
    //
    //     Parameters
    //     ----------
    //     tau : float
    //       Inverse reduced temperature Tc/T, [-]
    //     delta : float
    //       Reduced density rho/rhoc, [-]
    //
    //     Returns
    //     -------
    //     prop : dictionary with ideal adimensional helmholtz energy and deriv
    //       fio, [-]
    //       fiot: ∂fio/∂τ|δ
    //       fiod: ∂fio/∂δ|τ
    //       fiott: ∂²fio/∂τ²|δ
    //       fiodt: ∂²fio/∂τ∂δ
    //       fiodd: ∂²fio/∂δ²|τ
    //
    //     References
    //     ----------
    //     IAPWS, Revised Release on the IAPWS Formulation 1995 for the
    //     Thermodynamic Properties of Ordinary Water Substance for General and
    //     Scientific Use, September 2016, Table 4
    //     http://www.iapws.org/relguide/IAPWS-95.html
    //

    public virtual Dictionary<string, double> _phi0(object tau, object delta)
    {
      object t;
      object n;
      var Fi0 = this.Fi0;
      var fio = Fi0["ao_log"][0] * Math.Log(delta) + Fi0["ao_log"][1] * Math.Log(tau);
      var fiot = +Fi0["ao_log"][1] / tau;
      var fiott = -Fi0["ao_log"][1] / Math.Pow(tau, 2);
      var fiod = 1 / delta;
      var fiodd = -1 / Math.Pow(delta, 2);
      var fiodt = 0;

      foreach (var _tup_1 in zip(Fi0["ao_pow"], Fi0["pow"]))
      {
        n = _tup_1.Item1;
        t = _tup_1.Item2;
        fio += n * Math.Pow(tau, t);
        if(t != 0)
        {
          fiot += t * n * Math.Pow(tau, t - 1);
        }
        if(!new List<int>
        {
          0,
          1
        }.Contains(t))
        {
          fiott += n * t * (t - 1) * Math.Pow(tau, t - 2);
        }
      }

      foreach (var _tup_2 in zip(Fi0["ao_exp"], Fi0["titao"]))
      {
        n = _tup_2.Item1;
        t = _tup_2.Item2;
        fio += n * Math.Log(1 - Math.Exp(-tau * t));
        fiot += n * t * (Math.Pow(1 - Math.Exp(-t * tau), -1) - 1);
        fiott -= n * Math.Pow(t, 2) * Math.Exp(-t * tau) * Math.Pow(1 - Math.Exp(-t * tau), -2);
      }

      // Extension to especial terms of air
      if(Fi0.Contains("ao_exp2"))
      {
        foreach (var _tup_3 in zip(Fi0["ao_exp2"], Fi0["titao2"], Fi0["sum2"]))
        {
          n = _tup_3.Item1;
          var g = _tup_3.Item2;
          var C = _tup_3.Item3;
          fio += n * Math.Log(C + Math.Exp(g * tau));
          fiot += n * g / (C * Math.Exp(-g * tau) + 1);
          fiott += C * n * Math.Pow(g, 2) * Math.Exp(-g * tau) / Math.Pow(C * Math.Exp(-g * tau) + 1, 2);
        }
      }

      var prop = new Dictionary<string, double> { };

      prop["fio"] = fio;
      prop["fiot"] = fiot;
      prop["fiott"] = fiott;
      prop["fiod"] = fiod;
      prop["fiodd"] = fiodd;
      prop["fiodt"] = fiodt;

      return prop;
    }

    // Residual contribution to the free Helmholtz energy
    //
    //     Parameters
    //     ----------
    //     tau : float
    //       Inverse reduced temperature Tc/T, [-]
    //     delta : float
    //       Reduced density rho/rhoc, [-]
    //
    //     Returns
    //     -------
    //     prop : dict
    //       Dictionary with residual adimensional helmholtz energy and deriv:
    //       * fir
    //       * firt: ∂fir/∂τ|δ,x
    //       * fird: ∂fir/∂δ|τ,x
    //       * firtt: ∂²fir/∂τ²|δ,x
    //       * firdt: ∂²fir/∂τ∂δ|x
    //       * firdd: ∂²fir/∂δ²|τ,x
    //
    //     References
    //     ----------
    //     IAPWS, Revised Release on the IAPWS Formulation 1995 for the
    //     Thermodynamic Properties of Ordinary Water Substance for General and
    //     Scientific Use, September 2016, Table 5
    //     http://www.iapws.org/relguide/IAPWS-95.html
    //

    public virtual Dictionary<string, double> _phir(object tau, object delta)
    {
      object Deltadd;
      object b;
      object a;
      object g;
      object t;
      object d;
      object n;

      var fir = 0;
      // Polinomial terms
      var nr1 = this._constants.get("nr1", new List<object>());
      var d1 = this._constants.get("d1", new List<object>());
      var t1 = this._constants.get("t1", new List<object>());

      foreach (var _tup_1 in zip(nr1, d1, t1))
      {
        n = _tup_1.Item1;
        d = _tup_1.Item2;
        t = _tup_1.Item3;
        fir += n * Math.Pow(delta, d) * Math.Pow(tau, t);
        fird += n * d * Math.Pow(delta, d - 1) * Math.Pow(tau, t);
        firdd += n * d * (d - 1) * Math.Pow(delta, d - 2) * Math.Pow(tau, t);
        firt += n * t * Math.Pow(delta, d) * Math.Pow(tau, t - 1);
        firtt += n * t * (t - 1) * Math.Pow(delta, d) * Math.Pow(tau, t - 2);
        firdt += n * t * d * Math.Pow(delta, d - 1) * Math.Pow(tau, t - 1);
      }

      // Exponential terms
      var nr2 = this._constants.get("nr2", new List<object>());
      var d2 = this._constants.get("d2", new List<object>());
      var g2 = this._constants.get("gamma2", new List<object>());
      var t2 = this._constants.get("t2", new List<object>());
      var c2 = this._constants.get("c2", new List<object>());

      foreach (var _tup_2 in zip(nr2, d2, g2, t2, c2))
      {
        n = _tup_2.Item1;
        d = _tup_2.Item2;
        g = _tup_2.Item3;
        t = _tup_2.Item4;
        var c = _tup_2.Item5;
        fir += n * Math.Pow(delta, d) * Math.Pow(tau, t) * Math.Exp(-g * Math.Pow(delta, c));
        fird += n * Math.Exp(-g * Math.Pow(delta, c)) * Math.Pow(delta, d - 1) * Math.Pow(tau, t) * (d - g * c * Math.Pow(delta, c));
        firdd += n * Math.Exp(-g * Math.Pow(delta, c)) * Math.Pow(delta, d - 2) * Math.Pow(tau, t) * ((d - g * c * Math.Pow(delta, c)) * (d - 1 - g * c * Math.Pow(delta, c)) - Math.Pow(g, 2) * Math.Pow(c, 2) * Math.Pow(delta, c));
        firt += n * t * Math.Pow(delta, d) * Math.Pow(tau, t - 1) * Math.Exp(-g * Math.Pow(delta, c));
        firtt += n * t * (t - 1) * Math.Pow(delta, d) * Math.Pow(tau, t - 2) * Math.Exp(-g * Math.Pow(delta, c));
        firdt += n * t * Math.Pow(delta, d - 1) * Math.Pow(tau, t - 1) * (d - g * c * Math.Pow(delta, c)) * Math.Exp(-g * Math.Pow(delta, c));
      }

      // Gaussian terms
      var nr3 = this._constants.get("nr3", new List<object>());
      var d3 = this._constants.get("d3", new List<object>());
      var t3 = this._constants.get("t3", new List<object>());
      var a3 = this._constants.get("alfa3", new List<object>());
      var e3 = this._constants.get("epsilon3", new List<object>());
      var b3 = this._constants.get("beta3", new List<object>());
      var g3 = this._constants.get("gamma3", new List<object>());

      foreach (var _tup_3 in zip(nr3, d3, t3, a3, e3, b3, g3))
      {
        n = _tup_3.Item1;
        d = _tup_3.Item2;
        t = _tup_3.Item3;
        a = _tup_3.Item4;
        var e = _tup_3.Item5;
        b = _tup_3.Item6;
        g = _tup_3.Item7;
        fir += n * Math.Pow(delta, d) * Math.Pow(tau, t) * Math.Exp(-a * Math.Pow(delta - e, 2) - b * Math.Pow(tau - g, 2));
        fird += n * Math.Pow(delta, d) * Math.Pow(tau, t) * Math.Exp(-a * Math.Pow(delta - e, 2) - b * Math.Pow(tau - g, 2)) * (d / delta - 2 * a * (delta - e));
        firdd += n * Math.Pow(tau, t) * Math.Exp(-a * Math.Pow(delta - e, 2) - b * Math.Pow(tau - g, 2)) * (-2 * a * Math.Pow(delta, d) + 4 * Math.Pow(a, 2) * Math.Pow(delta, d) * Math.Pow(delta - e, 2) - 4 * d * a * Math.Pow(delta, d - 1) * (delta - e) + d * (d - 1) * Math.Pow(delta, d - 2));
        firt += n * Math.Pow(delta, d) * Math.Pow(tau, t) * Math.Exp(-a * Math.Pow(delta - e, 2) - b * Math.Pow(tau - g, 2)) * (t / tau - 2 * b * (tau - g));
        firtt += n * Math.Pow(delta, d) * Math.Pow(tau, t) * Math.Exp(-a * Math.Pow(delta - e, 2) - b * Math.Pow(tau - g, 2)) * (Math.Pow(t / tau - 2 * b * (tau - g), 2) - t / Math.Pow(tau, 2) - 2 * b);
        firdt += n * Math.Pow(delta, d) * Math.Pow(tau, t) * Math.Exp(-a * Math.Pow(delta - e, 2) - b * Math.Pow(tau - g, 2)) * (t / tau - 2 * b * (tau - g)) * (d / delta - 2 * a * (delta - e));
      }

      // Non analitic terms
      var nr4 = this._constants.get("nr4", new List<object>());
      var a4 = this._constants.get("a4", new List<object>());
      var b4 = this._constants.get("b4", new List<object>());
      var Ai = this._constants.get("A", new List<object>());
      var Bi = this._constants.get("B", new List<object>());
      var Ci = this._constants.get("C", new List<object>());
      var Di = this._constants.get("D", new List<object>());
      var bt4 = this._constants.get("beta4", new List<object>());

      foreach (var _tup_4 in zip(nr4, a4, b4, Ai, Bi, Ci, Di, bt4))
      {
        n = _tup_4.Item1;
        a = _tup_4.Item2;
        b = _tup_4.Item3;
        var A = _tup_4.Item4;
        var B = _tup_4.Item5;
        var C = _tup_4.Item6;
        var D = _tup_4.Item7;
        var bt = _tup_4.Item8;
        var Tita = 1 - tau + A * Math.Pow(Math.Pow(delta - 1, 2), 0.5 / bt);
        var F = Math.Exp(-C * Math.Pow(delta - 1, 2) - D * Math.Pow(tau - 1, 2));
        var Fd = -2 * C * F * (delta - 1);
        var Fdd = 2 * C * F * (2 * C * Math.Pow(delta - 1, 2) - 1);
        var Ft = -2 * D * F * (tau - 1);
        var Ftt = 2 * D * F * (2 * D * Math.Pow(tau - 1, 2) - 1);
        var Fdt = 4 * C * D * F * (delta - 1) * (tau - 1);
        var Delta = Math.Pow(Tita, 2) + B * Math.Pow(Math.Pow(delta - 1, 2), a);
        var Deltad = (delta - 1) * (A * Tita * 2 / bt * Math.Pow(Math.Pow(delta - 1, 2), 0.5 / bt - 1) + 2 * B * a * Math.Pow(Math.Pow(delta - 1, 2), a - 1));
        if(delta == 1)
        {
          Deltadd = 0;
        }
        else
        {
          Deltadd = Deltad / (delta - 1) + Math.Pow(delta - 1, 2) * (4 * B * a * (a - 1) * Math.Pow(Math.Pow(delta - 1, 2), a - 2) + 2 * Math.Pow(A, 2) / Math.Pow(bt, 2) * Math.Pow(Math.Pow(Math.Pow(delta - 1, 2), 0.5 / bt - 1), 2) + A * Tita * 4 / bt * (0.5 / bt - 1) * Math.Pow(Math.Pow(delta - 1, 2), 0.5 / bt - 2));
        }
        var DeltaBd = b * Math.Pow(Delta, b - 1) * Deltad;
        var DeltaBdd = b * (Math.Pow(Delta, b - 1) * Deltadd + (b - 1) * Math.Pow(Delta, b - 2) * Math.Pow(Deltad, 2));
        var DeltaBt = -2 * Tita * b * Math.Pow(Delta, b - 1);
        var DeltaBtt = 2 * b * Math.Pow(Delta, b - 1) + 4 * Math.Pow(Tita, 2) * b * (b - 1) * Math.Pow(Delta, b - 2);
        var DeltaBdt = -A * b * 2 / bt * Math.Pow(Delta, b - 1) * (delta - 1) * Math.Pow(Math.Pow(delta - 1, 2), 0.5 / bt - 1) - 2 * Tita * b * (b - 1) * Math.Pow(Delta, b - 2) * Deltad;
        fir += n * Math.Pow(Delta, b) * delta * F;
        fird += n * (Math.Pow(Delta, b) * (F + delta * Fd) + DeltaBd * delta * F);
        firdd += n * (Math.Pow(Delta, b) * (2 * Fd + delta * Fdd) + 2 * DeltaBd * (F + delta * Fd) + DeltaBdd * delta * F);
        firt += n * delta * (DeltaBt * F + Math.Pow(Delta, b) * Ft);
        firtt += n * delta * (DeltaBtt * F + 2 * DeltaBt * Ft + Math.Pow(Delta, b) * Ftt);
        firdt += n * (Math.Pow(Delta, b) * (Ft + delta * Fdt) + delta * DeltaBd * Ft + DeltaBt * (F + delta * Fd) + DeltaBdt * delta * F);
      }

      var prop = new Dictionary<string, double> { };
      prop["fir"] = fir;
      prop["firt"] = firt;
      prop["firtt"] = firtt;
      prop["fird"] = fird;
      prop["firdd"] = firdd;
      prop["firdt"] = firdt;

      return prop;
    }

    // Virial coefficient
    //
    //     Parameters
    //     ----------
    //     T : float
    //       Temperature [K]
    //
    //     Returns
    //     -------
    //     prop : dict
    //       Dictionary with residual adimensional helmholtz energy:
    //         * B: ∂fir/∂δ|δ->0
    //         * C: ∂²fir/∂δ²|δ->0
    //

    public virtual Dictionary<string, double> _virial(object T)
    {
      object b;
      object a;
      object g;
      object t;
      object d;
      object n;
      var Tc = this._constants.get("Tref", this.Tc);
      var tau = Tc / T;
      var B = 0;
      var delta = 1E-200;
      // Polinomial terms
      var nr1 = this._constants.get("nr1", new List<object>());
      var d1 = this._constants.get("d1", new List<object>());
      var t1 = this._constants.get("t1", new List<object>());

      foreach (var _tup_1 in zip(nr1, d1, t1))
      {
        n = _tup_1.Item1;
        d = _tup_1.Item2;
        t = _tup_1.Item3;
        B += n * d * Math.Pow(delta, d - 1) * Math.Pow(tau, t);
        C += n * d * (d - 1) * Math.Pow(delta, d - 2) * Math.Pow(tau, t);
      }

      // Exponential terms
      var nr2 = this._constants.get("nr2", new List<object>());
      var d2 = this._constants.get("d2", new List<object>());
      var g2 = this._constants.get("gamma2", new List<object>());
      var t2 = this._constants.get("t2", new List<object>());
      var c2 = this._constants.get("c2", new List<object>());

      foreach (var _tup_2 in zip(nr2, d2, g2, t2, c2))
      {
        n = _tup_2.Item1;
        d = _tup_2.Item2;
        g = _tup_2.Item3;
        t = _tup_2.Item4;
        var c = _tup_2.Item5;
        B += n * Math.Exp(-g * Math.Pow(delta, c)) * Math.Pow(delta, d - 1) * Math.Pow(tau, t) * (d - g * c * Math.Pow(delta, c));
        C += n * Math.Exp(-g * Math.Pow(delta, c)) * (Math.Pow(delta, d - 2) * Math.Pow(tau, t) * ((d - g * c * Math.Pow(delta, c)) * (d - 1 - g * c * Math.Pow(delta, c)) - Math.Pow(g, 2) * Math.Pow(c, 2) * Math.Pow(delta, c)));
      }

      // Gaussian terms
      var nr3 = this._constants.get("nr3", new List<object>());
      var d3 = this._constants.get("d3", new List<object>());
      var t3 = this._constants.get("t3", new List<object>());
      var a3 = this._constants.get("alfa3", new List<object>());
      var e3 = this._constants.get("epsilon3", new List<object>());
      var b3 = this._constants.get("beta3", new List<object>());
      var g3 = this._constants.get("gamma3", new List<object>());

      foreach (var _tup_3 in zip(nr3, d3, t3, a3, e3, b3, g3))
      {
        n = _tup_3.Item1;
        d = _tup_3.Item2;
        t = _tup_3.Item3;
        a = _tup_3.Item4;
        var e = _tup_3.Item5;
        b = _tup_3.Item6;
        g = _tup_3.Item7;
        B += n * Math.Pow(delta, d) * Math.Pow(tau, t) * Math.Exp(-a * Math.Pow(delta - e, 2) - b * Math.Pow(tau - g, 2)) * (d / delta - 2 * a * (delta - e));
        C += n * Math.Pow(tau, t) * Math.Exp(-a * Math.Pow(delta - e, 2) - b * Math.Pow(tau - g, 2)) * (-2 * a * Math.Pow(delta, d) + 4 * Math.Pow(a, 2) * Math.Pow(delta, d) * Math.Pow(delta - e, 2) - 4 * d * a * Math.Pow(delta, 2) * (delta - e) + d * 2 * delta);
      }

      // Non analitic terms
      var nr4 = this._constants.get("nr4", new List<object>());
      var a4 = this._constants.get("a4", new List<object>());
      var b4 = this._constants.get("b4", new List<object>());
      var Ai = this._constants.get("A", new List<object>());
      var Bi = this._constants.get("B", new List<object>());
      var Ci = this._constants.get("C", new List<object>());
      var Di = this._constants.get("D", new List<object>());
      var bt4 = this._constants.get("beta4", new List<object>());

      foreach (var _tup_4 in zip(nr4, a4, b4, Ai, Bi, Ci, Di, bt4))
      {
        n = _tup_4.Item1;
        a = _tup_4.Item2;
        b = _tup_4.Item3;
        var A = _tup_4.Item4;
        var B_ = _tup_4.Item5;
        var C_ = _tup_4.Item6;
        var D = _tup_4.Item7;
        var bt = _tup_4.Item8;
        var Tita = 1 - tau + A * Math.Pow(Math.Pow(delta - 1, 2), 0.5 / bt);
        var Delta = Math.Pow(Tita, 2) + B_ * Math.Pow(Math.Pow(delta - 1, 2), a);
        var Deltad = (delta - 1) * (A * Tita * 2 / bt * Math.Pow(Math.Pow(delta - 1, 2), 0.5 / bt - 1) + 2 * B_ * a * Math.Pow(Math.Pow(delta - 1, 2), a - 1));
        var Deltadd = Deltad / (delta - 1) + Math.Pow(delta - 1, 2) * (4 * B_ * a * (a - 1) * Math.Pow(Math.Pow(delta - 1, 2), a - 2) + 2 * Math.Pow(A, 2) / Math.Pow(bt, 2) * Math.Pow(Math.Pow(Math.Pow(delta - 1, 2), 0.5 / bt - 1), 2) + A * Tita * 4 / bt * (0.5 / bt - 1) * Math.Pow(Math.Pow(delta - 1, 2), 0.5 / bt - 2));
        var DeltaBd = b * Math.Pow(Delta, b - 1) * Deltad;
        var DeltaBdd = b * (Math.Pow(Delta, b - 1) * Deltadd + (b - 1) * Math.Pow(Delta, b - 2) * Math.Pow(Deltad, 2));
        var F = Math.Exp(-C_ * Math.Pow(delta - 1, 2) - D * Math.Pow(tau - 1, 2));
        var Fd = -2 * C_ * F * (delta - 1);
        var Fdd = 2 * C_ * F * (2 * C_ * Math.Pow(delta - 1, 2) - 1);
        B += n * (Math.Pow(Delta, b) * (F + delta * Fd) + DeltaBd * delta * F);
        C += n * (Math.Pow(Delta, b) * (2 * Fd + delta * Fdd) + 2 * DeltaBd * (F + delta * Fd) + DeltaBdd * delta * F);
      }

      var prop = new Dictionary<string, double> { };
      prop["B"] = B;
      prop["C"] = C;

      return prop;
    }

    // Calcule the dimensional form or Helmholtz free energy derivatives
    //
    //     Parameters
    //     ----------
    //     rho : float
    //       Density, [kg/m³]
    //     T : float
    //       Temperature, [K]
    //
    //     Returns
    //     -------
    //     prop : dict
    //       Dictionary with residual helmholtz energy and derivatives:
    //
    //         * fir, [kJ/kg]
    //         * firt: ∂fir/∂T|ρ, [kJ/kgK]
    //         * fird: ∂fir/∂ρ|T, [kJ/m³kg²]
    //         * firtt: ∂²fir/∂T²|ρ, [kJ/kgK²]
    //         * firdt: ∂²fir/∂T∂ρ, [kJ/m³kg²K]
    //         * firdd: ∂²fir/∂ρ²|T, [kJ/m⁶kg]
    //
    //     References
    //     ----------
    //     IAPWS, Guideline on an Equation of State for Humid Air in Contact with
    //     Seawater and Ice, Consistent with the IAPWS Formulation 2008 for the
    //     Thermodynamic Properties of Seawater, Table 7,
    //     http://www.iapws.org/relguide/SeaAir.html
    //

    public virtual Dictionary<string, double> _derivDimensional(object rho, object T)
    {
      object prop;
      if(!rho)
      {
        prop = new Dictionary<string, double> { };
        prop["fir"] = 0;
        prop["firt"] = 0;
        prop["fird"] = 0;
        prop["firtt"] = 0;
        prop["firdt"] = 0;
        prop["firdd"] = 0;
        return prop;
      }

      var R = this._constants.get("R") / this._constants.get("M", this.M);
      var rhoc = this._constants.get("rhoref", this.rhoc);
      var Tc = this._constants.get("Tref", this.Tc);
      var delta = rho / rhoc;
      var tau = Tc / T;
      var ideal = this._phi0(tau, delta);
      var fio = ideal["fio"];
      var fiot = ideal["fiot"];
      var fiott = ideal["fiott"];
      var fiod = ideal["fiod"];
      var fiodd = ideal["fiodd"];
      var res = this._phir(tau, delta);
      var fir = res["fir"];
      var firt = res["firt"];
      var firtt = res["firtt"];
      var fird = res["fird"];
      var firdd = res["firdd"];
      var firdt = res["firdt"];

      prop = new Dictionary<string, double> { };
      prop["fir"] = R * T * (fio + fir);
      prop["firt"] = R * (fio + fir - (fiot + firt) * tau);
      prop["fird"] = R * T / rhoc * (fiod + fird);
      prop["firtt"] = R * Math.Pow(tau, 2) / T * (fiott + firtt);
      prop["firdt"] = R / rhoc * (fiod + fird - firdt * tau);
      prop["firdd"] = R * T / Math.Pow(rhoc, 2) * (fiodd + firdd);

      return prop;
    }

    // Generic equation for the surface tension
    //
    //     Parameters
    //     ----------
    //     T : float
    //       Temperature, [K]
    //
    //     Returns
    //     -------
    //     σ : float
    //       Surface tension, [N/m]
    //
    //     Notes
    //     -----
    //     Need a _surf dict in the derived class with the parameters keys:
    //       sigma: coefficient
    //       exp: exponent
    //

    public virtual double _surface(double T)
    {
      var tau = 1 - T / this.Tc;
      var sigma = 0;
      foreach (var _tup_1 in zip(this._surf["sigma"], this._surf["exp"]))
      {
        var n = _tup_1.Item1;
        var t = _tup_1.Item2;
        sigma += n * Math.Pow(tau, t);
      }
      return sigma;
    }

    // Auxiliary equation for the vapour pressure
    //
    //     Parameters
    //     ----------
    //     T : float
    //       Temperature, [K]
    //
    //     Returns
    //     -------
    //     Pv : float
    //       Vapour pressure, [Pa]
    //
    //     References
    //     ----------
    //     IAPWS, Revised Supplementary Release on Saturation Properties of
    //     Ordinary Water Substance September 1992,
    //     http://www.iapws.org/relguide/Supp-sat.html, Eq.1
    //

    public static double _Vapor_Pressure(Dictionary<string, double> cls, double T)
    {
      var Tita = 1 - T / cls.Tc;
      var suma = 0;
      foreach (var _tup_1 in zip(cls._Pv["ao"], cls._Pv["exp"]))
      {
        var n = _tup_1.Item1;
        var x = _tup_1.Item2;
        suma += n * Math.Pow(Tita, x);
      }
      var Pr = Math.Exp(cls.Tc / T * suma);
      var Pv = Pr * cls.Pc;
      return Pv;
    }

    // Auxiliary equation for the density of saturated liquid
    //
    //     Parameters
    //     ----------
    //     T : float
    //       Temperature, [K]
    //
    //     Returns
    //     -------
    //     rho : float
    //       Saturated liquid density, [kg/m³]
    //
    //     References
    //     ----------
    //     IAPWS, Revised Supplementary Release on Saturation Properties of
    //     Ordinary Water Substance September 1992,
    //     http://www.iapws.org/relguide/Supp-sat.html, Eq.2
    //

    public static double _Liquid_Density(Dictionary<string, double> cls, double T)
    {
      var eq = cls._rhoL["eq"];
      var Tita = 1 - T / cls.Tc;
      if(eq == 2)
      {
        Tita = Math.Pow(Tita, 1.0 / 3);
      }
      var suma = 0;
      foreach (var _tup_1 in zip(cls._rhoL["ao"], cls._rhoL["exp"]))
      {
        var n = _tup_1.Item1;
        var x = _tup_1.Item2;
        suma += n * Math.Pow(Tita, x);
      }
      var Pr = suma + 1;
      var rho = Pr * cls.rhoc;
      return rho;
    }

    // Auxiliary equation for the density of saturated vapor
    //
    //     Parameters
    //     ----------
    //     T : float
    //       Temperature, [K]
    //
    //     Returns
    //     -------
    //     rho : float
    //       Saturated vapor density, [kg/m³]
    //
    //     References
    //     ----------
    //     IAPWS, Revised Supplementary Release on Saturation Properties of
    //     Ordinary Water Substance September 1992,
    //     http://www.iapws.org/relguide/Supp-sat.html, Eq.3
    //

    public static double _Vapor_Density(Dictionary<string, double> cls, double T)
    {
      var eq = cls._rhoG["eq"];
      var Tita = 1 - T / cls.Tc;
      if(eq == 4)
      {
        Tita = Math.Pow(Tita, 1.0 / 3);
      }
      var suma = 0;
      foreach (var _tup_1 in zip(cls._rhoG["ao"], cls._rhoG["exp"]))
      {
        var n = _tup_1.Item1;
        var x = _tup_1.Item2;
        suma += n * Math.Pow(Tita, x);
      }
      var Pr = Math.Exp(suma);
      var rho = Pr * cls.rhoc;
      return rho;
    }

    // Auxiliary equation for the dP/dT along saturation line
    //
    //     Parameters
    //     ----------
    //     T : float
    //       Temperature, [K]
    //
    //     Returns
    //     -------
    //     dPdT : float
    //       dPdT, [MPa/K]
    //
    //     References
    //     ----------
    //     IAPWS, Revised Supplementary Release on Saturation Properties of
    //     Ordinary Water Substance September 1992,
    //     http://www.iapws.org/relguide/Supp-sat.html, derived from Eq.1
    //

    public static double _dPdT_sat(Dictionary<string, double> cls, double T)
    {
      var Tita = 1 - T / cls.Tc;
      var suma1 = 0;
      var suma2 = 0;
      foreach (var _tup_1 in zip(cls._Pv["ao"], cls._Pv["exp"]))
      {
        var n = _tup_1.Item1;
        var x = _tup_1.Item2;
        suma1 -= n * x * Math.Pow(Tita, x - 1) / cls.Tc;
        suma2 += n * Math.Pow(Tita, x);
      }
      var Pr = (cls.Tc * suma1 / T - cls.Tc / Math.Pow(T, 2) * suma2) * Math.Exp(cls.Tc / T * suma2);
      var dPdT = Pr * cls.Pc;
      return dPdT;
    }
  }

  // Implementation of IAPWS Formulation 1995 for ordinary water substance,
  //   (revised release of 2016), for internal procedures, see MEoS base class
  //
  //   Examples
  //   --------
  //   >>> water=IAPWS95(T=300, rho=996.5560)
  //   >>> water.P, water.cv, water.w, water.s
  //   0.0992418350 4.13018112 1501.51914 0.393062643
  //
  //   >>> water=IAPWS95(T=500, rho=0.435)
  //   >>> water.P, water.cv, water.w, water.s
  //   0.0999679423 1.50817541 548.31425 7.944882714
  //
  //   >>> water=IAPWS95(T=900., P=700)
  //   >>> water.rho, water.cv, water.w, water.s
  //   870.7690 2.66422350 2019.33608 4.17223802
  //
  //   >>> water=IAPWS95(T=300., P=0.1)
  //   >>> water.P, water.rho, water.h, water.s, water.cp, water.w, water.virialB
  //   0.10000 996.56 112.65 0.39306 4.1806 1501.5 -0.066682
  //
  //   >>> water=IAPWS95(T=500., P=0.1)
  //   >>> water.P, water.rho, water.h, water.s, water.cp, water.w, water.virialB
  //   0.10000 0.43514 2928.6 7.9447 1.9813 548.31 -0.0094137
  //
  //   >>> water=IAPWS95(T=450., x=0.5)
  //   >>> water.T, water.P, water.rho, water.h, water.s, water.virialB
  //   450.00 0.93220 9.5723 1761.8 4.3589 -0.013028
  //
  //   >>> water=IAPWS95(P=1.5, rho=1000.)
  //   >>> water.T, water.rho, water.h, water.s, water.cp, water.w, water.virialB
  //   286.44 1000.0 57.253 0.19931 4.1855 1462.1 -0.085566
  //
  //   >>> water=IAPWS95(h=3000, s=8.)
  //   >>> water.T, water.P, water.h, water.s, water.cp, water.w, water.virialB
  //   536.24 0.11970 3000.0 8.0000 1.9984 567.04 -0.0076606
  //
  //   >>> water=IAPWS95(h=150, s=0.4)
  //   >>> water.T, water.P, water.rho, water.h, water.s, water.cp, water.w
  //   301.27 35.50549 1011.48 150.00 0.40000 4.0932 1564.1
  //
  //   >>> water=IAPWS95(T=450., rho=300)
  //   >>> water.T, water.P, water.rho, water.h, water.s, water.x, water.virialB
  //   450.00 0.93220 300.00 770.82 2.1568 0.010693 -0.013028
  //
  //   >>> water=IAPWS95(rho=300., P=0.1)
  //   >>> water.T, water.P, water.rho, water.h, water.s, water.x, water.virialB
  //   372.76 0.10000 300.00 420.56 1.3110 0.0013528 -0.025144
  //
  //   >>> water=IAPWS95(h=1500., P=0.1)
  //   >>> water.T, water.P, water.rho, water.h, water.s, water.x, water.virialB
  //   372.76 0.10000 1.2303 1500.0 4.2068 0.47952 -0.025144
  //
  //   >>> water=IAPWS95(s=5., P=3.5)
  //   >>> water.T, water.P, water.rho, water.h, water.s, water.x, water.virialB
  //   515.71 3.5000 25.912 2222.8 5.0000 0.66921 -0.0085877
  //
  //   >>> water=IAPWS95(T=500., u=900)
  //   >>> water.P, water.rho, water.u, water.h, water.s, water.cp, water.w
  //   108.21 903.62 900.00 1019.8 2.4271 4.1751 1576.0
  //
  //   >>> water=IAPWS95(P=0.3, u=1550.)
  //   >>> water.T, water.P, water.rho, water.u, water.h, water.s, water.x
  //   406.67 0.30000 3.3029 1550.0 1640.8 4.3260 0.49893
  //
  //   >>> water=IAPWS95(rho=300, h=1000.)
  //   >>> water.T, water.P, water.rho, water.u, water.h, water.s, water.x
  //   494.92 2.3991 300.00 992.00 1000.0 2.6315 0.026071
  //
  //   >>> water=IAPWS95(rho=30, s=8.)
  //   >>> water.T, water.P, water.rho, water.u, water.h, water.s, water.cp
  //   1562.42 21.671 30.000 4628.5 5350.9 8.0000 2.7190
  //
  //   >>> water=IAPWS95(rho=30, s=4.)
  //   >>> water.T, water.P, water.rho, water.u, water.h, water.s, water.x
  //   495.00 2.4029 30.000 1597.3 1677.4 4.0000 0.39218
  //
  //   >>> water=IAPWS95(rho=300, u=1000.)
  //   >>> water.T, water.P, water.rho, water.u, water.h, water.s, water.x
  //   496.44 2.4691 300.000 1000.0 1008.2 2.6476 0.02680
  //
  //   >>> water=IAPWS95(s=3., h=1000.)
  //   >>> water.T, water.P, water.rho, water.u, water.h, water.s, water.x
  //   345.73 0.034850 0.73526 952.60 1000.0 3.0000 0.29920
  //
  //   >>> water=IAPWS95(u=995., h=1000.)
  //   >>> water.T, water.P, water.rho, water.u, water.h, water.s, water.x
  //   501.89 2.7329 546.58 995.00 1000.0 2.6298 0.00866
  //
  //   >>> water=IAPWS95(u=1000., s=3.)
  //   >>> water.T, water.P, water.rho, water.u, water.h, water.s, water.x
  //   371.24 0.094712 1.99072 1000.00 1047.6 3.0000 0.28144
  //
  //   References
  //   ----------
  //   IAPWS, Revised Release on the IAPWS Formulation 1995 for the Thermodynamic
  //   Properties of Ordinary Water Substance for General and Scientific Use,
  //   September 2016, http://www.iapws.org/relguide/IAPWS-95.html
  //
  //   IAPWS, Revised Supplementary Release on Saturation Properties of Ordinary
  //   Water Substance September 1992, http://www.iapws.org/relguide/Supp-sat.html
  //
  //   IAPWS, Guideline on a Low-Temperature Extension of the IAPWS-95 Formulation
  //   for Water Vapor, http://www.iapws.org/relguide/LowT.html
  //
  //   IAPWS, Revised Advisory Note No. 3: Thermodynamic Derivatives from IAPWS
  //   Formulations, http://www.iapws.org/relguide/Advise3.pdf
  //

  public class IAPWS95 : MEoS
  {
    public string name = "water";
    public string CASNumber = "7732-18-5";
    public string formula = "H2O";
    public string synonym = "R-718";
    public double Tc = Tc;
    public double rhoc = rhoc;
    public double Pc = Pc;
    public double M = M;
    public double Tt = 273.16;
    public double Tb = 373.1243;
    public double f_acent = 0.3443;
    public double momentoDipolar = 1.855;

    public Dictionary<string, object> Fi0 = new Dictionary<object, object>
    {
      { "ao_log", new List<double> { 1, 3.00632 }},
      { "pow", new List<int> { 0, 1 }},
      { "ao_pow", new List<double> { -8.3204464837497, 6.6832105275932 }},
      { "ao_exp", new List<double> { 0.012436, 0.97315, 1.2795, 0.96956, 0.24873 }},
      { "titao", new List<double> { 1.28728967, 3.53734222, 7.74073708, 9.24437796, 27.5075105 }}
    };

    public Dictionary<string, object> _constants = new Dictionary<object, object>
    {
      { "R", 8.314371357587},
      { "nr1", new List<double> { 0.012533547935523, 7.8957634722828, -8.7803203303561, 0.31802509345418, -0.26145533859358, -0.0078199751687981, 0.0088089493102134 }},
      { "d1", new List<int> { 1, 1, 1, 2, 2, 3, 4 }},
      { "t1", new List<double> { -0.5, 0.875, 1, 0.5, 0.75, 0.375, 1 }},
      { "nr2", new List<double>
        {
          -0.66856572307965,
          0.20433810950965,
          -6.6212605039687E-05,
          -0.19232721156002,
          -0.25709043003438,
          0.16074868486251,
          -0.04009282892587,
          3.9343422603254E-07,
          -7.5941377088144E-06,
          0.00056250979351888,
          -1.5608652257135E-05,
          1.1537996422951E-09,
          3.6582165144204E-07,
          -1.3251180074668E-12,
          -6.2639586912454E-10,
          -0.10793600908932,
          0.017611491008752,
          0.22132295167546,
          -0.40247669763528,
          0.58083399985759,
          0.0049969146990806,
          -0.031358700712549,
          -0.74315929710341,
          0.4780732991548,
          0.020527940895948,
          -0.13636435110343,
          0.014180634400617,
          0.0083326504880713,
          -0.029052336009585,
          0.038615085574206,
          -0.020393486513704,
          -0.0016554050063734,
          0.0019955571979541,
          0.00015870308324157,
          -1.638856834253E-05,
          0.043613615723811,
          0.034994005463765,
          -0.076788197844621,
          0.022446277332006,
          -6.2689710414685E-05,
          -5.5711118565645E-10,
          -0.19905718354408,
          0.31777497330738,
          -0.11841182425981
        }},
      { "c2", new List<int>
        {
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
          1,
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
          2,
          2,
          2,
          2,
          2,
          2,
          2,
          2,
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
          6,
          6,
          6,
          6
        }},
      { "d2", new List<int>
        {
          1,
          1,
          1,
          2,
          2,
          3,
          4,
          4,
          5,
          7,
          9,
          10,
          11,
          13,
          15,
          1,
          2,
          2,
          2,
          3,
          4,
          4,
          4,
          5,
          6,
          6,
          7,
          9,
          9,
          9,
          9,
          9,
          10,
          10,
          12,
          3,
          4,
          4,
          5,
          14,
          3,
          6,
          6,
          6
        }},
      { "t2", new List<int>
        {
          4,
          6,
          12,
          1,
          5,
          4,
          2,
          13,
          9,
          3,
          4,
          11,
          4,
          13,
          1,
          7,
          1,
          9,
          10,
          10,
          3,
          7,
          10,
          10,
          6,
          10,
          10,
          1,
          2,
          3,
          4,
          8,
          6,
          9,
          8,
          16,
          22,
          23,
          23,
          10,
          50,
          44,
          46,
          50
        }},
      { "gamma2", new List<int> { 1 } * 44},
      { "nr3", new List<double> { -31.306260323435, 31.546140237781, -2521.3154341695 }},
      { "d3", new List<int> { 3 } * 3},
      { "t3", new List<int> { 0, 1, 4 }},
      { "alfa3", new List<int> { 20 } * 3},
      { "beta3", new List<int> { 150, 150, 250 }},
      { "gamma3", new List<double> { 1.21, 1.21, 1.25 }},
      { "epsilon3", new List<double> { 1.0 } * 3},
      { "nr4", new List<double> { -0.14874640856724, 0.31806110878444 }},
      { "a4", new List<double> { 3.5, 3.5 }},
      { "b4", new List<double> { 0.85, 0.95 }},
      { "B", new List<double> { 0.2, 0.2 }},
      { "C", new List<int> { 28, 32 }},
      { "D", new List<int> { 700, 800 }},
      { "A", new List<double> { 0.32, 0.32 }},
      { "beta4", new List<double> { 0.3, 0.3
        }}};

    public Dictionary<string, List<double>> _Pv = new Dictionary<object, object>
    {
      { "ao", new List<double>
      { -7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502 }},
      { "exp", new List<double>
      { 1, 1.5, 3, 3.5, 4, 7.5 }}
    };

    public Dictionary<string, object> _rhoL = new Dictionary<object, object>
    {
      { "eq", 2},
      { "ao", new List<double>
      { 1.99274064, 1.09965342, -0.510839303, -1.75493479, -45.5170352, -674694.45 }},
      { "exp", new List<int>
      { 1, 2, 5, 16, 43, 110 }}
    };

    public Dictionary<string, object> _rhoG = new Dictionary<object, object>
    {
      { "eq", 4},
      { "ao", new List<double>
      { -2.0315024, -2.6830294, -5.38626492, -17.2991605, -44.7586581, -63.9201063 }},
      { "exp", new List<double>
      { 1, 2, 4, 9, 18.5, 35.5 }}
    };

    // Low temperature extension of the IAPWS-95
    public virtual object _phi0(object tau, object delta)
    {
      var prop = MEoS._phi0(this, tau, delta);
      var T = this.Tc / tau;
      if(50 <= T < 130)
      {
        var _tup_1 = this._phiex(T);
        var fex = _tup_1.Item1;
        var fext = _tup_1.Item2;
        var fextt = _tup_1.Item3;
        prop["fio"] += fex;
        prop["fiot"] += fext;
        prop["fiott"] += fextt;
      }
      return prop;
    }

    // Low temperature extension
    public virtual object _phiex(object T)
    {
      var tau = this.Tc / T;
      var E = 0.278296458178592;
      var ep = this.Tc / 130;
      var fex = E * (-1 / 2 / tau - 3 / Math.Pow(ep, 2) * (tau + ep) * Math.Log(tau / ep) - 9 / 2 / ep + 9 * tau / 2 / Math.Pow(ep, 2) + Math.Pow(tau, 2) / 2 / Math.Pow(ep, 3));
      var fext = E * (1 / 2 / Math.Pow(tau, 2) - 3 / tau / ep - 3 / Math.Pow(ep, 2) * Math.Log(tau / ep) + 3 / 2 / Math.Pow(ep, 2) + tau / Math.Pow(ep, 3));
      var fextt = E * Math.Pow(-1 / tau + 1 / ep, 3);
      return Tuple.Create(fex, fext, fextt);
    }

    // Auxiliary equation for the alfa coefficient for calculate the
    //     enthalpy along the saturation line
    //
    //     Parameters
    //     ----------
    //     T : float
    //       Temperature, [K]
    //
    //     Returns
    //     -------
    //     alfa : float
    //       alfa coefficient, [kJ/kg]
    //
    //     References
    //     ----------
    //     IAPWS, Revised Supplementary Release on Saturation Properties of
    //     Ordinary Water Substance September 1992,
    //     http://www.iapws.org/relguide/Supp-sat.html, Eq.4
    //

    public static object _alfa_sat(object cls, object T)
    {
      var di = new List<double>
      { -1135.905627715, -5.65134998E-08, 2690.66631, 127.287297, -135.003439, 0.981825814 };
      var expi = new List<double>
      { 0, -19, 1, 4.5, 5, 54.5 };
      var Tita = T / cls.Tc;
      var alfa = 0;
      foreach (var _tup_1 in zip(di, expi))
      {
        var d = _tup_1.Item1;
        var x = _tup_1.Item2;
        alfa += d * Math.Pow(Tita, x);
      }
      return alfa;
    }

    // Auxiliary equation for the phi coefficient for calculate the
    //     entropy along the saturation line
    //
    //     Parameters
    //     ----------
    //     T : float
    //       Temperature, [K]
    //
    //     Returns
    //     -------
    //     phi : float
    //       phi coefficient, [kJ/kgK]
    //
    //     References
    //     ----------
    //     IAPWS, Revised Supplementary Release on Saturation Properties of
    //     Ordinary Water Substance September 1992,
    //     http://www.iapws.org/relguide/Supp-sat.html, Eq.5
    //

    public static object _phi_sat(object cls, object T)
    {
      var di = new List<double>
      { 2319.5246, -5.65134998E-08 * 19 / 20, 2690.66631, 127.287297 * 9 / 7, -135.003439 * 5 / 4, 0.981825814 * 109 / 107 };
      var expi = new List<double>
      { 0, -20, null, 3.5, 4, 53.5 };
      var Tita = T / cls.Tc;
      var suma = 0;
      foreach(var _tup_1 in zip(di, expi))
      {
        var d = _tup_1.Item1;
        var x = _tup_1.Item2;
        if(x == null)
        { suma += d * Math.Log(Tita); }
        else
        { suma += d * Math.Pow(Tita, x); }
      }
      var phi = suma / cls.Tc;
      return phi;
    }

    // Auxiliary equation for the specific enthalpy for saturated liquid
    //
    //     Parameters
    //     ----------
    //     T : float
    //       Temperature, [K]
    //
    //     Returns
    //     -------
    //     h : float
    //       Saturated liquid enthalpy, [kJ/kg]
    //
    //     References
    //     ----------
    //     IAPWS, Revised Supplementary Release on Saturation Properties of
    //     Ordinary Water Substance September 1992,
    //     http://www.iapws.org/relguide/Supp-sat.html, Eq.6
    //

    public static object _Liquid_Enthalpy(object cls, object T)
    {
      var alfa = cls._alfa_sat(T);
      var rho = cls._Liquid_Density(T);
      var dpdT = cls._dPdT_sat(T);
      var h = alfa + T / rho * dpdT * 1000;
      return h;
    }

    // Auxiliary equation for the specific enthalpy for saturated vapor
    //
    //     Parameters
    //     ----------
    //     T : float
    //       Temperature, [K]
    //
    //     Returns
    //     -------
    //     h : float
    //       Saturated vapor enthalpy, [kJ/kg]
    //
    //     References
    //     ----------
    //     IAPWS, Revised Supplementary Release on Saturation Properties of
    //     Ordinary Water Substance September 1992,
    //     http://www.iapws.org/relguide/Supp-sat.html, Eq.7
    //

    public static object _Vapor_Enthalpy(object cls, object T)
    {
      var alfa = cls._alfa_sat(T);
      var rho = cls._Vapor_Density(T);
      var dpdT = cls._dPdT_sat(T);
      var h = alfa + T / rho * dpdT * 1000;
      return h;
    }

    // Auxiliary equation for the specific entropy for saturated liquid
    //
    //     Parameters
    //     ----------
    //     T : float
    //       Temperature, [K]
    //
    //     Returns
    //     -------
    //     s : float
    //       Saturated liquid entropy, [kJ/kgK]
    //
    //     References
    //     ----------
    //     IAPWS, Revised Supplementary Release on Saturation Properties of
    //     Ordinary Water Substance September 1992,
    //     http://www.iapws.org/relguide/Supp-sat.html, Eq.8
    //

    public static object _Liquid_Entropy(object cls, object T)
    {
      var phi = cls._phi_sat(T);
      var rho = cls._Liquid_Density(T);
      var dpdT = cls._dPdT_sat(T);
      var s = phi + dpdT / rho * 1000;
      return s;
    }

    // Auxiliary equation for the specific entropy for saturated vapor
    //
    //     Parameters
    //     ----------
    //     T : float
    //       Temperature, [K]
    //
    //     Returns
    //     -------
    //     s : float
    //       Saturated liquid entropy, [kJ/kgK]
    //
    //     References
    //     ----------
    //     IAPWS, Revised Supplementary Release on Saturation Properties of
    //     Ordinary Water Substance September 1992,
    //     http://www.iapws.org/relguide/Supp-sat.html, Eq.9
    //

    public static object _Vapor_Entropy(object cls, object T)
    {
      var phi = cls._phi_sat(T);
      var rho = cls._Vapor_Density(T);
      var dpdT = cls._dPdT_sat(T);
      var s = phi + dpdT / rho * 1000;
      return s;
    }

    public virtual object _visco(object rho, object T, object fase)
    {
      var @ref = new IAPWS95();
      var st = @ref._Helmholtz(rho, 1.5 * Tc);
      var delta = rho / rhoc;
      var drho = 1000.0 / this.R / 1.5 / Tc / (1 + 2 * delta * st["fird"] + Math.Pow(delta, 2) * st["firdd"]);
      return _iapws._Viscosity(rho, T, fase, drho);
    }

    public virtual object _thermo(object rho, object T, object fase)
    {
      var @ref = new IAPWS95();
      var st = @ref._Helmholtz(rho, 1.5 * Tc);
      var delta = rho / rhoc;
      var drho = 1000.0 / this.R / 1.5 / Tc / (1 + 2 * delta * st["fird"] + Math.Pow(delta, 2) * st["firdd"]);
      return _iapws._ThCond(rho, T, fase, drho);
    }

    public virtual object _surface(object T)
    {
      var s = _iapws._Tension(T);
      return s;
    }
  }

  public class IAPWS95_PT : IAPWS95 { public IAPWS95_PT(object P, object T) : base(T: T, P: P) { } } // Derivated class for direct P and T input
  public class IAPWS95_Ph : IAPWS95 { public IAPWS95_Ph(object P, object h) : base(P: P, h: h) { } } // Derivated class for direct P and h input
  public class IAPWS95_Ps : IAPWS95 { public IAPWS95_Ps(object P, object s) : base(P: P, s: s) { } } // Derivated class for direct P and s input
  public class IAPWS95_Px : IAPWS95 { public IAPWS95_Px(object P, object x) : base(P: P, x: x) { } } // Derivated class for direct P and v input
  public class IAPWS95_Tx : IAPWS95 { public IAPWS95_Tx(object T, object x) : base(T: T, x: x) { } } // Derivated class for direct T and x input

}

