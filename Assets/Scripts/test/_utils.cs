using System;
using System.Collections.Generic;

public static class _utils
{

  static _utils()
  {
    //Miscelaneous internal utilities. This module include:
    //
    //  * :func:`getphase`: Get phase string of state
    //  * :class:`_fase`: Base class to define a phase state
    //  * :func:`deriv_H`: Calculate generic partial derivative with a fundamental
    //  Helmholtz free energy equation of state
    //  * :func:`deriv_G`: Calculate generic partial derivative with a fundamental
    //  Gibbs free energy equation of state
  }

  //!/usr/bin/python
  // -*- coding: utf-8 -*-
  // Return fluid phase string name
  //
  //   Parameters
  //   ----------
  //   Tc : float
  //     Critical temperature, [K]
  //   Pc : float
  //     Critical pressure, [MPa]
  //   T : float
  //     Temperature, [K]
  //   P : float
  //     Pressure, [MPa]
  //   x : float
  //     Quality, [-]
  //   region: int
  //     Region number, used only for IAPWS97 region definition
  //
  //   Returns
  //   -------
  //   phase : str
  //     Phase name
  //

  public static object getphase( object Tc, object Pc, object T, object P, object x, object region)
  {
    // Avoid round problem
    P = round(P, 8);
    T = round(T, 8);
    if(P > Pc && T > Tc) var phase = "Supercritical fluid";
    else if(T > Tc) phase = "Gas";
    else if(P > Pc) phase = "Compressible liquid";
    else if(P == Pc && T == Tc) phase = "Critical point";
    else if(region == 4 && x == 1) phase = "Saturated vapor";
    else if(region == 4 && x == 0) phase = "Saturated liquid";
    else if(region == 4) phase = "Two phases";
    else if(x == 1) phase = "Vapour";
    else if(x == 0) phase = "Liquid";
    return phase;
  }

  // Class to implement a null phase

  public class _fase : object
  {
    public None v = null;
    public None rho = null;
    public None h = null;
    public None s = null;
    public None u = null;
    public None a = null;
    public None g = null;
    public None cp = null;
    public None cv = null;
    public None cp_cv = null;
    public None w = null;
    public None Z = null;
    public None fi = null;
    public None f = null;
    public None mu = null;
    public None k = null;
    public None nu = null;
    public None Prandt = null;
    public None epsilon = null;
    public None alfa = null;
    public None n = null;
    public None alfap = null;
    public None betap = null;
    public None joule = null;
    public None Gruneisen = null;
    public None alfav = null;
    public None kappa = null;
    public None betas = null;
    public None gamma = null;
    public None Kt = null;
    public None kt = null;
    public None Ks = null;
    public None ks = null;
    public None dpdT_rho = null;
    public None dpdrho_T = null;
    public None drhodT_P = null;
    public None drhodP_T = null;
    public None dhdT_rho = null;
    public None dhdT_P = null;
    public None dhdrho_T = null;
    public None dhdrho_P = null;
    public None dhdP_T = null;
    public None dhdP_rho = null;
    public None Z_rho = null;
    public None IntP = null;
    public None hInput = null;
  }

  // Calculate generic partial derivative
  //   :math:`\left.\frac{\partial z}{\partial x}\right|_{y}` from a fundamental
  //   helmholtz free energy equation of state
  //
  //   Parameters
  //   ----------
  //   state : any python object
  //     Only need to define P and T properties, non phase specific properties
  //   z : str
  //     Name of variables in numerator term of derivatives
  //   x : str
  //     Name of variables in denominator term of derivatives
  //   y : str
  //     Name of constant variable in partial derivaritive
  //   fase : any python object
  //     Define phase specific properties (v, cv, alfap, s, betap)
  //
  //   Notes
  //   -----
  //   x, y and z can be the following values:
  //
  //     * P: Pressure
  //     * T: Temperature
  //     * v: Specific volume
  //     * rho: Density
  //     * u: Internal Energy
  //     * h: Enthalpy
  //     * s: Entropy
  //     * g: Gibbs free energy
  //     * a: Helmholtz free energy
  //
  //   Returns
  //   -------
  //   deriv : float
  //     ∂z/∂x|y
  //
  //   References
  //   ----------
  //   IAPWS, Revised Advisory Note No. 3: Thermodynamic Derivatives from IAPWS
  //   Formulations, http://www.iapws.org/relguide/Advise3.pdf
  //

  public static object deriv_H(object state, object z, object x, object y, object fase)
  {
    // We use the relation between rho and v and his partial derivative
    // ∂v/∂b|c = -1/ρ² ∂ρ/∂b|c
    // ∂a/∂v|c = -ρ² ∂a/∂ρ|c
    var mul = 1;
    if(z == "rho")
    {
      mul = -Math.Pow(fase.rho, 2);
      z = "v";
    }
    if(x == "rho")
    {
      mul = -1 / Math.Pow(fase.rho, 2);
      x = "v";
    }
    if(y == "rho")
    {
      y = "v";
    }
    var dT = new Dictionary<object, object>
    {
      { "P", state.P * 1000 * fase.alfap},
      { "T", 1},
      { "v", 0},
      { "u", fase.cv},
      { "h", fase.cv + state.P * 1000 * fase.v * fase.alfap},
      { "s", fase.cv / state.T},
      { "g", state.P * 1000 * fase.v * fase.alfap - fase.s},
      { "a", -fase.s}
    };
    var dv = new Dictionary<object, object>
    {
      { "P", -state.P * 1000 * fase.betap},
      { "T", 0},
      { "v", 1},
      { "u", state.P * 1000 * (state.T * fase.alfap - 1)},
      { "h", state.P * 1000 * (state.T * fase.alfap - fase.v * fase.betap)},
      { "s", state.P * 1000 * fase.alfap},
      { "g", -state.P * 1000 * fase.v * fase.betap},
      { "a", -state.P * 1000}
    };
    var deriv = (dv[z] * dT[y] - dT[z] * dv[y]) / (dv[x] * dT[y] - dT[x] * dv[y]);
    return mul * deriv;
  }

  // Calculate generic partial derivative
  //   :math:`\left.\frac{\partial z}{\partial x}\right|_{y}` from a fundamental
  //   Gibbs free energy equation of state
  //
  //   Parameters
  //   ----------
  //   state : any python object
  //     Only need to define P and T properties, non phase specific properties
  //   z : str
  //     Name of variables in numerator term of derivatives
  //   x : str
  //     Name of variables in denominator term of derivatives
  //   y : str
  //     Name of constant variable in partial derivaritive
  //   fase : any python object
  //     Define phase specific properties (v, cp, alfav, s, xkappa)
  //
  //   Notes
  //   -----
  //   x, y and z can be the following values:
  //
  //     * P: Pressure
  //     * T: Temperature
  //     * v: Specific volume
  //     * rho: Density
  //     * u: Internal Energy
  //     * h: Enthalpy
  //     * s: Entropy
  //     * g: Gibbs free energy
  //     * a: Helmholtz free energy
  //
  //   Returns
  //   -------
  //   deriv : float
  //     ∂z/∂x|y
  //
  //   References
  //   ----------
  //   IAPWS, Revised Advisory Note No. 3: Thermodynamic Derivatives from IAPWS
  //   Formulations, http://www.iapws.org/relguide/Advise3.pdf
  //

  public static object deriv_G( object state, object z, object x, object y, object fase)
  {
    var mul = 1;
    if(z == "rho")
    {
      mul = -Math.Pow(fase.rho, 2);
      z = "v";
    }
    if(x == "rho")
    {
      mul = -1 / Math.Pow(fase.rho, 2);
      x = "v";
    }
    var dT = new Dictionary<object, object>
    {
      { "P", 0},
      { "T", 1},
      { "v", fase.v * fase.alfav},
      { "u", fase.cp - state.P * 1000 * fase.v * fase.alfav},
      { "h", fase.cp},
      { "s", fase.cp / state.T},
      { "g", -fase.s},
      { "a", -state.P * 1000 * fase.v * fase.alfav - fase.s}
    };
    var dP = new Dictionary<object, object>
    {
      { "P", 1},
      { "T", 0},
      { "v", -fase.v * fase.xkappa},
      { "u", fase.v * (state.P * 1000 * fase.xkappa - state.T * fase.alfav)},
      { "h", fase.v * (1 - state.T * fase.alfav)},
      { "s", -fase.v * fase.alfav},
      { "g", fase.v},
      { "a", state.P * 1000 * fase.v * fase.xkappa}
    };
    var deriv = (dP[z] * dT[y] - dT[z] * dP[y]) / (dP[x] * dT[y] - dT[x] * dP[y]);
    return mul * deriv;
  }

}

