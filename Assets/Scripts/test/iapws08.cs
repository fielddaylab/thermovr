using System.Collections.Generic;
using System.Collections;
using System;

public static class iapws08
{

  static iapws08()
  {
    //IAPWS standard for Seawater IAPWS08 and related functionality. The module
    //include:
    //
    //:class:`SeaWater`: Global module class with all the functionality integrated
    //
    //Other functionality:
    //   * :func:`_Tb`: Boiling temperature of seawater
    //   * :func:`_Tf`: Freezing temperature of seawater
    //   * :func:`_Triple`: Triple point properties of seawater
    //   * :func:`_OsmoticPressure`: Osmotic pressure of seawater
    //   * :func:`_ThCond_SeaWater`: Thermal conductivity of seawater
    //   * :func:`_solNa2SO4`: Solubility of sodium sulfate in aqueous mixtures of
    //   sodium chloride and sulfuric acid
    //   * :func:`_critNaCl`: Critical locus of aqueous solutions of sodium chloride
  }

  public static double Rm = 8.314472;
  public static double Sn = 0.03516504;
  public static double S_ = Sn * 40 / 35;
  public static double Ms = 31.4038218;
  public static int T_ = 40;
  public static int P_ = 100;
  public static double Po = 0.101325;
  public static double To = 273.15;

  //
  //   Class to model seawater with standard IAPWS-08
  //
  //   Parameters
  //   ----------
  //   T : float
  //     Temperature, [K]
  //   P : float
  //     Pressure, [MPa]
  //   S : float
  //     Salinity, [kg/kg]
  //
  //   fast : bool, default False
  //     Use the Supplementary release SR7-09 to speed up the calculation
  //   IF97 : bool, default False
  //     Use the Advisory Note No. 5 with industrial formulation
  //
  //   Returns
  //   -------
  //   rho : float
  //     Density, [kg/m³]
  //   v : float
  //     Specific volume, [m³/kg]
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //   s : float
  //     Specific entropy, [kJ/kg·K]
  //   u : float
  //     Specific internal energy, [kJ/kg]
  //   g : float
  //     Specific Gibbs free energy, [kJ/kg]
  //   a : float
  //     Specific Helmholtz free energy, [kJ/kg]
  //   cp : float
  //     Specific isobaric heat capacity, [kJ/kg·K]
  //   cv : float
  //     Specific isochoric heat capacity, [kJ/kg·K]
  //   gt : float
  //     Derivative Gibbs energy with temperature, [kJ/kg·K]
  //   gp : float
  //     Derivative Gibbs energy with pressure, [m³/kg]
  //   gtt : float
  //     Derivative Gibbs energy with temperature square, [kJ/kg·K²]
  //   gtp : float
  //     Derivative Gibbs energy with pressure and temperature, [m³/kg·K]
  //   gpp : float
  //     Derivative Gibbs energy with temperature square, [m³/kg·MPa]
  //   gs : float
  //     Derivative Gibbs energy with salinity, [kJ/kg]
  //   gsp : float
  //     Derivative Gibbs energy with salinity and pressure, [m³/kg]
  //   alfav : float
  //     Thermal expansion coefficient, [1/K]
  //   betas : float
  //     Isentropic temperature-pressure coefficient, [K/MPa]
  //   xkappa : float
  //     Isothermal compressibility, [1/MPa]
  //   ks : float
  //     Isentropic compressibility, [1/MPa]
  //   w : float
  //     Sound Speed, [m/s]
  //
  //   m : float
  //     Molality of seawater, [mol/kg]
  //   mu : float
  //     Relative chemical potential, [kJ/kg]
  //   muw : float
  //     Chemical potential of H2O, [kJ/kg]
  //   mus : float
  //     Chemical potential of sea salt, [kJ/kg]
  //   osm : float
  //     Osmotic coefficient, [-]
  //   haline : float
  //     Haline contraction coefficient, [kg/kg]
  //
  //   Notes
  //   ------
  //   :class:`Warning` if input isn't in limit:
  //
  //     * 261 ≤ T ≤ 353
  //     * 0 < P ≤ 100
  //     * 0 ≤ S ≤ 0.12
  //
  //   References
  //   ----------
  //   IAPWS, Release on the IAPWS Formulation 2008 for the Thermodynamic
  //   Properties of Seawater, http://www.iapws.org/relguide/Seawater.html
  //
  //   IAPWS, Supplementary Release on a Computationally Efficient Thermodynamic
  //   Formulation for Liquid Water for Oceanographic Use,
  //   http://www.iapws.org/relguide/OceanLiquid.html
  //
  //   IAPWS, Guideline on the Thermal Conductivity of Seawater,
  //   http://www.iapws.org/relguide/Seawater-ThCond.html
  //
  //   IAPWS, Revised Advisory Note No. 3: Thermodynamic Derivatives from IAPWS
  //   Formulations, http://www.iapws.org/relguide/Advise3.pdf
  //
  //   IAPWS,  Advisory Note No. 5: Industrial Calculation of the Thermodynamic
  //   Properties of Seawater, http://www.iapws.org/relguide/Advise5.html
  //
  //   Examples
  //   --------
  //   >>> salt = iapws.SeaWater(T=300, P=1, S=0.04)
  //   >>> salt.rho
  //   1026.7785717245113
  //   >>> salt.gs
  //   88.56221805501536
  //   >>> salt.haline
  //   0.7311487666026304
  //

  public class SeaWater : object
  {
    public int a;
    public object alfav;
    public object betas;
    public object cp;
    public int cv;
    public object h;
    public None haline;
    public object k;
    public int ks;
    public double m;
    public None mu;
    public None mus;
    public None muw;
    public None osm;
    public object P;
    public double rho;
    public object s;
    public object T;
    public int u;
    public object v;
    public double w;
    public object xkappa;

    public Dictionary<string, object> kwargs = new Dictionary<string, double>
    {
      { "T", 0.0},
      { "P", 0.0},
      { "S", null},
      { "fast", false},
      { "IF97", false}};

    public int status = 0;
    public string msg = "Undefined";

    public SeaWater(Hashtable kwargs)
    {
      this.kwargs = SeaWater.kwargs.copy();
      this.@__call__(kwargs);
    }

    // Make instance callable to can add input parameter one to one
    public virtual object @__call__(Hashtable kwargs)
    {
      this.kwargs.update(kwargs);
      if(this.kwargs["T"] && this.kwargs["P"] && this.kwargs["S"] != null)
      {
        this.status = 1;
        this.calculo();
        this.msg = "";
      }
    }

    // Calculate procedure
    public virtual object calculo()
    {
      object kw;
      var T = this.kwargs["T"];
      var P = this.kwargs["P"];
      var S = this.kwargs["S"];
      this.m = S / (1 - S) / Ms;
      if(this.kwargs["fast"] && T <= 313.15)
      {
        var pw = this._waterSupp(T, P);
      }
      else if(this.kwargs["IF97"])
      {
        pw = this._waterIF97(T, P);
      }
      else
      {
        pw = this._water(T, P);
      }
      var ps = this._saline(T, P, S);
      var prop = new Dictionary<object, object>
      {
      };
      foreach (var key in ps)
      {
        prop[key] = pw[key] + ps[key];
        this.@__setattr__(key, prop[key]);
      }
      this.T = T;
      this.P = P;
      this.rho = 1.0 / prop["gp"];
      this.v = prop["gp"];
      this.s = -prop["gt"];
      this.cp = -T * prop["gtt"];
      this.cv = T * (Math.Pow(prop["gtp"], 2) / prop["gpp"] - prop["gtt"]);
      this.h = prop["g"] - T * prop["gt"];
      this.u = prop["g"] - T * prop["gt"] - P * 1000 * prop["gp"];
      this.a = prop["g"] - P * 1000 * prop["gp"];
      this.alfav = prop["gtp"] / prop["gp"];
      this.betas = -prop["gtp"] / prop["gtt"];
      this.xkappa = -prop["gpp"] / prop["gp"];
      this.ks = (Math.Pow(prop["gtp"], 2) - prop["gtt"] * prop["gpp"]) / prop["gp"] / prop["gtt"];
      this.w = prop["gp"] * Math.Pow(prop["gtt"] * 1000 / (Math.Pow(prop["gtp"], 2) - prop["gtt"] * 1000 * prop["gpp"] * 1E-06), 0.5);
      if(pw.Contains("thcond"))
      {
        kw = pw["thcond"];
      }
      else
      {
        kw = _iapws._ThCond(1 / pw["gp"], T);
      }
      try
      {
        this.k = _iapws._ThCond_SeaWater(T, P, S) + kw;
      } catch (NotImplementedError)
      {
        this.k = null;
      }
      if(S)
      {
        this.mu = prop["gs"];
        this.muw = prop["g"] - S * prop["gs"];
        this.mus = prop["g"] + (1 - S) * prop["gs"];
        this.osm = -ps["g"] - S * prop["gs"] / this.m / Rm / T;
        this.haline = -prop["gsp"] / prop["gp"];
      }
      else
      {
        this.mu = null;
        this.muw = null;
        this.mus = null;
        this.osm = null;
        this.haline = null;
      }
    }

    // Wrapper derivative for custom derived properties
    //     where x, y, z can be: P, T, v, u, h, s, g, a
    public virtual object derivative(object z, object x, object y)
    {
      return _utils.deriv_G(this, z, x, y, this);
    }

    // Get properties of pure water, Table4 pag 8
    public static object _water(object cls, object T, object P)
    {
      var water = iapws95.IAPWS95(P: P, T: T);
      var prop = new Dictionary<string, double> { };
      prop["g"] = water.h - T * water.s;
      prop["gt"] = -water.s;
      prop["gp"] = 1.0 / water.rho;
      prop["gtt"] = -water.cp / T;
      prop["gtp"] = water.betas * water.cp / T;
      prop["gpp"] = -1000000.0 / Math.Pow(water.rho * water.w, 2) - Math.Pow(water.betas, 2) * 1000.0 * water.cp / T;
      prop["gs"] = 0;
      prop["gsp"] = 0;
      prop["thcond"] = water.k;
      return prop;
    }

    public static object _waterIF97(object cls, object T, object P)
    {
      var water = iapws97.IAPWS97(P: P, T: T);
      var betas = water.derivative("T", "P", "s", water);
      var prop = new Dictionary<string, double> { };
      prop["g"] = water.h - T * water.s;
      prop["gt"] = -water.s;
      prop["gp"] = 1.0 / water.rho;
      prop["gtt"] = -water.cp / T;
      prop["gtp"] = betas * water.cp / T;
      prop["gpp"] = -1000000.0 / Math.Pow(water.rho * water.w, 2) - Math.Pow(betas, 2) * 1000.0 * water.cp / T;
      prop["gs"] = 0;
      prop["gsp"] = 0;
      return prop;
    }

    // Get properties of pure water using the supplementary release SR7-09,
    //     Table4 pag 6
    public static object _waterSupp(object cls, object T, object P)
    {
      var tau = (T - 273.15) / 40;
      var pi = (P - 0.101325) / 100;
      var J = new List<int>
      {
        0,
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
        3,
        3,
        4,
        4,
        4,
        4,
        4,
        5,
        5,
        5,
        5,
        5,
        6,
        6,
        6,
        6,
        7,
        7
      };
      var K = new List<int>
      {
        0,
        1,
        2,
        3,
        4,
        5,
        6,
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
        4,
        0,
        1,
        2,
        3,
        4,
        0,
        1,
        2,
        3,
        0,
        1
      };
      var G = new List<double>
      {
        101.342743139674,
        100015.695367145,
        -2544.5765420363,
        284.517778446287,
        -33.3146754253611,
        4.20263108803084,
        -0.546428511471039,
        5.90578347909402,
        -270.983805184062,
        776.153611613101,
        -196.51255088122,
        28.9796526294175,
        -2.13290083518327,
        -12357.785933039,
        1455.0364540468,
        -756.558385769359,
        273.479662323528,
        -55.5604063817218,
        4.34420671917197,
        736.741204151612,
        -672.50778314507,
        499.360390819152,
        -239.545330654412,
        48.8012518593872,
        -1.66307106208905,
        -148.185936433658,
        397.968445406972,
        -301.815380621876,
        152.196371733841,
        -26.3748377232802,
        58.0259125842571,
        -194.618310617595,
        120.520654902025,
        -55.2723052340152,
        6.48190668077221,
        -18.9843846514172,
        63.5113936641785,
        -22.2897317140459,
        8.17060541818112,
        3.05081646487967,
        -9.63108119393062
      };
      var g = 0;
      var gt = 0;
      var gp = 0;
      var gtt = 0;
      var gtp = 0;
      var gpp = 0;
      foreach (var _tup_1 in zip(J, K, G))
      {
        var j = _tup_1.Item1;
        var k = _tup_1.Item2;
        var gi = _tup_1.Item3;
        g += gi * Math.Pow(tau, j) * Math.Pow(pi, k);
        if(j >= 1)
        {
          gt += gi * j * Math.Pow(tau, j - 1) * Math.Pow(pi, k);
        }
        if(k >= 1)
        {
          gp += k * gi * Math.Pow(tau, j) * Math.Pow(pi, k - 1);
        }
        if(j >= 2)
        {
          gtt += j * (j - 1) * gi * Math.Pow(tau, j - 2) * Math.Pow(pi, k);
        }
        if(j >= 1 && k >= 1)
        {
          gtp += j * k * gi * Math.Pow(tau, j - 1) * Math.Pow(pi, k - 1);
        }
        if(k >= 2)
        {
          gpp += k * (k - 1) * gi * Math.Pow(tau, j) * Math.Pow(pi, k - 2);
        }
      }
      var prop = new Dictionary<object, object>
      {
      };
      prop["g"] = g * 0.001;
      prop["gt"] = gt / 40 * 0.001;
      prop["gp"] = gp / 100 * 1E-06;
      prop["gtt"] = gtt / Math.Pow(40, 2) * 0.001;
      prop["gtp"] = gtp / 40 / 100 * 1E-06;
      prop["gpp"] = gpp / Math.Pow(100, 2) * 1E-06;
      prop["gs"] = 0;
      prop["gsp"] = 0;
      return prop;
    }

    // Eq 4
    public static object _saline(object cls, object T, object P, object S)
    {
      // Check input in range of validity
      if(T <= 261 || T > 353 || P <= 0 || P > 100 || S < 0 || S > 0.12)
      {
        Debug.Log("Incoming out of bound");
      }
      var S_ = 0.03516504 * 40 / 35;
      var X = Math.Pow(S / S_, 0.5);
      var tau = (T - 273.15) / 40;
      var pi = (P - 0.101325) / 100;
      var I = new List<int>
      {
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        1,
        2,
        3,
        4,
        5,
        6,
        2,
        3,
        4,
        2,
        3,
        4,
        2,
        3,
        4,
        2,
        4,
        2,
        2,
        3,
        4,
        5,
        2,
        3,
        4,
        2,
        3,
        2,
        3,
        2,
        3,
        2,
        3,
        4,
        2,
        3,
        2,
        3,
        2,
        2,
        2,
        3,
        4,
        2,
        3,
        2,
        3,
        2,
        2,
        2,
        3,
        2,
        2,
        2,
        2,
        2,
        2
      };
      var J = new List<int>
      {
        0,
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
        3,
        3,
        3,
        4,
        4,
        4,
        5,
        5,
        6,
        0,
        0,
        0,
        0,
        1,
        1,
        1,
        2,
        2,
        3,
        3,
        4,
        4,
        0,
        0,
        0,
        1,
        1,
        2,
        2,
        3,
        4,
        0,
        0,
        0,
        1,
        1,
        2,
        2,
        3,
        4,
        0,
        0,
        1,
        2,
        3,
        0,
        1,
        2
      };
      var K = new List<int>
      {
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
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
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        4,
        4,
        4,
        4,
        4,
        5,
        5,
        5
      };
      var G = new List<double>
      {
        5812.81456626732,
        1416.27648484197,
        -2432.14662381794,
        2025.80115603697,
        -1091.66841042967,
        374.60123787784,
        -48.5891069025409,
        851.226734946706,
        168.072408311545,
        -493.407510141682,
        543.835333000098,
        -196.028306689776,
        36.7571622995805,
        880.031352997204,
        -43.0664675978042,
        -68.5572509204491,
        -225.267649263401,
        -10.0227370861875,
        49.3667694856254,
        91.4260447751259,
        0.875600661808945,
        -17.1397577419788,
        -21.6603240875311,
        2.49697009569508,
        2.13016970847183,
        -3310.49154044839,
        199.459603073901,
        -54.7919133532887,
        36.0284195611086,
        729.116529735046,
        -175.292041186547,
        -22.6683558512829,
        -860.764303783977,
        383.058066002476,
        694.244814133268,
        -460.319931801257,
        -297.728741987187,
        234.565187611355,
        384.794152978599,
        -52.2940909281335,
        -4.08193978912261,
        -343.956902961561,
        83.1923927801819,
        337.409530269367,
        -54.1917262517112,
        -204.889641964903,
        74.726141138756,
        -96.5324320107458,
        68.0444942726459,
        -30.1755111971161,
        124.687671116248,
        -29.483064349429,
        -178.314556207638,
        25.6398487389914,
        113.561697840594,
        -36.4872919001588,
        15.8408172766824,
        -3.41251932441282,
        -31.656964386073,
        44.2040358308,
        -11.1282734326413,
        -2.62480156590992,
        7.04658803315449,
        -7.92001547211682
      };
      var g = 0;
      var gt = 0;
      var gp = 0;
      var gtt = 0;
      var gtp = 0;
      var gpp = 0;
      var gs = 0;
      var gsp = 0;
      // Calculate only for some salinity
      if(S != 0)
      {
        foreach (var _tup_1 in zip(I, J, K, G))
        {
          var i = _tup_1.Item1;
          var j = _tup_1.Item2;
          var k = _tup_1.Item3;
          var gi = _tup_1.Item4;
          if(i == 1)
          {
            g += gi * Math.Pow(X, 2) * Math.Log(X) * Math.Pow(tau, j) * Math.Pow(pi, k);
            gs += gi * (2 * Math.Log(X) + 1) * Math.Pow(tau, j) * Math.Pow(pi, k);
          }
          else
          {
            g += gi * Math.Pow(X, i) * Math.Pow(tau, j) * Math.Pow(pi, k);
            gs += i * gi * Math.Pow(X, i - 2) * Math.Pow(tau, j) * Math.Pow(pi, k);
          }
          if(j >= 1)
          {
            if(i == 1) { gt += gi * Math.Pow(X, 2) * Math.Log(X) * j * Math.Pow(tau, j - 1) * Math.Pow(pi, k); }
            else { gt += gi * Math.Pow(X, i) * j * Math.Pow(tau, j - 1) * Math.Pow(pi, k); }
          }
          if(k >= 1)
          {
            gp += k * gi * Math.Pow(X, i) * Math.Pow(tau, j) * Math.Pow(pi, k - 1);
            gsp += i * k * gi * Math.Pow(X, i - 2) * Math.Pow(tau, j) * Math.Pow(pi, k - 1);
          }
          if(j >= 2) { gtt += j * (j - 1) * gi * Math.Pow(X, i) * Math.Pow(tau, j - 2) * Math.Pow(pi, k); }
          if(j >= 1 && k >= 1) { gtp += j * k * gi * Math.Pow(X, i) * Math.Pow(tau, j - 1) * Math.Pow(pi, k - 1); }
          if(k >= 2) { gpp += k * (k - 1) * gi * Math.Pow(X, i) * Math.Pow(tau, j) * Math.Pow(pi, k - 2); }
        }
      }

      var prop = new Dictionary<string, double> { };
      prop["g"] = g * 0.001;
      prop["gt"] = gt / 40 * 0.001;
      prop["gp"] = gp / 100 * 1E-06;
      prop["gtt"] = gtt / Math.Pow(40, 2) * 0.001;
      prop["gtp"] = gtp / 40 / 100 * 1E-06;
      prop["gpp"] = gpp / Math.Pow(100, 2) * 1E-06;
      prop["gs"] = gs / S_ / 2 * 0.001;
      prop["gsp"] = gsp / S_ / 2 / 100 * 1E-06;
      return prop;
    }
  }

  // Procedure to calculate the boiling temperature of seawater
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //   S : float
  //     Salinity, [kg/kg]
  //
  //   Returns
  //   -------
  //   Tb : float
  //     Boiling temperature, [K]
  //
  //   References
  //   ----------
  //   IAPWS,  Advisory Note No. 5: Industrial Calculation of the Thermodynamic
  //   Properties of Seawater, http://www.iapws.org/relguide/Advise5.html, Eq 7
  //

  public static object _Tb(object P, object S)
  {
    Func<object, object> f = T =>
    {
      var pw = iapws97._Region1(T, P);
      var gw = pw["h"] - T * pw["s"];
      var pv = iapws97._Region2(T, P);
      var gv = pv["h"] - T * pv["s"];
      var ps = SeaWater._saline(T, P, S);
      return -ps["g"] + S * ps["gs"] - gw + gv;
    };
    var Tb = fsolve(f, 300)[0];
    return Tb;
  }

  // Procedure to calculate the freezing temperature of seawater
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //   S : float
  //     Salinity, [kg/kg]
  //
  //   Returns
  //   -------
  //   Tf : float
  //     Freezing temperature, [K]
  //
  //   References
  //   ----------
  //   IAPWS,  Advisory Note No. 5: Industrial Calculation of the Thermodynamic
  //   Properties of Seawater, http://www.iapws.org/relguide/Advise5.html, Eq 12
  //

  public static object _Tf(object P, object S)
  {
    Func<object, object> f = T =>
    {
      T = (float)T;
      var pw = iapws97._Region1(T, P);
      var gw = pw["h"] - T * pw["s"];
      var gih = _iapws._Ice(T, P)["g"];
      var ps = SeaWater._saline(T, P, S);
      return -ps["g"] + S * ps["gs"] - gw + gih;
    };
    var Tf = fsolve(f, 300)[0];
    return Tf;
  }

  // Procedure to calculate the triple point pressure and temperature for
  //   seawater
  //
  //   Parameters
  //   ----------
  //   S : float
  //     Salinity, [kg/kg]
  //
  //   Returns
  //   -------
  //   prop : dict
  //     Dictionary with the triple point properties:
  //
  //       * Tt: Triple point temperature, [K]
  //       * Pt: Triple point pressure, [MPa]
  //
  //   References
  //   ----------
  //   IAPWS,  Advisory Note No. 5: Industrial Calculation of the Thermodynamic
  //   Properties of Seawater, http://www.iapws.org/relguide/Advise5.html, Eq 7
  //

  public static object _Triple(object S)
  {
    Func<object, object> f = parr =>
    {
      var _tup_1 = parr;
      var T = _tup_1.Item1;
      var P = _tup_1.Item2;
      var pw = iapws97._Region1(T, P);
      var gw = pw["h"] - T * pw["s"];
      var pv = iapws97._Region2(T, P);
      var gv = pv["h"] - T * pv["s"];
      var gih = _iapws._Ice(T, P)["g"];
      var ps = SeaWater._saline(T, P, S);
      return Tuple.Create(-ps["g"] + S * ps["gs"] - gw + gih, -ps["g"] + S * ps["gs"] - gw + gv);
    };
    var _tup_1 = fsolve(f, new List<double>
    {
      273,
      0.0006
    });
    var Tt = _tup_1.Item1;
    var Pt = _tup_1.Item2;

    var prop = new Dictionary<object, object> { };
    prop["Tt"] = Tt;
    prop["Pt"] = Pt;

    return prop;
  }

  // Procedure to calculate the osmotic pressure of seawater
  //
  //   Parameters
  //   ----------
  //   T : float
  //     Tmperature, [K]
  //   P : float
  //     Pressure, [MPa]
  //   S : float
  //     Salinity, [kg/kg]
  //
  //   Returns
  //   -------
  //   Posm : float
  //     Osmotic pressure, [MPa]
  //
  //   References
  //   ----------
  //   IAPWS,  Advisory Note No. 5: Industrial Calculation of the Thermodynamic
  //   Properties of Seawater, http://www.iapws.org/relguide/Advise5.html, Eq 15
  //

  public static object _OsmoticPressure(object T, object P, object S)
  {
    var pw = iapws97._Region1(T, P);
    var gw = pw["h"] - T * pw["s"];
    Func<object, object> f = Posm =>
    {
      var pw2 = iapws97._Region1(T, P + Posm);
      var gw2 = pw2["h"] - T * pw2["s"];
      var ps = SeaWater._saline(T, P + Posm, S);
      return -ps["g"] + S * ps["gs"] - gw + gw2;
    };
    var Posm = fsolve(f, 0)[0];
    return Posm;
  }

  // Equation for the thermal conductivity of seawater
  //
  //   Parameters
  //   ----------
  //   T : float
  //     Temperature, [K]
  //   P : float
  //     Pressure, [MPa]
  //   S : float
  //     Salinity, [kg/kg]
  //
  //   Returns
  //   -------
  //   k : float
  //     Thermal conductivity excess relative to that of the pure water, [W/mK]
  //
  //   Notes
  //   ------
  //   Raise :class:`NotImplementedError` if input isn't in limit:
  //
  //     * 273.15 ≤ T ≤ 523.15
  //     * 0 ≤ P ≤ 140
  //     * 0 ≤ S ≤ 0.17
  //
  //   Examples
  //   --------
  //   >>> _ThCond_Seawater(293.15, 0.1, 0.035)
  //   -0.00418604
  //
  //   References
  //   ----------
  //   IAPWS, Guideline on the Thermal Conductivity of Seawater,
  //   http://www.iapws.org/relguide/Seawater-ThCond.html
  //

  public static object _ThCond_SeaWater(object T, object P, object S)
  {
    // Check input parameters
    if(T < 273.15 || T > 523.15 || P < 0 || P > 140 || S < 0 || S > 0.17)
    {
      throw new NotImplementedException("Incoming out of bound");
    }
    // Eq 4
    var a1 = -7.180891E-05 + 1.831971E-07 * P;
    var a2 = 0.001048077 - 4.494722E-06 * P;
    // Eq 5
    var b1 = 0.1463375 + 0.0009208586 * P;
    var b2 = -0.003086908 + 1.798489E-05 * P;
    var a = a1 * Math.Exp(a2 * (T - 273.15));
    var b = b1 * Math.Exp(b2 * (T - 273.15));
    // Eq 1
    var DL = a * Math.Pow(1000 * S, 1 + b);
    return DL;
  }

  // Equation for the solubility of sodium sulfate in aqueous mixtures of
  //   sodium chloride and sulfuric acid
  //
  //   Parameters
  //   ----------
  //   T : float
  //     Temperature, [K]
  //   mH2SO4 : float
  //     Molality of sufuric acid, [mol/kg(water)]
  //   mNaCl : float
  //     Molality of sodium chloride, [mol/kg(water)]
  //
  //   Returns
  //   -------
  //   S : float
  //     Molal solutility of sodium sulfate, [mol/kg(water)]
  //
  //   Notes
  //   ------
  //   Raise :class:`NotImplementedError` if input isn't in limit:
  //
  //     * 523.15 ≤ T ≤ 623.15
  //     * 0 ≤ mH2SO4 ≤ 0.75
  //     * 0 ≤ mNaCl ≤ 2.25
  //
  //   Examples
  //   --------
  //   >>> _solNa2SO4(523.15, 0.25, 0.75)
  //   2.68
  //
  //   References
  //   ----------
  //   IAPWS, Solubility of Sodium Sulfate in Aqueous Mixtures of Sodium Chloride
  //   and Sulfuric Acid from Water to Concentrated Solutions,
  //   http://www.iapws.org/relguide/na2so4.pdf
  //

  public static object _solNa2SO4(object T, object mH2SO4, object mNaCl)
  {
    // Check input parameters
    if(T < 523.15 || T > 623.15 || mH2SO4 < 0 || mH2SO4 > 0.75 || mNaCl < 0 || mNaCl > 2.25)
    {
      throw new NotImplementedException("Incoming out of bound");
    }
    var A00 = -0.8085987 * T + 81.4613752 + 0.10537803 * T * Math.Log(T);
    var A10 = 3.4636364 * T - 281.63322 - 0.46779874 * T * Math.Log(T);
    var A20 = -6.0029634 * T + 480.60108 + 0.81382854 * T * Math.Log(T);
    var A30 = 4.4540258 * T - 359.36872 - 0.60306734 * T * Math.Log(T);
    var A01 = 0.4909061 * T - 46.556271 - 0.064612393 * T * Math.Log(T);
    var A02 = -0.002781314 * T + 1.722695 + 1.3319698E-06 * T * Math.Log(T);
    var A03 = -0.014074108 * T + 0.99020227 + 0.0019397832 * T * Math.Log(T);
    var A11 = -0.87146573 * T + 71.808756 + 0.11749585 * T * Math.Log(T);
    var S = A00 + A10 * mH2SO4 + A20 * Math.Pow(mH2SO4, 2) + A30 * Math.Pow(mH2SO4, 3) + A01 * mNaCl + A02 * Math.Pow(mNaCl, 2) + A03 * Math.Pow(mNaCl, 3) + A11 * mH2SO4 * mNaCl;
    return S;
  }

  // Equation for the critical locus of aqueous solutions of sodium chloride
  //
  //   Parameters
  //   ----------
  //   x : float
  //     Mole fraction of NaCl, [-]
  //
  //   Returns
  //   -------
  //   prop : dict
  //     A dictionary withe the properties:
  //
  //       * Tc: critical temperature, [K]
  //       * Pc: critical pressure, [MPa]
  //       * rhoc: critical density, [kg/m³]
  //
  //   Notes
  //   ------
  //   Raise :class:`NotImplementedError` if input isn't in limit:
  //
  //     * 0 ≤ x ≤ 0.12
  //
  //   Examples
  //   --------
  //   >>> _critNaCl(0.1)
  //   975.571016
  //
  //   References
  //   ----------
  //   IAPWS, Revised Guideline on the Critical Locus of Aqueous Solutions of
  //   Sodium Chloride, http://www.iapws.org/relguide/critnacl.html
  //

  public static object _critNaCl(object x)
  {
    // Check input parameters
    if(x < 0 || x > 0.12)
    {
      throw new NotImplementedException("Incoming out of bound");
    }
    var T1 = _iapws.Tc * (1 + 23.0 * x - 330.0 * Math.Pow(x, 1.5) - 1800.0 * Math.Pow(x, 2));
    var T2 = _iapws.Tc * (1 + 17.57 * x - 302.6 * Math.Pow(x, 1.5) + 2838.0 * Math.Pow(x, 2) - 13490.0 * Math.Pow(x, 2.5) + 32780.0 * Math.Pow(x, 3) - 36740.0 * Math.Pow(x, 3.5) + 14370.0 * Math.Pow(x, 4));
    var f1 = (abs(10000 * x - 10 - 1) - abs(10000 * x - 10 + 1)) / 4 + 0.5;
    var f2 = (abs(10000 * x - 10 + 1) - abs(10000 * x - 10 - 1)) / 4 + 0.5;
    // Eq 1
    var tc = f1 * T1 + f2 * T2;
    // Eq 7
    var rc = _iapws.rhoc * (1 + 176.07 * x - 2969.3 * Math.Pow(x, 1.5) + 24886.0 * Math.Pow(x, 2) - 113770.0 * Math.Pow(x, 2.5) + 288470.0 * Math.Pow(x, 3) - 381950.0 * Math.Pow(x, 3.5) + 206330.0 * Math.Pow(x, 4));
    // Eq 8
    var DT = tc - _iapws.Tc;
    var pc = _iapws.Pc * (1 + 0.0091443 * DT + 5.1636E-05 * Math.Pow(DT, 2) - 2.536E-07 * Math.Pow(DT, 3) + 3.6494E-10 * Math.Pow(DT, 4));
    var prop = new Dictionary<string, double> { };
    prop["Tc"] = tc;
    prop["rhoc"] = rc;
    prop["Pc"] = pc;
    return prop;
  }

}

