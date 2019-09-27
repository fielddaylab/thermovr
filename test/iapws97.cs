/*
using R = _iapws.R;
using Tc = _iapws.Tc;
using Pc = _iapws.Pc;
using rhoc = _iapws.rhoc;
using Tt = _iapws.Tt;
using Pt = _iapws.Pt;
using Tb = _iapws.Tb;
*/
using System;
using System.Collections.Generic;
using System.Linq;
using System.Collections;

public static class iapws97
{

  static iapws97()
  {
    //IAPWS-IF97 standard implementation
    //
    //.. image:: https://raw.githubusercontent.com/jjgomera/iapws/master/images/iapws97.png
    //  :alt: iapws97
    //
    //The module implement the fundamental equation for the five regions (rectangular
    //boxes) and the backward equation (marked in grey).
    //
    //:class:`IAPWS97`: Global module class with all the functionality integrated
    //
    //Fundamental equations:
    //   * :func:`_Region1`
    //   * :func:`_Region2`
    //   * :func:`_Region3`
    //   * :func:`_Region4`
    //   * :func:`_TSat_P`
    //   * :func:`_PSat_T`
    //   * :func:`_Region5`
    //
    //Backward equations:
    //   * :func:`_Backward1_T_Ph`
    //   * :func:`_Backward1_T_Ps`
    //   * :func:`_Backward1_P_hs`
    //   * :func:`_Backward2_T_Ph`
    //   * :func:`_Backward2_T_Ps`
    //   * :func:`_Backward2_P_hs`
    //   * :func:`_Backward3_T_Ph`
    //   * :func:`_Backward3_T_Ps`
    //   * :func:`_Backward3_P_hs`
    //   * :func:`_Backward3_v_Ph`
    //   * :func:`_Backward3_v_Ps`
    //   * :func:`_Backward3_v_PT`
    //   * :func:`_Backward4_T_hs`
    //
    //Boundary equations:
    //   * :func:`_h13_s`
    //   * :func:`_h3a_s`
    //   * :func:`_h1_s`
    //   * :func:`_t_hs`
    //   * :func:`_PSat_h`
    //   * :func:`_h2ab_s`
    //   * :func:`_h_3ab`
    //   * :func:`_h2c3b_s`
    //   * :func:`_hab_s`
    //   * :func:`_hbc_P`
    //
    //
    //References:
    //
    //IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
    //Thermodynamic Properties of Water and Steam August 2007,
    //http://www.iapws.org/relguide/IF97-Rev.html
    //
    //IAPWS, Revised Supplementary Release on Backward Equations for Pressure
    //as a Function of Enthalpy and Entropy p(h,s) for Regions 1 and 2 of the IAPWS
    //Industrial Formulation 1997 for the Thermodynamic Properties of Water and
    //Steam, http://www.iapws.org/relguide/Supp-PHS12-2014.pdf
    //
    //IAPWS, Revised Supplementary Release on Backward Equations for the
    //Functions T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS
    //Industrial Formulation 1997 for the Thermodynamic Properties of Water and
    //Steam, http://www.iapws.org/relguide/Supp-Tv%28ph,ps%293-2014.pdf
    //
    //IAPWS, Revised Supplementary Release on Backward Equations p(h,s) for
    //Region 3, Equations as a Function of h and s for the Region Boundaries, and an
    //Equation Tsat(h,s) for Region 4 of the IAPWS Industrial Formulation 1997 for
    //the Thermodynamic Properties of Water and Steam,
    //http://www.iapws.org/relguide/Supp-phs3-2014.pdf
    //
    //IAPWS, Revised Supplementary Release on Backward Equations for Specific
    //Volume as a Function of Pressure and Temperature v(p,T) for Region 3 of the
    //IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and
    //Steam, http://www.iapws.org/relguide/Supp-VPT3-2016.pdf
    //
    //IAPWS, Revised Advisory Note No. 3: Thermodynamic Derivatives from IAPWS
    //Formulations, http://www.iapws.org/relguide/Advise3.pdf
    //
    //Wagner, W; Kretzschmar, H-J: International Steam Tables: Properties of
    //Water and Steam Based on the Industrial Formulation IAPWS-IF97; Springer, 2008;
    //doi: 10.1007/978-3-540-74234-0
  }

  public static double sc = 4.41202148223476;
  public static double hc = 2087.54684511715;
  public static double Pmin = 0.000611212677444;
  public static double Ps_623 = 16.5291642526;

  // Boundary Region1-Region3
  // Define the boundary between Region 1 and 3, h=f(s)
  //
  //   Parameters
  //   ----------
  //   s : float
  //     Specific entropy, [kJ/kgK]
  //
  //   Returns
  //   -------
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //
  //   Notes
  //   ------
  //   Raise :class:`NotImplementedError` if input isn't in limit:
  //
  //     * s(100MPa,623.15K) ≤ s ≤ s'(623.15K)
  //
  //   References
  //   ----------
  //   IAPWS, Revised Supplementary Release on Backward Equations p(h,s) for
  //   Region 3, Equations as a Function of h and s for the Region Boundaries, and
  //   an Equation Tsat(h,s) for Region 4 of the IAPWS Industrial Formulation 1997
  //   for the Thermodynamic Properties of Water and Steam,
  //   http://www.iapws.org/relguide/Supp-phs3-2014.pdf. Eq 7
  //
  //   Examples
  //   --------
  //   >>> _h13_s(3.7)
  //   1632.525047
  //   >>> _h13_s(3.5)
  //   1566.104611
  //

  public static object _h13_s(object s)
  {
    // Check input parameters
    if(s < 3.397782955 || s > 3.77828134) { throw new NotImplementedException("Incoming out of bound"); }
    var sigma = s / 3.8;
    var I = new List<int> { 0, 1, 1, 3, 5, 6 };
    var J = new List<int> { 0, -2, 2, -12, -4, -3 };
    var n = new List<double> { 0.913965547600543, -4.30944856041991E-05, 60.3235694765419, 1.17518273082168E-18, 0.220000904781292, -69.0815545851641 };
    var suma = 0;
    foreach (var _tup_1 in zip(I, J, n))
    {
      var i = _tup_1.Item1;
      var j = _tup_1.Item2;
      var ni = _tup_1.Item3;
      suma += ni * Math.Pow(sigma - 0.884, i) * Math.Pow(sigma - 0.864, j);
    }
    return 1700 * suma;
  }

  // Boundary Region2-Region3
  // Define the boundary between Region 2 and 3, P=f(T)
  //
  //   Parameters
  //   ----------
  //   T : float
  //     Temperature, [K]
  //
  //   Returns
  //   -------
  //   P : float
  //     Pressure, [MPa]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  //   Thermodynamic Properties of Water and Steam August 2007,
  //   http://www.iapws.org/relguide/IF97-Rev.html, Eq 5
  //
  //   Examples
  //   --------
  //   >>> _P23_T(623.15)
  //   16.52916425
  //

  public static object _P23_T(object T)
  {
    var n = new List<double> { 348.05185628969, -1.1671859879975, 0.0010192970039326 };
    return n[0] + n[1] * T + n[2] * Math.Pow(T, 2);
  }

  // Define the boundary between Region 2 and 3, T=f(P)
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //
  //   Returns
  //   -------
  //   T : float
  //     Temperature, [K]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  //   Thermodynamic Properties of Water and Steam August 2007,
  //   http://www.iapws.org/relguide/IF97-Rev.html, Eq 5
  //
  //   Examples
  //   --------
  //   >>> _t_P(16.52916425)
  //   623.15
  //

  public static object _t_P(object P)
  {
    var n = new List<double> { 0.0010192970039326, 572.54459862746, 13.9188397787 };
    return n[1] + Math.Pow((P - n[2]) / n[0], 0.5);
  }

  // Define the boundary between Region 2 and 3, T=f(h,s)
  //
  //   Parameters
  //   ----------
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //   s : float
  //     Specific entropy, [kJ/kgK]
  //
  //   Returns
  //   -------
  //   T : float
  //     Temperature, [K]
  //
  //   Notes
  //   ------
  //   Raise :class:`NotImplementedError` if input isn't in limit:
  //
  //     * 5.048096828 ≤ s ≤ 5.260578707
  //     * 2.563592004e3 ≤ h ≤ 2.812942061e3
  //
  //   References
  //   ----------
  //   IAPWS, Revised Supplementary Release on Backward Equations p(h,s) for
  //   Region 3, Equations as a Function of h and s for the Region Boundaries, and
  //   an Equation Tsat(h,s) for Region 4 of the IAPWS Industrial Formulation 1997
  //   for the Thermodynamic Properties of Water and Steam,
  //   http://www.iapws.org/relguide/Supp-phs3-2014.pdf. Eq 8
  //
  //   Examples
  //   --------
  //   >>> _t_hs(2600, 5.1)
  //   713.5259364
  //   >>> _t_hs(2800, 5.2)
  //   817.6202120
  //

  public static object _t_hs(object h, object s)
  {
    // Check input parameters
    if(s < 5.048096828 || s > 5.260578707 || h < 2563.592004 || h > 2812.942061)
    {
      throw new NotImplementedException("Incoming out of bound");
    }
    var nu = h / 3000;
    var sigma = s / 5.3;
    var I = new List<int>
    {
      -12,
      -10,
      -8,
      -4,
      -3,
      -2,
      -2,
      -2,
      -2,
      0,
      1,
      1,
      1,
      3,
      3,
      5,
      6,
      6,
      8,
      8,
      8,
      12,
      12,
      14,
      14
    };
    var J = new List<int>
    {
      10,
      8,
      3,
      4,
      3,
      -6,
      2,
      3,
      4,
      0,
      -3,
      -2,
      10,
      -2,
      -1,
      -5,
      -6,
      -3,
      -8,
      -2,
      -1,
      -12,
      -1,
      -12,
      1
    };
    var n = new List<double>
    {
      0.00062909626082981,
      -0.000823453502583165,
      5.15446951519474E-08,
      -1.17565945784945,
      3.48519684726192,
      -5.07837382408313E-12,
      -2.84637670005479,
      -2.36092263939673,
      6.01492324973779,
      1.48039650824546,
      0.000360075182221907,
      -0.0126700045009952,
      -1221843.32521413,
      0.149276502463272,
      0.698733471798484,
      -0.0252207040114321,
      0.0147151930985213,
      -1.08618917681849,
      -0.000936875039816322,
      81.9877897570217,
      -182.041861521835,
      2.61907376402688E-06,
      -29162.6417025961,
      1.40660774926165E-05,
      7832370.62349385
    };
    var suma = 0;
    foreach (var _tup_1 in zip(I, J, n))
    {
      var i = _tup_1.Item1;
      var j = _tup_1.Item2;
      var ni = _tup_1.Item3;
      suma += ni * Math.Pow(nu - 0.727, i) * Math.Pow(sigma - 0.864, j);
    }
    return 900 * suma;
  }

  // Saturated line
  // Define the saturated line, P=f(T)
  //
  //   Parameters
  //   ----------
  //   T : float
  //     Temperature, [K]
  //
  //   Returns
  //   -------
  //   P : float
  //     Pressure, [MPa]
  //
  //   Notes
  //   ------
  //   Raise :class:`NotImplementedError` if input isn't in limit:
  //
  //     * 273.15 ≤ T ≤ 647.096
  //
  //   References
  //   ----------
  //   IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  //   Thermodynamic Properties of Water and Steam August 2007,
  //   http://www.iapws.org/relguide/IF97-Rev.html, Eq 30
  //
  //   Examples
  //   --------
  //   >>> _PSat_T(500)
  //   2.63889776
  //

  public static object _PSat_T(object T)
  {
    // Check input parameters
    if(T < 273.15 || T > Tc) { throw new NotImplementedException("Incoming out of bound"); }
    var n = new List<double>
    {
      0,
      1167.0521452767,
      -724213.16703206,
      -17.073846940092,
      12020.82470247,
      -3232555.0322333,
      14.91510861353,
      -4823.2657361591,
      405113.40542057,
      -0.23855557567849,
      650.17534844798
    };
    var tita = T + n[9] / (T - n[10]);
    var A = Math.Pow(tita, 2) + n[1] * tita + n[2];
    var B = n[3] * Math.Pow(tita, 2) + n[4] * tita + n[5];
    var C = n[6] * Math.Pow(tita, 2) + n[7] * tita + n[8];
    return Math.Pow(2 * C / (-B + Math.Pow(Math.Pow(B, 2) - 4 * A * C, 0.5)), 4);
  }

  // Define the saturated line, T=f(P)
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //
  //   Returns
  //   -------
  //   T : float
  //     Temperature, [K]
  //
  //   Notes
  //   ------
  //   Raise :class:`NotImplementedError` if input isn't in limit:
  //
  //     * 0.00061121 ≤ P ≤ 22.064
  //
  //   References
  //   ----------
  //   IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  //   Thermodynamic Properties of Water and Steam August 2007,
  //   http://www.iapws.org/relguide/IF97-Rev.html, Eq 31
  //
  //   Examples
  //   --------
  //   >>> _TSat_P(10)
  //   584.149488
  //

  public static object _TSat_P(object P)
  {
    // Check input parameters
    if(P < 611.212677 / 1000000.0 || P > 22.064)
    {
      throw new NotImplementedException("Incoming out of bound");
    }
    var n = new List<double>
    {
      0,
      1167.0521452767,
      -724213.16703206,
      -17.073846940092,
      12020.82470247,
      -3232555.0322333,
      14.91510861353,
      -4823.2657361591,
      405113.40542057,
      -0.23855557567849,
      650.17534844798
    };
    var beta = Math.Pow(P, 0.25);
    var E = Math.Pow(beta, 2) + n[3] * beta + n[6];
    var F = n[1] * Math.Pow(beta, 2) + n[4] * beta + n[7];
    var G = n[2] * Math.Pow(beta, 2) + n[5] * beta + n[8];
    var D = 2 * G / (-F - Math.Pow(Math.Pow(F, 2) - 4 * E * G, 0.5));
    return (n[10] + D - Math.Pow(Math.Pow(n[10] + D, 2) - 4 * (n[9] + n[10] * D), 0.5)) / 2;
  }

  // Define the saturated line, P=f(h) for region 3
  //
  //   Parameters
  //   ----------
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //
  //   Returns
  //   -------
  //   P : float
  //     Pressure, [MPa]
  //
  //   Notes
  //   ------
  //   Raise :class:`NotImplementedError` if input isn't in limit:
  //
  //     * h'(623.15K) ≤ h ≤ h''(623.15K)
  //
  //   References
  //   ----------
  //   IAPWS, Revised Supplementary Release on Backward Equations for the
  //   Functions T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS
  //   Industrial Formulation 1997 for the Thermodynamic Properties of Water and
  //   Steam, http://www.iapws.org/relguide/Supp-Tv%28ph,ps%293-2014.pdf, Eq 10
  //
  //   Examples
  //   --------
  //   >>> _PSat_h(1700)
  //   17.24175718
  //   >>> _PSat_h(2400)
  //   20.18090839
  //

  public static object _PSat_h(object h)
  {
    // Check input parameters
    var hmin_Ps3 = _Region1(623.15, Ps_623)["h"];
    var hmax_Ps3 = _Region2(623.15, Ps_623)["h"];
    if(h < hmin_Ps3 || h > hmax_Ps3)
    {
      throw new NotImplementedException("Incoming out of bound");
    }
    var nu = h / 2600;
    var I = new List<int>
    {
      0,
      1,
      1,
      1,
      1,
      5,
      7,
      8,
      14,
      20,
      22,
      24,
      28,
      36
    };
    var J = new List<int>
    {
      0,
      1,
      3,
      4,
      36,
      3,
      0,
      24,
      16,
      16,
      3,
      18,
      8,
      24
    };
    var n = new List<double>
    {
      0.600073641753024,
      -9.36203654849857,
      24.6590798594147,
      -107.014222858224,
      -91582131580576.8,
      -8623.32011700662,
      -23.5837344740032,
      2.52304969384128E+17,
      -3.89718771997719E+18,
      -3.33775713645296E+22,
      35649946963.6328,
      -1.48547544720641E+26,
      3.30611514838798E+18,
      8.13641294467829E+37
    };
    var suma = 0;
    foreach (var _tup_1 in zip(I, J, n))
    {
      var i = _tup_1.Item1;
      var j = _tup_1.Item2;
      var ni = _tup_1.Item3;
      suma += ni * Math.Pow(nu - 1.02, i) * Math.Pow(nu - 0.608, j);
    }
    return 22 * suma;
  }

  // Define the saturated line, P=f(s) for region 3
  //
  //   Parameters
  //   ----------
  //   s : float
  //     Specific entropy, [kJ/kgK]
  //
  //   Returns
  //   -------
  //   P : float
  //     Pressure, [MPa]
  //
  //   Notes
  //   ------
  //   Raise :class:`NotImplementedError` if input isn't in limit:
  //
  //     * s'(623.15K) ≤ s ≤ s''(623.15K)
  //
  //   References
  //   ----------
  //   IAPWS, Revised Supplementary Release on Backward Equations for the
  //   Functions T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS
  //   Industrial Formulation 1997 for the Thermodynamic Properties of Water and
  //   Steam, http://www.iapws.org/relguide/Supp-Tv%28ph,ps%293-2014.pdf, Eq 11
  //
  //   Examples
  //   --------
  //   >>> _PSat_s(3.8)
  //   16.87755057
  //   >>> _PSat_s(5.2)
  //   16.68968482
  //

  public static object _PSat_s(object s)
  {
    // Check input parameters
    var smin_Ps3 = _Region1(623.15, Ps_623)["s"];
    var smax_Ps3 = _Region2(623.15, Ps_623)["s"];
    if(s < smin_Ps3 || s > smax_Ps3)
    {
      throw new NotImplementedException("Incoming out of bound");
    }
    var sigma = s / 5.2;
    var I = new List<int>
    {
      0,
      1,
      1,
      4,
      12,
      12,
      16,
      24,
      28,
      32
    };
    var J = new List<int>
    {
      0,
      1,
      32,
      7,
      4,
      14,
      36,
      10,
      0,
      18
    };
    var n = new List<double>
    {
      0.639767553612785,
      -12.9727445396014,
      -2.24595125848403E+15,
      1774667.41801846,
      7170793495.71538,
      -3.78829107169011E+17,
      -9.55586736431328E+34,
      1.87269814676188E+23,
      119254746466.473,
      1.10649277244882E+36
    };
    var suma = 0;
    foreach (var _tup_1 in zip(I, J, n))
    {
      var i = _tup_1.Item1;
      var j = _tup_1.Item2;
      var ni = _tup_1.Item3;
      suma += ni * Math.Pow(sigma - 1.03, i) * Math.Pow(sigma - 0.699, j);
    }
    return 22 * suma;
  }

  // Define the saturated line boundary between Region 1 and 4, h=f(s)
  //
  //   Parameters
  //   ----------
  //   s : float
  //     Specific entropy, [kJ/kgK]
  //
  //   Returns
  //   -------
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //
  //   Notes
  //   ------
  //   Raise :class:`NotImplementedError` if input isn't in limit:
  //
  //     * s'(273.15K) ≤ s ≤ s'(623.15K)
  //
  //   References
  //   ----------
  //   IAPWS, Revised Supplementary Release on Backward Equations p(h,s) for
  //   Region 3, Equations as a Function of h and s for the Region Boundaries, and
  //   an Equation Tsat(h,s) for Region 4 of the IAPWS Industrial Formulation 1997
  //   for the Thermodynamic Properties of Water and Steam,
  //   http://www.iapws.org/relguide/Supp-phs3-2014.pdf. Eq 3
  //
  //   Examples
  //   --------
  //   >>> _h1_s(1)
  //   308.5509647
  //   >>> _h1_s(3)
  //   1198.359754
  //

  public static object _h1_s(object s)
  {
    // Check input parameters
    if(s < -0.0001545495919 || s > 3.77828134) { throw new NotImplementedException("Incoming out of bound"); }

    var sigma = s / 3.8;
    var I = new List<int>
    {
      0,
      0,
      1,
      1,
      2,
      2,
      3,
      3,
      4,
      4,
      4,
      5,
      5,
      7,
      8,
      12,
      12,
      14,
      14,
      16,
      20,
      20,
      22,
      24,
      28,
      32,
      32
    };
    var J = new List<int>
    {
      14,
      36,
      3,
      16,
      0,
      5,
      4,
      36,
      4,
      16,
      24,
      18,
      24,
      1,
      4,
      2,
      4,
      1,
      22,
      10,
      12,
      28,
      8,
      3,
      0,
      6,
      8
    };
    var n = new List<double>
    {
      0.332171191705237,
      0.000611217706323496,
      -8.82092478906822,
      -0.45562819254325,
      -2.63483840850452E-05,
      -22.3949661148062,
      -4.28398660164013,
      -0.616679338856916,
      -14.682303110404,
      284.523138727299,
      -113.398503195444,
      1156.71380760859,
      395.551267359325,
      -1.54891257229285,
      19.4486637751291,
      -3.57915139457043,
      -3.35369414148819,
      -0.66442679633246,
      32332.1885383934,
      3317.66744667084,
      -22350.1257931087,
      5739538.75852936,
      173.226193407919,
      -0.0363968822121321,
      8.34596332878346E-07,
      5.03611916682674,
      65.5444787064505
    };
    var suma = 0;
    foreach (var _tup_1 in zip(I, J, n))
    {
      var i = _tup_1.Item1;
      var j = _tup_1.Item2;
      var ni = _tup_1.Item3;
      suma += ni * Math.Pow(sigma - 1.09, i) * Math.Pow(sigma + 3.66E-05, j);
    }
    return 1700 * suma;
  }

  // Define the saturated line boundary between Region 4 and 3a, h=f(s)
  //
  //   Parameters
  //   ----------
  //   s : float
  //     Specific entropy, [kJ/kgK]
  //
  //   Returns
  //   -------
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //
  //   Notes
  //   ------
  //   Raise :class:`NotImplementedError` if input isn't in limit:
  //
  //     * s'(623.15K) ≤ s ≤ sc
  //
  //   References
  //   ----------
  //   IAPWS, Revised Supplementary Release on Backward Equations p(h,s) for
  //   Region 3, Equations as a Function of h and s for the Region Boundaries, and
  //   an Equation Tsat(h,s) for Region 4 of the IAPWS Industrial Formulation 1997
  //   for the Thermodynamic Properties of Water and Steam,
  //   http://www.iapws.org/relguide/Supp-phs3-2014.pdf. Eq 4
  //
  //   Examples
  //   --------
  //   >>> _h3a_s(3.8)
  //   1685.025565
  //   >>> _h3a_s(4.2)
  //   1949.352563
  //

  public static object _h3a_s(object s)
  {
    // Check input parameters
    if(s < 3.77828134 || s > 4.41202148223476) { throw new NotImplementedException("Incoming out of bound"); }

    var sigma = s / 3.8;
    var I = new List<int>
    {
      0,
      0,
      0,
      0,
      2,
      3,
      4,
      4,
      5,
      5,
      6,
      7,
      7,
      7,
      10,
      10,
      10,
      32,
      32
    };
    var J = new List<int>
    {
      1,
      4,
      10,
      16,
      1,
      36,
      3,
      16,
      20,
      36,
      4,
      2,
      28,
      32,
      14,
      32,
      36,
      0,
      6
    };
    var n = new List<double>
    {
      0.822673364673336,
      0.181977213534479,
      -0.0112000260313624,
      -0.000746778287048033,
      -0.179046263257381,
      0.0424220110836657,
      -0.341355823438768,
      -2.09881740853565,
      -8.22477343323596,
      -4.99684082076008,
      0.191413958471069,
      0.0581062241093136,
      -1655.05498701029,
      1588.70443421201,
      -85.0623535172818,
      -31771.4386511207,
      -94589.0406632871,
      -1.3927384708869E-06,
      0.63105253224098
    };
    var suma = 0;
    foreach (var _tup_1 in zip(I, J, n))
    {
      var i = _tup_1.Item1;
      var j = _tup_1.Item2;
      var ni = _tup_1.Item3;
      suma += ni * Math.Pow(sigma - 1.09, i) * Math.Pow(sigma + 3.66E-05, j);
    }
    return 1700 * suma;
  }

  // Define the saturated line boundary between Region 4 and 2a-2b, h=f(s)
  //
  //   Parameters
  //   ----------
  //   s : float
  //     Specific entropy, [kJ/kgK]
  //
  //   Returns
  //   -------
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //
  //   Notes
  //   ------
  //   Raise :class:`NotImplementedError` if input isn't in limit:
  //
  //     * 5.85 ≤ s ≤ s"(273.15K)
  //
  //   References
  //   ----------
  //   IAPWS, Revised Supplementary Release on Backward Equations p(h,s) for
  //   Region 3, Equations as a Function of h and s for the Region Boundaries, and
  //   an Equation Tsat(h,s) for Region 4 of the IAPWS Industrial Formulation 1997
  //   for the Thermodynamic Properties of Water and Steam,
  //   http://www.iapws.org/relguide/Supp-phs3-2014.pdf. Eq 5
  //
  //   Examples
  //   --------
  //   >>> _h2ab_s(7)
  //   2723.729985
  //   >>> _h2ab_s(9)
  //   2511.861477
  //

  public static object _h2ab_s(object s)
  {
    // Check input parameters
    if(s < 5.85 || s > 9.155759395) { throw new NotImplementedException("Incoming out of bound"); }

    var sigma1 = s / 5.21;
    var sigma2 = s / 9.2;
    var I = new List<int>
    {
      1,
      1,
      2,
      2,
      4,
      4,
      7,
      8,
      8,
      10,
      12,
      12,
      18,
      20,
      24,
      28,
      28,
      28,
      28,
      28,
      32,
      32,
      32,
      32,
      32,
      36,
      36,
      36,
      36,
      36
    };
    var J = new List<int>
    {
      8,
      24,
      4,
      32,
      1,
      2,
      7,
      5,
      12,
      1,
      0,
      7,
      10,
      12,
      32,
      8,
      12,
      20,
      22,
      24,
      2,
      7,
      12,
      14,
      24,
      10,
      12,
      20,
      22,
      28
    };
    var n = new List<double>
    {
      -524.581170928788,
      -9269472.18142218,
      -237.385107491666,
      21077015581.2776,
      -23.9494562010986,
      221.802480294197,
      -5104725.33393438,
      1249813.96109147,
      2000084369.96201,
      -815.158509791035,
      -157.612685637523,
      -11420042233.2791,
      6.62364680776872E+15,
      -2.27622818296144E+18,
      -1.71048081348406E+31,
      6.60788766938091E+15,
      1.66320055886021E+22,
      -2.18003784381501E+29,
      -7.87276140295618E+29,
      1.51062329700346E+31,
      7957321.70300541,
      1.31957647355347E+15,
      -3.2509706829914E+23,
      -4.18600611419248E+25,
      2.97478906557467E+34,
      -9.53588761745473E+19,
      1.66957699620939E+24,
      -1.75407764869978E+32,
      3.47581490626396E+34,
      -7.10971318427851E+38
    };
    var suma = 0;
    foreach (var _tup_1 in zip(I, J, n))
    {
      var i = _tup_1.Item1;
      var j = _tup_1.Item2;
      var ni = _tup_1.Item3;
      suma += ni * Math.Pow(1 / sigma1 - 0.513, i) * Math.Pow(sigma2 - 0.524, j);
    }
    return 2800 * Math.Exp(suma);
  }

  // Define the saturated line boundary between Region 4 and 2c-3b, h=f(s)
  //
  //   Parameters
  //   ----------
  //   s : float
  //     Specific entropy, [kJ/kgK]
  //
  //   Returns
  //   -------
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //
  //   Notes
  //   ------
  //   Raise :class:`NotImplementedError` if input isn't in limit:
  //
  //     * sc ≤ s ≤ 5.85
  //
  //   References
  //   ----------
  //   IAPWS, Revised Supplementary Release on Backward Equations p(h,s) for
  //   Region 3, Equations as a Function of h and s for the Region Boundaries, and
  //   an Equation Tsat(h,s) for Region 4 of the IAPWS Industrial Formulation 1997
  //   for the Thermodynamic Properties of Water and Steam,
  //   http://www.iapws.org/relguide/Supp-phs3-2014.pdf. Eq 6
  //
  //   Examples
  //   --------
  //   >>> _h2c3b_s(5.5)
  //   2687.693850
  //   >>> _h2c3b_s(4.5)
  //   2144.360448
  //

  public static object _h2c3b_s(object s)
  {
    // Check input parameters
    if(s < 4.41202148223476 || s > 5.85) { throw new NotImplementedException("Incoming out of bound"); }

    var sigma = s / 5.9;
    var I = new List<int>
    {
      0,
      0,
      0,
      1,
      1,
      5,
      6,
      7,
      8,
      8,
      12,
      16,
      22,
      22,
      24,
      36
    };
    var J = new List<int>
    {
      0,
      3,
      4,
      0,
      12,
      36,
      12,
      16,
      2,
      20,
      32,
      36,
      2,
      32,
      7,
      20
    };
    var n = new List<double>
    {
      1.04351280732769,
      -2.27807912708513,
      1.80535256723202,
      0.420440834792042,
      -105721.24483466,
      4.36911607493884E+24,
      -328032702839.753,
      -6.7868676080427E+15,
      7439.57464645363,
      -3.56896445355761E+19,
      1.67590585186801E+31,
      -3.55028625419105E+37,
      396611982166.538,
      -4.14716268484468E+40,
      3.59080103867382E+18,
      -1.16994334851995E+40
    };
    var suma = 0;
    foreach (var _tup_1 in zip(I, J, n))
    {
      var i = _tup_1.Item1;
      var j = _tup_1.Item2;
      var ni = _tup_1.Item3;
      suma += ni * Math.Pow(sigma - 1.02, i) * Math.Pow(sigma - 0.726, j);
    }
    return 2800 * Math.Pow(suma, 4);
  }

  // Region 1
  // Basic equation for region 1
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
  //     Dict with calculated properties. The available properties are:
  //
  //       * v: Specific volume, [m³/kg]
  //       * h: Specific enthalpy, [kJ/kg]
  //       * s: Specific entropy, [kJ/kgK]
  //       * cp: Specific isobaric heat capacity, [kJ/kgK]
  //       * cv: Specific isocoric heat capacity, [kJ/kgK]
  //       * w: Speed of sound, [m/s]
  //       * alfav: Cubic expansion coefficient, [1/K]
  //       * kt: Isothermal compressibility, [1/MPa]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  //   Thermodynamic Properties of Water and Steam August 2007,
  //   http://www.iapws.org/relguide/IF97-Rev.html, Eq 7
  //
  //   Examples
  //   --------
  //   >>> _Region1(300,3)["v"]
  //   0.00100215168
  //   >>> _Region1(300,3)["h"]
  //   115.331273
  //   >>> _Region1(300,3)["h"]-3000*_Region1(300,3)["v"]
  //   112.324818
  //   >>> _Region1(300,80)["s"]
  //   0.368563852
  //   >>> _Region1(300,80)["cp"]
  //   4.01008987
  //   >>> _Region1(300,80)["cv"]
  //   3.91736606
  //   >>> _Region1(500,3)["w"]
  //   1240.71337
  //   >>> _Region1(500,3)["alfav"]
  //   0.00164118128
  //   >>> _Region1(500,3)["kt"]
  //   0.00112892188
  //

  public static Dictionary<string, double> _Region1(double T, double P)
  {
    if(P < 0) { P = Pmin; }

    var I = new List<int>
    {
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
      2,
      2,
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
      8,
      8,
      21,
      23,
      29,
      30,
      31,
      32
    };
    var J = new List<int>
    {
      -2,
      -1,
      0,
      1,
      2,
      3,
      4,
      5,
      -9,
      -7,
      -1,
      0,
      1,
      3,
      -3,
      0,
      1,
      3,
      17,
      -4,
      0,
      6,
      -5,
      -2,
      10,
      -8,
      -11,
      -6,
      -29,
      -31,
      -38,
      -39,
      -40,
      -41
    };
    var n = new List<double>
    {
      0.14632971213167,
      -0.84548187169114,
      -3.756360367204,
      3.3855169168385,
      -0.95791963387872,
      0.15772038513228,
      -0.016616417199501,
      0.00081214629983568,
      0.00028319080123804,
      -0.00060706301565874,
      -0.018990068218419,
      -0.032529748770505,
      -0.021841717175414,
      -5.283835796993E-05,
      -0.00047184321073267,
      -0.00030001780793026,
      4.7661393906987E-05,
      -4.4141845330846E-06,
      -7.2694996297594E-16,
      -3.1679644845054E-05,
      -2.8270797985312E-06,
      -8.5205128120103E-10,
      -2.2425281908E-06,
      -6.5171222895601E-07,
      -1.4341729937924E-13,
      -4.0516996860117E-07,
      -1.2734301741641E-09,
      -1.7424871230634E-10,
      -6.8762131295531E-19,
      1.4478307828521E-20,
      2.6335781662795E-23,
      -1.1947622640071E-23,
      1.8228094581404E-24,
      -9.3537087292458E-26
    };
    var Tr = 1386 / T;
    var Pr = P / 16.53;
    var g = 0;
    foreach (var _tup_1 in zip(I, J, n))
    {
      var i = _tup_1.Item1;
      var j = _tup_1.Item2;
      var ni = _tup_1.Item3;
      g += ni * Math.Pow(7.1 - Pr, i) * Math.Pow(Tr - 1.222, j);
      gp -= ni * i * Math.Pow(7.1 - Pr, i - 1) * Math.Pow(Tr - 1.222, j);
      gpp += ni * i * (i - 1) * Math.Pow(7.1 - Pr, i - 2) * Math.Pow(Tr - 1.222, j);
      gt += ni * j * Math.Pow(7.1 - Pr, i) * Math.Pow(Tr - 1.222, j - 1);
      gtt += ni * j * (j - 1) * Math.Pow(7.1 - Pr, i) * Math.Pow(Tr - 1.222, j - 2);
      gpt -= ni * i * j * Math.Pow(7.1 - Pr, i - 1) * Math.Pow(Tr - 1.222, j - 1);
    }

    var propiedades = new Dictionary<string, double> { };
    propiedades["T"] = T;
    propiedades["P"] = P;
    propiedades["v"] = Pr * gp * R * T / P / 1000;
    propiedades["h"] = Tr * gt * R * T;
    propiedades["s"] = R * (Tr * gt - g);
    propiedades["cp"] = -R * Math.Pow(Tr, 2) * gtt;
    propiedades["cv"] = R * (-Math.Pow(Tr, 2) * gtt + Math.Pow(gp - Tr * gpt, 2) / gpp);
    propiedades["w"] = Math.Sqrt(R * T * 1000 * Math.Pow(gp, 2) / (Math.Pow(gp - Tr * gpt, 2) / (Math.Pow(Tr, 2) * gtt) - gpp));
    propiedades["alfav"] = (1 - Tr * gpt / gp) / T;
    propiedades["kt"] = -Pr * gpp / gp / P;
    propiedades["region"] = 1;
    propiedades["x"] = 0;
    return propiedades;
  }

  //
  //   Backward equation for region 1, T=f(P,h)
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //
  //   Returns
  //   -------
  //   T : float
  //     Temperature, [K]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  //   Thermodynamic Properties of Water and Steam August 2007,
  //   http://www.iapws.org/relguide/IF97-Rev.html, Eq 11
  //
  //   Examples
  //   --------
  //   >>> _Backward1_T_Ph(3,500)
  //   391.798509
  //   >>> _Backward1_T_Ph(80,1500)
  //   611.041229
  //

  public static object _Backward1_T_Ph(object P, object h)
  {
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
      1,
      2,
      2,
      3,
      3,
      4,
      5,
      6
    };
    var J = new List<int>
    {
      0,
      1,
      2,
      6,
      22,
      32,
      0,
      1,
      2,
      3,
      4,
      10,
      32,
      10,
      32,
      10,
      32,
      32,
      32,
      32
    };
    var n = new List<double>
    {
      -238.72489924521,
      404.21188637945,
      113.49746881718,
      -5.8457616048039,
      -0.0001528548241314,
      -1.0866707695377E-06,
      -13.391744872602,
      43.211039183559,
      -54.010067170506,
      30.535892203916,
      -6.5964749423638,
      0.0093965400878363,
      1.157364750534E-07,
      -2.5858641282073E-05,
      -4.0644363084799E-09,
      6.6456186191635E-08,
      8.0670734103027E-11,
      -9.3477771213947E-13,
      5.8265442020601E-15,
      -1.5020185953503E-17
    };
    var Pr = P / 1;
    var nu = h / 2500;
    var T = 0;
    foreach (var _tup_1 in zip(I, J, n))
    {
      var i = _tup_1.Item1;
      var j = _tup_1.Item2;
      var ni = _tup_1.Item3;
      T += ni * Math.Pow(Pr, i) * Math.Pow(nu + 1, j);
    }
    return T;
  }

  // Backward equation for region 1, T=f(P,s)
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //   s : float
  //     Specific entropy, [kJ/kgK]
  //
  //   Returns
  //   -------
  //   T : float
  //     Temperature, [K]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  //   Thermodynamic Properties of Water and Steam August 2007,
  //   http://www.iapws.org/relguide/IF97-Rev.html, Eq 13
  //
  //   Examples
  //   --------
  //   >>> _Backward1_T_Ps(3,0.5)
  //   307.842258
  //   >>> _Backward1_T_Ps(80,3)
  //   565.899909
  //

  public static object _Backward1_T_Ps(object P, object s)
  {
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
      3,
      3,
      4
    };
    var J = new List<int>
    {
      0,
      1,
      2,
      3,
      11,
      31,
      0,
      1,
      2,
      3,
      12,
      31,
      0,
      1,
      2,
      9,
      31,
      10,
      32,
      32
    };
    var n = new List<double>
    {
      174.78268058307,
      34.806930892873,
      6.5292584978455,
      0.33039981775489,
      -1.9281382923196E-07,
      -2.4909197244573E-23,
      -0.26107636489332,
      0.22592965981586,
      -0.064256463395226,
      0.0078876289270526,
      3.5672110607366E-10,
      1.7332496994895E-24,
      0.00056608900654837,
      -0.00032635483139717,
      4.4778286690632E-05,
      -5.1322156908507E-10,
      -4.2522657042207E-26,
      2.6400441360689E-13,
      7.8124600459723E-29,
      -3.0732199903668E-31
    };
    var Pr = P / 1;
    var sigma = s / 1;
    var T = 0;
    foreach (var _tup_1 in zip(I, J, n))
    {
      var i = _tup_1.Item1;
      var j = _tup_1.Item2;
      var ni = _tup_1.Item3;
      T += ni * Math.Pow(Pr, i) * Math.Pow(sigma + 2, j);
    }
    return T;
  }

  // Backward equation for region 1, P=f(h,s)
  //
  //   Parameters
  //   ----------
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //   s : float
  //     Specific entropy, [kJ/kgK]
  //
  //   Returns
  //   -------
  //   P : float
  //     Pressure, [MPa]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Supplementary Release on Backward Equations for Pressure
  //   as a Function of Enthalpy and Entropy p(h,s) for Regions 1 and 2 of the
  //   IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of
  //   Water and Steam, http://www.iapws.org/relguide/Supp-PHS12-2014.pdf, Eq 1
  //
  //   Examples
  //   --------
  //   >>> _Backward1_P_hs(0.001,0)
  //   0.0009800980612
  //   >>> _Backward1_P_hs(90,0)
  //   91.92954727
  //   >>> _Backward1_P_hs(1500,3.4)
  //   58.68294423
  //

  public static object _Backward1_P_hs(object h, object s)
  {
    var I = new List<int>
    {
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
      2,
      2,
      2,
      3,
      4,
      4,
      5
    };
    var J = new List<int>
    {
      0,
      1,
      2,
      4,
      5,
      6,
      8,
      14,
      0,
      1,
      4,
      6,
      0,
      1,
      10,
      4,
      1,
      4,
      0
    };
    var n = new List<double>
    {
      -0.691997014660582,
      -18.361254878756,
      -9.28332409297335,
      65.9639569909906,
      -16.2060388912024,
      450.620017338667,
      854.68067822417,
      6075.23214001162,
      32.6487682621856,
      -26.9408844582931,
      -319.9478483343,
      -928.35430704332,
      30.3634537455249,
      -65.0540422444146,
      -4309.9131651613,
      -747.512324096068,
      730.000345529245,
      1142.84032569021,
      -436.407041874559
    };
    var nu = h / 3400;
    var sigma = s / 7.6;
    var P = 0;
    foreach (var _tup_1 in zip(I, J, n))
    {
      var i = _tup_1.Item1;
      var j = _tup_1.Item2;
      var ni = _tup_1.Item3;
      P += ni * Math.Pow(nu + 0.05, i) * Math.Pow(sigma + 0.05, j);
    }
    return 100 * P;
  }

  // Region 2
  // Basic equation for region 2
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
  //     Dict with calculated properties. The available properties are:
  //
  //       * v: Specific volume, [m³/kg]
  //       * h: Specific enthalpy, [kJ/kg]
  //       * s: Specific entropy, [kJ/kgK]
  //       * cp: Specific isobaric heat capacity, [kJ/kgK]
  //       * cv: Specific isocoric heat capacity, [kJ/kgK]
  //       * w: Speed of sound, [m/s]
  //       * alfav: Cubic expansion coefficient, [1/K]
  //       * kt: Isothermal compressibility, [1/MPa]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  //   Thermodynamic Properties of Water and Steam August 2007,
  //   http://www.iapws.org/relguide/IF97-Rev.html, Eq 15-17
  //
  //   Examples
  //   --------
  //   >>> _Region2(700,30)["v"]
  //   0.00542946619
  //   >>> _Region2(700,30)["h"]
  //   2631.49474
  //   >>> _Region2(700,30)["h"]-30000*_Region2(700,30)["v"]
  //   2468.61076
  //   >>> _Region2(700,0.0035)["s"]
  //   10.1749996
  //   >>> _Region2(700,0.0035)["cp"]
  //   2.08141274
  //   >>> _Region2(700,0.0035)["cv"]
  //   1.61978333
  //   >>> _Region2(300,0.0035)["w"]
  //   427.920172
  //   >>> _Region2(300,0.0035)["alfav"]
  //   0.00337578289
  //   >>> _Region2(300,0.0035)["kt"]
  //   286.239651
  //

  public static Dictionary<string, double> _Region2(double T, double P)
  {
    if(P < 0) { P = Pmin; }

    var Tr = 540 / T;
    var Pr = P / 1;
    var _tup_1 = Region2_cp0(Tr, Pr);
    var go = _tup_1.Item1;
    var gop = _tup_1.Item2;
    var gopp = _tup_1.Item3;
    var got = _tup_1.Item4;
    var gott = _tup_1.Item5;
    var gopt = _tup_1.Item6;
    var Ir = new List<int>
    {
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
      3,
      3,
      3,
      4,
      4,
      4,
      5,
      6,
      6,
      6,
      7,
      7,
      7,
      8,
      8,
      9,
      10,
      10,
      10,
      16,
      16,
      18,
      20,
      20,
      20,
      21,
      22,
      23,
      24,
      24,
      24
    };
    var Jr = new List<int>
    {
      0,
      1,
      2,
      3,
      6,
      1,
      2,
      4,
      7,
      36,
      0,
      1,
      3,
      6,
      35,
      1,
      2,
      3,
      7,
      3,
      16,
      35,
      0,
      11,
      25,
      8,
      36,
      13,
      4,
      10,
      14,
      29,
      50,
      57,
      20,
      35,
      48,
      21,
      53,
      39,
      26,
      40,
      58
    };
    var nr = new List<double>
    {
      -0.0017731742473213,
      -0.017834862292358,
      -0.045996013696365,
      -0.057581259083432,
      -0.05032527872793,
      -3.3032641670203E-05,
      -0.00018948987516315,
      -0.0039392777243355,
      -0.043797295650573,
      -2.6674547914087E-05,
      2.0481737692309E-08,
      4.3870667284435E-07,
      -3.227767723857E-05,
      -0.0015033924542148,
      -0.040668253562649,
      -7.8847309559367E-10,
      1.2790717852285E-08,
      4.8225372718507E-07,
      2.2922076337661E-06,
      -1.6714766451061E-11,
      -0.0021171472321355,
      -23.895741934104,
      -5.905956432427E-18,
      -1.2621808899101E-06,
      -0.038946842435739,
      1.1256211360459E-11,
      -8.2311340897998,
      1.9809712802088E-08,
      1.0406965210174E-19,
      -1.0234747095929E-13,
      -1.0018179379511E-09,
      -8.0882908646985E-11,
      0.10693031879409,
      -0.33662250574171,
      8.9185845355421E-25,
      3.0629316876232E-13,
      -4.2002467698208E-06,
      -5.9056029685639E-26,
      3.7826947613457E-06,
      -1.2768608934681E-15,
      7.3087610595061E-29,
      5.5414715350778E-17,
      -9.436970724121E-07
    };
    var gr = 0;
    foreach (var _tup_2 in zip(Ir, Jr, nr))
    {
      var i = _tup_2.Item1;
      var j = _tup_2.Item2;
      var ni = _tup_2.Item3;
      gr += ni * Math.Pow(Pr, i) * Math.Pow(Tr - 0.5, j);
      grp += ni * i * Math.Pow(Pr, i - 1) * Math.Pow(Tr - 0.5, j);
      grpp += ni * i * (i - 1) * Math.Pow(Pr, i - 2) * Math.Pow(Tr - 0.5, j);
      grt += ni * j * Math.Pow(Pr, i) * Math.Pow(Tr - 0.5, j - 1);
      grtt += ni * j * (j - 1) * Math.Pow(Pr, i) * Math.Pow(Tr - 0.5, j - 2);
      grpt += ni * i * j * Math.Pow(Pr, i - 1) * Math.Pow(Tr - 0.5, j - 1);
    }
    var propiedades = new Dictionary<string, double> { };
    propiedades["T"] = T;
    propiedades["P"] = P;
    propiedades["v"] = Pr * (gop + grp) * R * T / P / 1000;
    propiedades["h"] = Tr * (got + grt) * R * T;
    propiedades["s"] = R * (Tr * (got + grt) - (go + gr));
    propiedades["cp"] = -R * Math.Pow(Tr, 2) * (gott + grtt);
    propiedades["cv"] = R * (-Math.Pow(Tr, 2) * (gott + grtt) - Math.Pow(1 + Pr * grp - Tr * Pr * grpt, 2) / (1 - Math.Pow(Pr, 2) * grpp));
    propiedades["w"] = Math.Pow(R * T * 1000 * (1 + 2 * Pr * grp + Math.Pow(Pr, 2) * Math.Pow(grp, 2)) / (1 - Math.Pow(Pr, 2) * grpp + Math.Pow(1 + Pr * grp - Tr * Pr * grpt, 2) / Math.Pow(Tr, 2) / (gott + grtt)), 0.5);
    propiedades["alfav"] = (1 + Pr * grp - Tr * Pr * grpt) / (1 + Pr * grp) / T;
    propiedades["kt"] = (1 - Math.Pow(Pr, 2) * grpp) / (1 + Pr * grp) / P;
    propiedades["region"] = 2;
    propiedades["x"] = 1;
    return propiedades;
  }

  // Ideal properties for Region 2
  //
  //   Parameters
  //   ----------
  //   Tr : float
  //     Reduced temperature, [-]
  //   Pr : float
  //     Reduced pressure, [-]
  //
  //   Returns
  //   -------
  //   prop : array
  //     Array with ideal Gibbs energy partial derivatives:
  //
  //       * g: Ideal Specific Gibbs energy [kJ/kg]
  //       * gp: ∂g/∂P|T
  //       * gpp: ∂²g/∂P²|T
  //       * gt: ∂g/∂T|P
  //       * gtt: ∂²g/∂T²|P
  //       * gpt: ∂²g/∂T∂P
  //
  //   References
  //   ----------
  //   IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  //   Thermodynamic Properties of Water and Steam August 2007,
  //   http://www.iapws.org/relguide/IF97-Rev.html, Eq 16
  //
  //

  public static object Region2_cp0(object Tr, object Pr)
  {
    var Jo = new List<int> { 0, 1, -5, -4, -3, -2, -1, 2, 3 };
    var no = new List<double> { -9.6927686500217, 10.086655968018, -0.005608791128302, 0.071452738081455, -0.40710498223928, 1.4240819171444, -4.383951131945, -0.28408632460772, 0.021268463753307 };

    var go = Math.Log(Pr);
    var gop = Math.Pow(Pr, -1);
    var gopp = -Math.Pow(Pr, -2);
    var got = 0;

    foreach (var _tup_1 in zip(Jo, no))
    {
      var j = _tup_1.Item1;
      var ni = _tup_1.Item2;
      go += ni * Math.Pow(Tr, j);
      got += ni * j * Math.Pow(Tr, j - 1);
      gott += ni * j * (j - 1) * Math.Pow(Tr, j - 2);
    }
    return Tuple.Create(go, gop, gopp, got, gott, gopt);
  }

  // Define the boundary between Region 2b and 2c, P=f(h)
  //
  //   Parameters
  //   ----------
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //
  //   Returns
  //   -------
  //   P : float
  //     Pressure, [MPa]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  //   Thermodynamic Properties of Water and Steam August 2007,
  //   http://www.iapws.org/relguide/IF97-Rev.html, Eq 20
  //
  //   Examples
  //   --------
  //   >>> _P_2bc(3516.004323)
  //   100.0
  //

  public static object _P_2bc(object h)
  {
    return 905.84278514723 - 0.67955786399241 * h + 0.00012809002730136 * Math.Pow(h, 2);
  }

  // Define the boundary between Region 2b and 2c, h=f(P)
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //
  //   Returns
  //   -------
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  //   Thermodynamic Properties of Water and Steam August 2007,
  //   http://www.iapws.org/relguide/IF97-Rev.html, Eq 21
  //
  //   Examples
  //   --------
  //   >>> _hbc_P(100)
  //   3516.004323
  //

  public static object _hbc_P(object P)
  {
    return 2652.6571908428 + Math.Pow((P - 4.5257578905948) / 0.00012809002730136, 0.5);
  }

  // Define the boundary between Region 2a and 2b, h=f(s)
  //
  //   Parameters
  //   ----------
  //   s : float
  //     Specific entropy, [kJ/kgK]
  //
  //   Returns
  //   -------
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Supplementary Release on Backward Equations for Pressure
  //   as a Function of Enthalpy and Entropy p(h,s) for Regions 1 and 2 of the
  //   IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of
  //   Water and Steam, http://www.iapws.org/relguide/Supp-PHS12-2014.pdf, Eq 2
  //
  //   Examples
  //   --------
  //   >>> _hab_s(7)
  //   3376.437884
  //

  public static object _hab_s(object s)
  {
    var smin = _Region2(_TSat_P(4), 4)["s"];
    var smax = _Region2(1073.15, 4)["s"];
    if(s < smin) { var h = 0; }
    else if(s > smax) { h = 5000; }
    else { h = -3498.98083432139 + 2575.60716905876 * s - 421.073558227969 * Math.Pow(s, 2) + 27.6349063799944 * Math.Pow(s, 3); }
    return h;
  }

  // Backward equation for region 2a, T=f(P,h)
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //
  //   Returns
  //   -------
  //   T : float
  //     Temperature, [K]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  //   Thermodynamic Properties of Water and Steam August 2007,
  //   http://www.iapws.org/relguide/IF97-Rev.html, Eq 22
  //
  //   Examples
  //   --------
  //   >>> _Backward2a_T_Ph(0.001,3000)
  //   534.433241
  //   >>> _Backward2a_T_Ph(3,4000)
  //   1010.77577
  //

  public static object _Backward2a_T_Ph(object P, object h)
  {
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
      3,
      3,
      4,
      4,
      4,
      5,
      5,
      5,
      6,
      6,
      7
    };
    var J = new List<int>
    {
      0,
      1,
      2,
      3,
      7,
      20,
      0,
      1,
      2,
      3,
      7,
      9,
      11,
      18,
      44,
      0,
      2,
      7,
      36,
      38,
      40,
      42,
      44,
      24,
      44,
      12,
      32,
      44,
      32,
      36,
      42,
      34,
      44,
      28
    };
    var n = new List<double>
    {
      1089.8952318288,
      849.51654495535,
      -107.81748091826,
      33.153654801263,
      -7.4232016790248,
      11.765048724356,
      1.844574935579,
      -4.1792700549624,
      6.2478196935812,
      -17.344563108114,
      -200.58176862096,
      271.96065473796,
      -455.11318285818,
      3091.9688604755,
      252266.40357872,
      -0.0061707422868339,
      -0.31078046629583,
      11.670873077107,
      128127984.04046,
      -985549096.23276,
      2822454697.3002,
      -3594897141.0703,
      1722734991.3197,
      -13551.334240775,
      12848734.66465,
      1.3865724283226,
      235988.32556514,
      -13105236.545054,
      7399.9835474766,
      -551966.9703006,
      3715408.5996233,
      19127.72923966,
      -415351.64835634,
      -62.459855192507
    };
    var Pr = P / 1;
    var nu = h / 2000;
    var T = 0;
    foreach (var _tup_1 in zip(I, J, n))
    {
      var i = _tup_1.Item1;
      var j = _tup_1.Item2;
      var ni = _tup_1.Item3;
      T += ni * Math.Pow(Pr, i) * Math.Pow(nu - 2.1, j);
    }
    return T;
  }

  // Backward equation for region 2b, T=f(P,h)
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //
  //   Returns
  //   -------
  //   T : float
  //     Temperature, [K]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  //   Thermodynamic Properties of Water and Steam August 2007,
  //   http://www.iapws.org/relguide/IF97-Rev.html, Eq 23
  //
  //   Examples
  //   --------
  //   >>> _Backward2b_T_Ph(5,4000)
  //   1015.31583
  //   >>> _Backward2b_T_Ph(25,3500)
  //   875.279054
  //

  public static object _Backward2b_T_Ph(object P, object h)
  {
    var I = new List<int>
    {
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
      4,
      5,
      5,
      5,
      6,
      7,
      7,
      9,
      9
    };
    var J = new List<int>
    {
      0,
      1,
      2,
      12,
      18,
      24,
      28,
      40,
      0,
      2,
      6,
      12,
      18,
      24,
      28,
      40,
      2,
      8,
      18,
      40,
      1,
      2,
      12,
      24,
      2,
      12,
      18,
      24,
      28,
      40,
      18,
      24,
      40,
      28,
      2,
      28,
      1,
      40
    };
    var n = new List<double>
    {
      1489.5041079516,
      743.07798314034,
      -97.708318797837,
      2.4742464705674,
      -0.63281320016026,
      1.1385952129658,
      -0.47811863648625,
      0.0085208123431544,
      0.93747147377932,
      3.3593118604916,
      3.3809355601454,
      0.16844539671904,
      0.73875745236695,
      -0.47128737436186,
      0.15020273139707,
      -0.002176411421975,
      -0.021810755324761,
      -0.10829784403677,
      -0.046333324635812,
      7.1280351959551E-05,
      0.00011032831789999,
      0.00018955248387902,
      0.0030891541160537,
      0.0013555504554949,
      2.8640237477456E-07,
      -1.0779857357512E-05,
      -7.6462712454814E-05,
      1.4052392818316E-05,
      -3.1083814331434E-05,
      -1.0302738212103E-06,
      2.821728163504E-07,
      1.2704902271945E-06,
      7.3803353468292E-08,
      -1.1030139238909E-08,
      -8.1456365207833E-14,
      -2.5180545682962E-11,
      -1.7565233969407E-18,
      8.6934156344163E-15
    };
    var Pr = P / 1;
    var nu = h / 2000;
    var T = 0;
    foreach (var _tup_1 in zip(I, J, n))
    {
      var i = _tup_1.Item1;
      var j = _tup_1.Item2;
      var ni = _tup_1.Item3;
      T += ni * Math.Pow(Pr - 2, i) * Math.Pow(nu - 2.6, j);
    }
    return T;
  }

  // Backward equation for region 2c, T=f(P,h)
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //
  //   Returns
  //   -------
  //   T : float
  //     Temperature, [K]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  //   Thermodynamic Properties of Water and Steam August 2007,
  //   http://www.iapws.org/relguide/IF97-Rev.html, Eq 24
  //
  //   Examples
  //   --------
  //   >>> _Backward2c_T_Ph(40,2700)
  //   743.056411
  //   >>> _Backward2c_T_Ph(60,3200)
  //   882.756860
  //

  public static object _Backward2c_T_Ph(object P, object h)
  {
    var I = new List<int>
    {
      -7,
      -7,
      -6,
      -6,
      -5,
      -5,
      -2,
      -2,
      -1,
      -1,
      0,
      0,
      1,
      1,
      2,
      6,
      6,
      6,
      6,
      6,
      6,
      6,
      6
    };
    var J = new List<int>
    {
      0,
      4,
      0,
      2,
      0,
      2,
      0,
      1,
      0,
      2,
      0,
      1,
      4,
      8,
      4,
      0,
      1,
      4,
      10,
      12,
      16,
      20,
      22
    };
    var n = new List<double>
    {
      -3236839855524.2,
      7326335090218.1,
      358250899454.47,
      -583401318515.9,
      -10783068217.47,
      20825544563.171,
      610747.83564516,
      859777.2253558,
      -25745.72360417,
      31081.088422714,
      1208.2315865936,
      482.19755109255,
      3.7966001272486,
      -10.842984880077,
      -0.04536417267666,
      1.4559115658698E-13,
      1.126159740723E-12,
      -1.7804982240686E-11,
      1.2324579690832E-07,
      -1.1606921130984E-06,
      2.7846367088554E-05,
      -0.00059270038474176,
      0.0012918582991878
    };
    var Pr = P / 1;
    var nu = h / 2000;
    var T = 0;
    foreach (var _tup_1 in zip(I, J, n))
    {
      var i = _tup_1.Item1;
      var j = _tup_1.Item2;
      var ni = _tup_1.Item3;
      T += ni * Math.Pow(Pr + 25, i) * Math.Pow(nu - 1.8, j);
    }
    return T;
  }

  // Backward equation for region 2, T=f(P,h)
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //
  //   Returns
  //   -------
  //   T : float
  //     Temperature, [K]
  //

  public static object _Backward2_T_Ph(object P, object h)
  {
    object T;
    if(P <= 4)
    {
      T = _Backward2a_T_Ph(P, h);
    }
    else if(4 < P <= 6.546699678)
    {
      T = _Backward2b_T_Ph(P, h);
    }
    else
    {
      var hf = _hbc_P(P);
      if(h >= hf)
      {
        T = _Backward2b_T_Ph(P, h);
      }
      else
      {
        T = _Backward2c_T_Ph(P, h);
      }
    }
    if(P <= 22.064)
    {
      var Tsat = _TSat_P(P);
      T = max(Tsat, T);
    }
    return T;
  }

  // Backward equation for region 2a, T=f(P,s)
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //   s : float
  //     Specific entropy, [kJ/kgK]
  //
  //   Returns
  //   -------
  //   T : float
  //     Temperature, [K]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  //   Thermodynamic Properties of Water and Steam August 2007,
  //   http://www.iapws.org/relguide/IF97-Rev.html, Eq 25
  //
  //   Examples
  //   --------
  //   >>> _Backward2a_T_Ps(0.1,7.5)
  //   399.517097
  //   >>> _Backward2a_T_Ps(2.5,8)
  //   1039.84917
  //

  public static object _Backward2a_T_Ps(object P, object s)
  {
    var I = new List<double>
    {
      -1.5,
      -1.5,
      -1.5,
      -1.5,
      -1.5,
      -1.5,
      -1.25,
      -1.25,
      -1.25,
      -1.0,
      -1.0,
      -1.0,
      -1.0,
      -1.0,
      -1.0,
      -0.75,
      -0.75,
      -0.5,
      -0.5,
      -0.5,
      -0.5,
      -0.25,
      -0.25,
      -0.25,
      -0.25,
      0.25,
      0.25,
      0.25,
      0.25,
      0.5,
      0.5,
      0.5,
      0.5,
      0.5,
      0.5,
      0.5,
      0.75,
      0.75,
      0.75,
      0.75,
      1.0,
      1.0,
      1.25,
      1.25,
      1.5,
      1.5
    };
    var J = new List<int>
    {
      -24,
      -23,
      -19,
      -13,
      -11,
      -10,
      -19,
      -15,
      -6,
      -26,
      -21,
      -17,
      -16,
      -9,
      -8,
      -15,
      -14,
      -26,
      -13,
      -9,
      -7,
      -27,
      -25,
      -11,
      -6,
      1,
      4,
      8,
      11,
      0,
      1,
      5,
      6,
      10,
      14,
      16,
      0,
      4,
      9,
      17,
      7,
      18,
      3,
      15,
      5,
      18
    };
    var n = new List<double>
    {
      -392359.83861984,
      515265.7382727,
      40482.443161048,
      -321.93790923902,
      96.961424218694,
      -22.867846371773,
      -449429.14124357,
      -5011.8336020166,
      0.35684463560015,
      44235.33584819,
      -13673.388811708,
      421632.60207864,
      22516.925837475,
      474.42144865646,
      -149.31130797647,
      -197811.26320452,
      -23554.39947076,
      -19070.616302076,
      55375.669883164,
      3829.3691437363,
      -603.91860580567,
      1936.3102620331,
      4266.064369861,
      -5978.0638872718,
      -704.01463926862,
      338.36784107553,
      20.862786635187,
      0.033834172656196,
      -4.3124428414893E-05,
      166.53791356412,
      -139.86292055898,
      -0.78849547999872,
      0.072132411753872,
      -0.0059754839398283,
      -1.2141358953904E-05,
      2.3227096733871E-07,
      -10.538463566194,
      2.0718925496502,
      -0.072193155260427,
      2.074988708112E-07,
      -0.018340657911379,
      2.9036272348696E-07,
      0.21037527893619,
      0.00025681239729999,
      -0.012799002933781,
      -8.2198102652018E-06
    };
    var Pr = P / 1;
    var sigma = s / 2;
    var T = 0;
    foreach (var _tup_1 in zip(I, J, n))
    {
      var i = _tup_1.Item1;
      var j = _tup_1.Item2;
      var ni = _tup_1.Item3;
      T += ni * Math.Pow(Pr, i) * Math.Pow(sigma - 2, j);
    }
    return T;
  }

  // Backward equation for region 2b, T=f(P,s)
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //   s : float
  //     Specific entropy, [kJ/kgK]
  //
  //   Returns
  //   -------
  //   T : float
  //     Temperature, [K]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  //   Thermodynamic Properties of Water and Steam August 2007,
  //   http://www.iapws.org/relguide/IF97-Rev.html, Eq 26
  //
  //   Examples
  //   --------
  //   >>> _Backward2b_T_Ps(8,6)
  //   600.484040
  //   >>> _Backward2b_T_Ps(90,6)
  //   1038.01126
  //

  public static object _Backward2b_T_Ps(object P, object s)
  {
    var I = new List<int>
    {
      -6,
      -6,
      -5,
      -5,
      -4,
      -4,
      -4,
      -3,
      -3,
      -3,
      -3,
      -2,
      -2,
      -2,
      -2,
      -1,
      -1,
      -1,
      -1,
      -1,
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
      5,
      5,
      5
    };
    var J = new List<int>
    {
      0,
      11,
      0,
      11,
      0,
      1,
      11,
      0,
      1,
      11,
      12,
      0,
      1,
      6,
      10,
      0,
      1,
      5,
      8,
      9,
      0,
      1,
      2,
      4,
      5,
      6,
      9,
      0,
      1,
      2,
      3,
      7,
      8,
      0,
      1,
      5,
      0,
      1,
      3,
      0,
      1,
      0,
      1,
      2
    };
    var n = new List<double>
    {
      316876.65083497,
      20.864175881858,
      -398593.99803599,
      -21.816058518877,
      223697.85194242,
      -2784.1703445817,
      9.920743607148,
      -75197.512299157,
      2970.8605951158,
      -3.4406878548526,
      0.38815564249115,
      17511.29508575,
      -1423.7112854449,
      1.0943803364167,
      0.89971619308495,
      -3375.9740098958,
      471.62885818355,
      -1.9188241993679,
      0.41078580492196,
      -0.33465378172097,
      1387.0034777505,
      -406.63326195838,
      41.72734715961,
      2.1932549434532,
      -1.0320050009077,
      0.35882943516703,
      0.0052511453726066,
      12.838916450705,
      -2.8642437219381,
      0.56912683664855,
      -0.099962954584931,
      -0.0032632037778459,
      0.00023320922576723,
      -0.1533480985745,
      0.029072288239902,
      0.00037534702741167,
      0.0017296691702411,
      -0.00038556050844504,
      -3.5017712292608E-05,
      -1.4566393631492E-05,
      5.6420857267269E-06,
      4.1286150074605E-08,
      -2.0684671118824E-08,
      1.6409393674725E-09
    };
    var Pr = P / 1;
    var sigma = s / 0.7853;
    var T = 0;
    foreach (var _tup_1 in zip(I, J, n))
    {
      var i = _tup_1.Item1;
      var j = _tup_1.Item2;
      var ni = _tup_1.Item3;
      T += ni * Math.Pow(Pr, i) * Math.Pow(10 - sigma, j);
    }
    return T;
  }

  // Backward equation for region 2c, T=f(P,s)
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //   s : float
  //     Specific entropy, [kJ/kgK]
  //
  //   Returns
  //   -------
  //   T : float
  //     Temperature, [K]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  //   Thermodynamic Properties of Water and Steam August 2007,
  //   http://www.iapws.org/relguide/IF97-Rev.html, Eq 27
  //
  //   Examples
  //   --------
  //   >>> _Backward2c_T_Ps(20,5.75)
  //   697.992849
  //   >>> _Backward2c_T_Ps(80,5.75)
  //   949.017998
  //

  public static object _Backward2c_T_Ps(object P, object s)
  {
    var I = new List<int>
    {
      -2,
      -2,
      -1,
      0,
      0,
      0,
      0,
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
      5,
      6,
      6,
      7,
      7,
      7,
      7,
      7
    };
    var J = new List<int>
    {
      0,
      1,
      0,
      0,
      1,
      2,
      3,
      0,
      1,
      3,
      4,
      0,
      1,
      2,
      0,
      1,
      5,
      0,
      1,
      4,
      0,
      1,
      2,
      0,
      1,
      0,
      1,
      3,
      4,
      5
    };
    var n = new List<double>
    {
      909.68501005365,
      2404.566708842,
      -591.6232638713,
      541.45404128074,
      -270.98308411192,
      979.76525097926,
      -469.66772959435,
      14.399274604723,
      -19.104204230429,
      5.3299167111971,
      -21.252975375934,
      -0.3114733441376,
      0.60334840894623,
      -0.042764839702509,
      0.0058185597255259,
      -0.014597008284753,
      0.0056631175631027,
      -7.6155864584577E-05,
      0.00022440342919332,
      -1.2561095013413E-05,
      6.3323132660934E-07,
      -2.0541989675375E-06,
      3.6405370390082E-08,
      -2.9759897789215E-09,
      1.0136618529763E-08,
      5.9925719692351E-12,
      -2.0677870105164E-11,
      -2.0874278181886E-11,
      1.0162166825089E-10,
      -1.6429828281347E-10
    };
    var Pr = P / 1;
    var sigma = s / 2.9251;
    var T = 0;
    foreach (var _tup_1 in zip(I, J, n))
    {
      var i = _tup_1.Item1;
      var j = _tup_1.Item2;
      var ni = _tup_1.Item3;
      T += ni * Math.Pow(Pr, i) * Math.Pow(2 - sigma, j);
    }
    return T;
  }

  // Backward equation for region 2, T=f(P,s)
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //   s : float
  //     Specific entropy, [kJ/kgK]
  //
  //   Returns
  //   -------
  //   T : float
  //     Temperature, [K]
  //

  public static object _Backward2_T_Ps(object P, object s)
  {
    object T;
    if(P <= 4)
    {
      T = _Backward2a_T_Ps(P, s);
    }
    else if(s >= 5.85)
    {
      T = _Backward2b_T_Ps(P, s);
    }
    else
    {
      T = _Backward2c_T_Ps(P, s);
    }
    if(P <= 22.064)
    {
      var Tsat = _TSat_P(P);
      T = max(Tsat, T);
    }
    return T;
  }

  // Backward equation for region 2a, P=f(h,s)
  //
  //   Parameters
  //   ----------
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //   s : float
  //     Specific entropy, [kJ/kgK]
  //
  //   Returns
  //   -------
  //   P : float
  //     Pressure, [MPa]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Supplementary Release on Backward Equations for Pressure
  //   as a Function of Enthalpy and Entropy p(h,s) for Regions 1 and 2 of the
  //   IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of
  //   Water and Steam, http://www.iapws.org/relguide/Supp-PHS12-2014.pdf, Eq 3
  //
  //   Examples
  //   --------
  //   >>> _Backward2a_P_hs(2800,6.5)
  //   1.371012767
  //   >>> _Backward2a_P_hs(2800,9.5)
  //   0.001879743844
  //   >>> _Backward2a_P_hs(4100,9.5)
  //   0.1024788997
  //

  public static object _Backward2a_P_hs(object h, object s)
  {
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
      3,
      3,
      4,
      5,
      5,
      6,
      7
    };
    var J = new List<int>
    {
      1,
      3,
      6,
      16,
      20,
      22,
      0,
      1,
      2,
      3,
      5,
      6,
      10,
      16,
      20,
      22,
      3,
      16,
      20,
      0,
      2,
      3,
      6,
      16,
      16,
      3,
      16,
      3,
      1
    };
    var n = new List<double>
    {
      -0.0182575361923032,
      -0.125229548799536,
      0.592290437320145,
      6.04769706185122,
      238.624965444474,
      -298.639090222922,
      0.051225081304075,
      -0.437266515606486,
      0.413336902999504,
      -5.16468254574773,
      -5.57014838445711,
      12.8555037824478,
      11.414410895329,
      -119.504225652714,
      -2847.7798596156,
      4317.57846408006,
      1.1289404080265,
      1974.09186206319,
      1516.12444706087,
      0.0141324451421235,
      0.585501282219601,
      -2.97258075863012,
      5.94567314847319,
      -6236.56565798905,
      9659.86235133332,
      6.81500934948134,
      -6332.07286824489,
      -5.5891922446576,
      0.0400645798472063
    };
    var nu = h / 4200;
    var sigma = s / 12;
    var suma = 0;
    foreach (var _tup_1 in zip(I, J, n))
    {
      var i = _tup_1.Item1;
      var j = _tup_1.Item2;
      var ni = _tup_1.Item3;
      suma += ni * Math.Pow(nu - 0.5, i) * Math.Pow(sigma - 1.2, j);
    }
    return 4 * Math.Pow(suma, 4);
  }

  // Backward equation for region 2b, P=f(h,s)
  //
  //   Parameters
  //   ----------
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //   s : float
  //     Specific entropy, [kJ/kgK]
  //
  //   Returns
  //   -------
  //   P : float
  //     Pressure, [MPa]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Supplementary Release on Backward Equations for Pressure
  //   as a Function of Enthalpy and Entropy p(h,s) for Regions 1 and 2 of the
  //   IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of
  //   Water and Steam, http://www.iapws.org/relguide/Supp-PHS12-2014.pdf, Eq 4
  //
  //   Examples
  //   --------
  //   >>> _Backward2b_P_hs(2800,6)
  //   4.793911442
  //   >>> _Backward2b_P_hs(3600,6)
  //   83.95519209
  //   >>> _Backward2b_P_hs(3600,7)
  //   7.527161441
  //

  public static object _Backward2b_P_hs(object h, object s)
  {
    var I = new List<int>
    {
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
      3,
      4,
      4,
      5,
      5,
      6,
      6,
      6,
      7,
      7,
      8,
      8,
      8,
      8,
      12,
      14
    };
    var J = new List<int>
    {
      0,
      1,
      2,
      4,
      8,
      0,
      1,
      2,
      3,
      5,
      12,
      1,
      6,
      18,
      0,
      1,
      7,
      12,
      1,
      16,
      1,
      12,
      1,
      8,
      18,
      1,
      16,
      1,
      3,
      14,
      18,
      10,
      16
    };
    var n = new List<double>
    {
      0.0801496989929495,
      -0.543862807146111,
      0.337455597421283,
      8.9055545115745,
      313.840736431485,
      0.797367065977789,
      -1.2161697355624,
      8.72803386937477,
      -16.9769781757602,
      -186.552827328416,
      95115.9274344237,
      -18.9168510120494,
      -4334.0703719484,
      543212633.012715,
      0.144793408386013,
      128.024559637516,
      -67230.9534071268,
      33697238.0095287,
      -586.63419676272,
      -22140322476.9889,
      1716.06668708389,
      -570817595.806302,
      -3121.09693178482,
      -2078413.8463301,
      3056059461577.86,
      3221.57004314333,
      326810259797.295,
      -1441.04158934487,
      410.694867802691,
      109077066873.024,
      -24796465425889.3,
      1888019068.65134,
      -123651009018773.0
    };
    var nu = h / 4100;
    var sigma = s / 7.9;
    var suma = 0;
    foreach (var _tup_1 in zip(I, J, n))
    {
      var i = _tup_1.Item1;
      var j = _tup_1.Item2;
      var ni = _tup_1.Item3;
      suma += ni * Math.Pow(nu - 0.6, i) * Math.Pow(sigma - 1.01, j);
    }
    return 100 * Math.Pow(suma, 4);
  }

  // Backward equation for region 2c, P=f(h,s)
  //
  //   Parameters
  //   ----------
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //   s : float
  //     Specific entropy, [kJ/kgK]
  //
  //   Returns
  //   -------
  //   P : float
  //     Pressure, [MPa]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Supplementary Release on Backward Equations for Pressure
  //   as a Function of Enthalpy and Entropy p(h,s) for Regions 1 and 2 of the
  //   IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of
  //   Water and Steam, http://www.iapws.org/relguide/Supp-PHS12-2014.pdf, Eq 5
  //
  //   Examples
  //   --------
  //   >>> _Backward2c_P_hs(2800,5.1)
  //   94.39202060
  //   >>> _Backward2c_P_hs(2800,5.8)
  //   8.414574124
  //   >>> _Backward2c_P_hs(3400,5.8)
  //   83.76903879
  //

  public static object _Backward2c_P_hs(object h, object s)
  {
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
      4,
      5,
      5,
      5,
      5,
      6,
      6,
      10,
      12,
      16
    };
    var J = new List<int>
    {
      0,
      1,
      2,
      3,
      4,
      8,
      0,
      2,
      5,
      8,
      14,
      2,
      3,
      7,
      10,
      18,
      0,
      5,
      8,
      16,
      18,
      18,
      1,
      4,
      6,
      14,
      8,
      18,
      7,
      7,
      10
    };
    var n = new List<double>
    {
      0.112225607199012,
      -3.39005953606712,
      -32.0503911730094,
      -197.5973051049,
      -407.693861553446,
      13294.3775222331,
      1.70846839774007,
      37.3694198142245,
      3581.44365815434,
      423014.446424664,
      -751071025.760063,
      52.3446127607898,
      -228.351290812417,
      -960652.417056937,
      -80705929.2526074,
      1626980172256.69,
      0.772465073604171,
      46392.9973837746,
      -13731788.5134128,
      1704703926305.12,
      -25110462818730.8,
      31774883083552.0,
      53.8685623675312,
      -55308.9094625169,
      -1028615.22421405,
      2042494187562.34,
      273918446.626977,
      -2.63963146312685E+15,
      -1078908541.08088,
      -29649262098.0124,
      -1.11754907323424E+15
    };
    var nu = h / 3500;
    var sigma = s / 5.9;
    var suma = 0;
    foreach (var _tup_1 in zip(I, J, n))
    {
      var i = _tup_1.Item1;
      var j = _tup_1.Item2;
      var ni = _tup_1.Item3;
      suma += ni * Math.Pow(nu - 0.7, i) * Math.Pow(sigma - 1.1, j);
    }
    return 100 * Math.Pow(suma, 4);
  }

  // Backward equation for region 2, P=f(h,s)
  //
  //   Parameters
  //   ----------
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //   s : float
  //     Specific entropy, [kJ/kgK]
  //
  //   Returns
  //   -------
  //   P : float
  //     Pressure, [MPa]
  //

  public static object _Backward2_P_hs(object h, object s)
  {
    var sfbc = 5.85;
    var hamin = _hab_s(s);
    if(h <= hamin)
    {
      var P = _Backward2a_P_hs(h, s);
    }
    else if(s >= sfbc)
    {
      P = _Backward2b_P_hs(h, s);
    }
    else
    {
      P = _Backward2c_P_hs(h, s);
    }
    return P;
  }

  // Region 3
  // Basic equation for region 3
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
  //   prop : dict
  //     Dict with calculated properties. The available properties are:
  //
  //       * v: Specific volume, [m³/kg]
  //       * h: Specific enthalpy, [kJ/kg]
  //       * s: Specific entropy, [kJ/kgK]
  //       * cp: Specific isobaric heat capacity, [kJ/kgK]
  //       * cv: Specific isocoric heat capacity, [kJ/kgK]
  //       * w: Speed of sound, [m/s]
  //       * alfav: Cubic expansion coefficient, [1/K]
  //       * kt: Isothermal compressibility, [1/MPa]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  //   Thermodynamic Properties of Water and Steam August 2007,
  //   http://www.iapws.org/relguide/IF97-Rev.html, Eq 28
  //
  //   Examples
  //   --------
  //   >>> _Region3(500,650)["P"]
  //   25.5837018
  //   >>> _Region3(500,650)["h"]
  //   1863.43019
  //   >>> p = _Region3(500, 650)
  //   >>> p["h"]-p["P"]*1000*p["v"]
  //   1812.26279
  //   >>> _Region3(200,650)["s"]
  //   4.85438792
  //   >>> _Region3(200,650)["cp"]
  //   44.6579342
  //   >>> _Region3(200,650)["cv"]
  //   4.04118076
  //   >>> _Region3(200,650)["w"]
  //   383.444594
  //   >>> _Region3(500,750)["alfav"]
  //   0.00441515098
  //   >>> _Region3(500,750)["kt"]
  //   0.00806710817
  //

  public static Dictionary<string, double> _Region3(double rho, double T)
  {
    var I = new List<int>
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
      4,
      4,
      4,
      4,
      5,
      5,
      5,
      6,
      6,
      6,
      7,
      8,
      9,
      9,
      10,
      10,
      11
    };
    var J = new List<int>
    {
      0,
      1,
      2,
      7,
      10,
      12,
      23,
      2,
      6,
      15,
      17,
      0,
      2,
      6,
      7,
      22,
      26,
      0,
      2,
      4,
      16,
      26,
      0,
      2,
      4,
      26,
      1,
      3,
      26,
      0,
      2,
      26,
      2,
      26,
      2,
      26,
      0,
      1,
      26
    };
    var n = new List<double>
    {
      -15.732845290239,
      20.944396974307,
      -7.6867707878716,
      2.6185947787954,
      -2.808078114862,
      1.2053369696517,
      -0.0084566812812502,
      -1.2654315477714,
      -1.1524407806681,
      0.88521043984318,
      -0.64207765181607,
      0.38493460186671,
      -0.85214708824206,
      4.8972281541877,
      -3.0502617256965,
      0.039420536879154,
      0.12558408424308,
      -0.2799932969871,
      1.389979956946,
      -2.018991502357,
      -0.0082147637173963,
      -0.47596035734923,
      0.0439840744735,
      -0.44476435428739,
      0.90572070719733,
      0.70522450087967,
      0.10770512626332,
      -0.32913623258954,
      -0.50871062041158,
      -0.022175400873096,
      0.094260751665092,
      0.16436278447961,
      -0.013503372241348,
      -0.014834345352472,
      0.00057922953628084,
      0.0032308904703711,
      8.0964802996215E-05,
      -0.00016557679795037,
      -4.4923899061815E-05
    };
    var d = rho / rhoc;
    var Tr = Tc / T;
    var g = 1.0658070028513 * Math.Log(d);
    var gd = 1.0658070028513 / d;
    var gdd = -1.0658070028513 / Math.Pow(d, 2);
    var gt = 0;
    foreach (var _tup_1 in zip(I, J, n))
    {
      var i = _tup_1.Item1;
      var j = _tup_1.Item2;
      var ni = _tup_1.Item3;
      g += ni * Math.Pow(d, i) * Math.Pow(Tr, j);
      gd += ni * i * Math.Pow(d, i - 1) * Math.Pow(Tr, j);
      gdd += ni * i * (i - 1) * Math.Pow(d, i - 2) * Math.Pow(Tr, j);
      gt += ni * j * Math.Pow(d, i) * Math.Pow(Tr, j - 1);
      gtt += ni * j * (j - 1) * Math.Pow(d, i) * Math.Pow(Tr, j - 2);
      gdt += ni * i * j * Math.Pow(d, i - 1) * Math.Pow(Tr, j - 1);
    }

    var propiedades = new Dictionary<string, double> { };
    propiedades["T"] = T;
    propiedades["P"] = d * gd * R * T * rho / 1000;
    propiedades["v"] = 1 / rho;
    propiedades["h"] = R * T * (Tr * gt + d * gd);
    propiedades["s"] = R * (Tr * gt - g);
    propiedades["cp"] = R * (-Math.Pow(Tr, 2) * gtt + Math.Pow(d * gd - d * Tr * gdt, 2) / (2 * d * gd + Math.Pow(d, 2) * gdd));
    propiedades["cv"] = -R * Math.Pow(Tr, 2) * gtt;
    propiedades["w"] = Math.Sqrt(R * T * 1000 * (2 * d * gd + Math.Pow(d, 2) * gdd - Math.Pow(d * gd - d * Tr * gdt, 2) / Math.Pow(Tr, 2) / gtt));
    propiedades["alfav"] = (gd - Tr * gdt) / (2 * gd + d * gdd) / T;
    propiedades["kt"] = 1 / (2 * d * gd + Math.Pow(d, 2) * gdd) / rho / R / T * 1000;
    propiedades["region"] = 3;
    propiedades["x"] = 1;
    return propiedades;
  }

  // Define the boundary between Region 3a-3b, h=f(P)
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //
  //   Returns
  //   -------
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //
  //   Examples
  //   --------
  //   >>> _h_3ab(25)
  //   2095.936454
  //

  public static object _h_3ab(object P)
  {
    return 2014.64004206875 + 3.74696550136983 * P - 0.0219921901054187 * Math.Pow(P, 2) + 8.7513168600995E-05 * Math.Pow(P, 3);
  }

  // Define the boundary between Region 3a-3b, T=f(P)
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //
  //   Returns
  //   -------
  //   T : float
  //     Temperature, [K]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Supplementary Release on Backward Equations for Specific
  //   Volume as a Function of Pressure and Temperature v(p,T) for Region 3 of the
  //   IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water
  //   and Steam, http://www.iapws.org/relguide/Supp-VPT3-2016.pdf, Eq. 2
  //
  //   Examples
  //   --------
  //   >>> _tab_P(40)
  //   693.0341408
  //

  public static object _tab_P(object P)
  {
    var I = new List<int>
    {
      0,
      1,
      2,
      -1,
      -2
    };
    var n = new List<double>
    {
      1547.93642129415,
      -187.661219490113,
      21.3144632222113,
      -1918.87498864292,
      918.419702359447
    };
    var Pr = P / 1;
    var T = 0;
    foreach (var _tup_1 in zip(I, n))
    {
      var i = _tup_1.Item1;
      var ni = _tup_1.Item2;
      T += ni * Math.Pow(Math.Log(Pr), i);
    }
    return T;
  }

  // Define the boundary between Region 3o-3p, T=f(P)
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //
  //   Returns
  //   -------
  //   T : float
  //     Temperature, [K]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Supplementary Release on Backward Equations for Specific
  //   Volume as a Function of Pressure and Temperature v(p,T) for Region 3 of the
  //   IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water
  //   and Steam, http://www.iapws.org/relguide/Supp-VPT3-2016.pdf, Eq. 2
  //
  //   Examples
  //   --------
  //   >>> _top_P(22.8)
  //   650.0106943
  //

  public static object _top_P(object P)
  {
    var I = new List<int>
    {
      0,
      1,
      2,
      -1,
      -2
    };
    var n = new List<double>
    {
      969.461372400213,
      -332.500170441278,
      64.2859598466067,
      773.845935768222,
      -1523.13732937084
    };
    var Pr = P / 1;
    var T = 0;
    foreach (var _tup_1 in zip(I, n))
    {
      var i = _tup_1.Item1;
      var ni = _tup_1.Item2;
      T += ni * Math.Pow(Math.Log(Pr), i);
    }
    return T;
  }

  // Define the boundary between Region 3w-3x, T=f(P)
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //
  //   Returns
  //   -------
  //   T : float
  //     Temperature, [K]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Supplementary Release on Backward Equations for Specific
  //   Volume as a Function of Pressure and Temperature v(p,T) for Region 3 of the
  //   IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water
  //   and Steam, http://www.iapws.org/relguide/Supp-VPT3-2016.pdf, Eq. 2
  //
  //   Examples
  //   --------
  //   >>> _twx_P(22.3)
  //   648.2049480
  //

  public static object _twx_P(object P)
  {
    var I = new List<int>
    {
      0,
      1,
      2,
      -1,
      -2
    };
    var n = new List<double>
    {
      7.2805260914538,
      97.3505869861952,
      14.7370491183191,
      329.196213998375,
      873.371668682417
    };
    var Pr = P / 1;
    var T = 0;
    foreach (var _tup_1 in zip(I, n))
    {
      var i = _tup_1.Item1;
      var ni = _tup_1.Item2;
      T += ni * Math.Pow(Math.Log(Pr), i);
    }
    return T;
  }

  // Define the boundary between Region 3e-3f, T=f(P)
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //
  //   Returns
  //   -------
  //   T : float
  //     Temperature, [K]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Supplementary Release on Backward Equations for Specific
  //   Volume as a Function of Pressure and Temperature v(p,T) for Region 3 of the
  //   IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water
  //   and Steam, http://www.iapws.org/relguide/Supp-VPT3-2016.pdf, Eq. 3
  //
  //   Examples
  //   --------
  //   >>> _tef_P(40)
  //   713.9593992
  //

  public static object _tef_P(object P)
  {
    return 3.727888004 * (P - 22.064) + 647.096;
  }

  // Define the boundary between 3x-3y, T=f(P)
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //   xy: string
  //     Subregions options: cd, gh, ij, jk, mn, qu, rx, uv
  //
  //   Returns
  //   -------
  //   T : float
  //     Temperature, [K]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Supplementary Release on Backward Equations for Specific
  //   Volume as a Function of Pressure and Temperature v(p,T) for Region 3 of the
  //   IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water
  //   and Steam, http://www.iapws.org/relguide/Supp-VPT3-2016.pdf, Eq. 1
  //
  //   Examples
  //   --------
  //   >>> _txx_P(25,"cd")
  //   649.3659208
  //   >>> _txx_P(23,"gh")
  //   649.8873759
  //   >>> _txx_P(23,"ij")
  //   651.5778091
  //   >>> _txx_P(23,"jk")
  //   655.8338344
  //   >>> _txx_P(22.8,"mn")
  //   649.6054133
  //   >>> _txx_P(22,"qu")
  //   645.6355027
  //   >>> _txx_P(22,"rx")
  //   648.2622754
  //   >>> _txx_P(22.3,"uv")
  //   647.7996121
  //

  public static object _txx_P(object P, object xy)
  {
    var ng = new Dictionary<object, object>
    {
      { "cd", new List<double> { 585.276966696349, 2.78233532206915, -0.0127283549295878, 0.000159090746562729 }},
      { "gh", new List<double> { -24928.4240900418, 4281.43584791546, -269.02917314013, 7.51608051114157, -0.0787105249910383 }},
      { "ij", new List<double> { 584.814781649163, -0.616179320924617, 0.260763050899562, -0.00587071076864459, 5.15308185433082E-05 }},
      { "jk", new List<double> { 617.229772068439, -7.70600270141675, 0.697072596851896, -0.0157391839848015, 0.000137897492684194 }},
      { "mn", new List<double> { 535.339483742384, 7.61978122720128, -0.158365725441648, 0.00192871054508108 }},
      { "qu", new List<double> { 565.603648239126, 5.29062258221222, -0.102020639611016, 0.00122240301070145 }},
      { "rx", new List<double> { 584.561202520006, -1.02961025163669, 0.243293362700452, -0.00294905044740799 }},
      { "uv", new List<double> { 528.199646263062, 8.90579602135307, -0.222814134903755, 0.00286791682263697 }}
    };
    var n = ng[xy];
    var Pr = P / 1;
    var T = 0;
    foreach (var _tup_1 in n.Select((_p_1,_p_2) => Tuple.Create(_p_2, _p_1)))
    {
      var i = _tup_1.Item1;
      var ni = _tup_1.Item2;
      T += ni * Math.Pow(Pr, i);
    }
    return T;
  }

  // Backward equation for region 3a, v=f(P,h)
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Supplementary Release on Backward Equations for the
  //   Functions T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS
  //   Industrial Formulation 1997 for the Thermodynamic Properties of Water and
  //   Steam, http://www.iapws.org/relguide/Supp-Tv%28ph,ps%293-2014.pdf, Eq 4
  //
  //   Returns
  //   -------
  //   v : float
  //     Specific volume, [m³/kg]
  //
  //   Examples
  //   --------
  //   >>> _Backward3a_v_Ph(20,1700)
  //   0.001749903962
  //   >>> _Backward3a_v_Ph(100,2100)
  //   0.001676229776
  //

  public static object _Backward3a_v_Ph(object P, object h)
  {
    var I = new List<int>
    {
      -12,
      -12,
      -12,
      -12,
      -10,
      -10,
      -10,
      -8,
      -8,
      -6,
      -6,
      -6,
      -4,
      -4,
      -3,
      -2,
      -2,
      -1,
      -1,
      -1,
      -1,
      0,
      0,
      1,
      1,
      1,
      2,
      2,
      3,
      4,
      5,
      8
    };
    var J = new List<int>
    {
      6,
      8,
      12,
      18,
      4,
      7,
      10,
      5,
      12,
      3,
      4,
      22,
      2,
      3,
      7,
      3,
      16,
      0,
      1,
      2,
      3,
      0,
      1,
      0,
      1,
      2,
      0,
      2,
      0,
      2,
      2,
      2
    };
    var n = new List<double>
    {
      0.00529944062966028,
      -0.170099690234461,
      11.1323814312927,
      -2178.98123145125,
      -0.000506061827980875,
      0.556495239685324,
      -9.43672726094016,
      -0.297856807561527,
      93.9353943717186,
      0.0192944939465981,
      0.421740664704763,
      -3689141.2628233,
      -0.00737566847600639,
      -0.354753242424366,
      -1.99768169338727,
      1.15456297059049,
      5683.6687581596,
      0.00808169540124668,
      0.172416341519307,
      1.04270175292927,
      -0.297691372792847,
      0.560394465163593,
      0.275234661176914,
      -0.148347894866012,
      -0.0651142513478515,
      -2.92468715386302,
      0.0664876096952665,
      3.52335014263844,
      -0.0146340792313332,
      -2.24503486668184,
      1.10533464706142,
      -0.0408757344495612
    };
    var Pr = P / 100;
    var nu = h / 2100;
    var suma = 0;
    foreach (var _tup_1 in zip(I, J, n))
    {
      var i = _tup_1.Item1;
      var j = _tup_1.Item2;
      var ni = _tup_1.Item3;
      suma += ni * Math.Pow(Pr + 0.128, i) * Math.Pow(nu - 0.727, j);
    }
    return 0.0028 * suma;
  }

  // Backward equation for region 3b, v=f(P,h)
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //
  //   Returns
  //   -------
  //   v : float
  //     Specific volume, [m³/kg]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Supplementary Release on Backward Equations for the
  //   Functions T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS
  //   Industrial Formulation 1997 for the Thermodynamic Properties of Water and
  //   Steam, http://www.iapws.org/relguide/Supp-Tv%28ph,ps%293-2014.pdf, Eq 5
  //
  //   Examples
  //   --------
  //   >>> _Backward3b_v_Ph(20,2500)
  //   0.006670547043
  //   >>> _Backward3b_v_Ph(100,2700)
  //   0.002404234998
  //

  public static object _Backward3b_v_Ph(object P, object h)
  {
    var I = new List<int>
    {
      -12,
      -12,
      -8,
      -8,
      -8,
      -8,
      -8,
      -8,
      -6,
      -6,
      -6,
      -6,
      -6,
      -6,
      -4,
      -4,
      -4,
      -3,
      -3,
      -2,
      -2,
      -1,
      -1,
      -1,
      -1,
      0,
      1,
      1,
      2,
      2
    };
    var J = new List<int>
    {
      0,
      1,
      0,
      1,
      3,
      6,
      7,
      8,
      0,
      1,
      2,
      5,
      6,
      10,
      3,
      6,
      10,
      0,
      2,
      1,
      2,
      0,
      1,
      4,
      5,
      0,
      0,
      1,
      2,
      6
    };
    var n = new List<double>
    {
      -2.25196934336318E-09,
      1.40674363313486E-08,
      2.3378408528056E-06,
      -3.31833715229001E-05,
      0.00107956778514318,
      -0.271382067378863,
      1.07202262490333,
      -0.853821329075382,
      -2.15214194340526E-05,
      0.00076965608822273,
      -0.00431136580433864,
      0.453342167309331,
      -0.507749535873652,
      -100.475154528389,
      -0.219201924648793,
      -3.21087965668917,
      607.567815637771,
      0.000557686450685932,
      0.18749904002955,
      0.00905368030448107,
      0.285417173048685,
      0.0329924030996098,
      0.239897419685483,
      4.82754995951394,
      -11.8035753702231,
      0.169490044091791,
      -0.0179967222507787,
      0.0371810116332674,
      -0.0536288335065096,
      1.6069710109252
    };
    var Pr = P / 100;
    var nu = h / 2800;
    var suma = 0;
    foreach (var _tup_1 in zip(I, J, n))
    {
      var i = _tup_1.Item1;
      var j = _tup_1.Item2;
      var ni = _tup_1.Item3;
      suma += ni * Math.Pow(Pr + 0.0661, i) * Math.Pow(nu - 0.72, j);
    }
    return 0.0088 * suma;
  }

  // Backward equation for region 3, v=f(P,h)
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //
  //   Returns
  //   -------
  //   v : float
  //     Specific volume, [m³/kg]
  //

  public static object _Backward3_v_Ph(object P, object h)
  {
    var hf = _h_3ab(P);
    if(h <= hf) { return _Backward3a_v_Ph(P, h); }
    else { return _Backward3b_v_Ph(P, h); }
  }

  // Backward equation for region 3a, T=f(P,h)
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //
  //   Returns
  //   -------
  //   T : float
  //     Temperature, [K]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Supplementary Release on Backward Equations for the
  //   Functions T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS
  //   Industrial Formulation 1997 for the Thermodynamic Properties of Water and
  //   Steam, http://www.iapws.org/relguide/Supp-Tv%28ph,ps%293-2014.pdf, Eq 2
  //
  //   Examples
  //   --------
  //   >>> _Backward3a_T_Ph(20,1700)
  //   629.3083892
  //   >>> _Backward3a_T_Ph(100,2100)
  //   733.6163014
  //

  public static object _Backward3a_T_Ph(object P, object h)
  {
    var I = new List<int>
    {
      -12,
      -12,
      -12,
      -12,
      -12,
      -12,
      -12,
      -12,
      -10,
      -10,
      -10,
      -8,
      -8,
      -8,
      -8,
      -5,
      -3,
      -2,
      -2,
      -2,
      -1,
      -1,
      0,
      0,
      1,
      3,
      3,
      4,
      4,
      10,
      12
    };
    var J = new List<int>
    {
      0,
      1,
      2,
      6,
      14,
      16,
      20,
      22,
      1,
      5,
      12,
      0,
      2,
      4,
      10,
      2,
      0,
      1,
      3,
      4,
      0,
      2,
      0,
      1,
      1,
      0,
      1,
      0,
      3,
      4,
      5
    };
    var n = new List<double>
    {
      -1.33645667811215E-07,
      4.55912656802978E-06,
      -1.46294640700979E-05,
      0.0063934131297008,
      372.783927268847,
      -7186.54377460447,
      573494.7521034,
      -2675693.29111439,
      -3.34066283302614E-05,
      -0.0245479214069597,
      47.8087847764996,
      7.64664131818904E-06,
      0.00128350627676972,
      0.0171219081377331,
      -8.51007304583213,
      -0.0136513461629781,
      -3.84460997596657E-06,
      0.00337423807911655,
      -0.551624873066791,
      0.72920227710747,
      -0.00992522757376041,
      -0.119308831407288,
      0.793929190615421,
      0.454270731799386,
      0.20999859125991,
      -0.00642109823904738,
      -0.023515586860454,
      0.00252233108341612,
      -0.00764885133368119,
      0.0136176427574291,
      -0.0133027883575669
    };
    var Pr = P / 100.0;
    var nu = h / 2300.0;
    var suma = 0;
    foreach (var _tup_1 in zip(I, J, n))
    {
      var i = _tup_1.Item1;
      var j = _tup_1.Item2;
      n = _tup_1.Item3;
      suma += n * Math.Pow(Pr + 0.24, i) * Math.Pow(nu - 0.615, j);
    }
    return 760 * suma;
  }

  // Backward equation for region 3b, T=f(P,h)
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //
  //   Returns
  //   -------
  //   T : float
  //     Temperature, [K]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Supplementary Release on Backward Equations for the
  //   Functions T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS
  //   Industrial Formulation 1997 for the Thermodynamic Properties of Water and
  //   Steam, http://www.iapws.org/relguide/Supp-Tv%28ph,ps%293-2014.pdf, Eq 3
  //
  //   Examples
  //   --------
  //   >>> _Backward3b_T_Ph(20,2500)
  //   641.8418053
  //   >>> _Backward3b_T_Ph(100,2700)
  //   842.0460876
  //

  public static object _Backward3b_T_Ph(object P, object h)
  {
    var I = new List<int>
    {
      -12,
      -12,
      -10,
      -10,
      -10,
      -10,
      -10,
      -8,
      -8,
      -8,
      -8,
      -8,
      -6,
      -6,
      -6,
      -4,
      -4,
      -3,
      -2,
      -2,
      -1,
      -1,
      -1,
      -1,
      -1,
      -1,
      0,
      0,
      1,
      3,
      5,
      6,
      8
    };
    var J = new List<int>
    {
      0,
      1,
      0,
      1,
      5,
      10,
      12,
      0,
      1,
      2,
      4,
      10,
      0,
      1,
      2,
      0,
      1,
      5,
      0,
      4,
      2,
      4,
      6,
      10,
      14,
      16,
      0,
      2,
      1,
      1,
      1,
      1,
      1
    };
    var n = new List<double>
    {
      3.2325457364492E-05,
      -0.000127575556587181,
      -0.000475851877356068,
      0.00156183014181602,
      0.105724860113781,
      -85.8514221132534,
      724.140095480911,
      0.00296475810273257,
      -0.00592721983365988,
      -0.0126305422818666,
      -0.115716196364853,
      84.9000969739595,
      -0.0108602260086615,
      0.0154304475328851,
      0.0750455441524466,
      0.0252520973612982,
      -0.0602507901232996,
      -3.07622221350501,
      -0.0574011959864879,
      5.03471360939849,
      -0.925081888584834,
      3.91733882917546,
      -77.314600713019,
      9493.08762098587,
      -1410437.19679409,
      8491662.30819026,
      0.861095729446704,
      0.32334644281172,
      0.873281936020439,
      -0.436653048526683,
      0.286596714529479,
      -0.131778331276228,
      0.00676682064330275
    };
    var Pr = P / 100.0;
    var nu = h / 2800.0;
    var suma = 0;
    foreach (var _tup_1 in zip(I, J, n))
    {
      var i = _tup_1.Item1;
      var j = _tup_1.Item2;
      n = _tup_1.Item3;
      suma += n * Math.Pow(Pr + 0.298, i) * Math.Pow(nu - 0.72, j);
    }
    return 860 * suma;
  }

  // Backward equation for region 3, T=f(P,h)
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //
  //   Returns
  //   -------
  //   T : float
  //     Temperature, [K]
  //

  public static object _Backward3_T_Ph(object P, object h)
  {
    object T;
    var hf = _h_3ab(P);
    if(h <= hf) { T = _Backward3a_T_Ph(P, h); }
    else { T = _Backward3b_T_Ph(P, h); }
    return T;
  }

  // Backward equation for region 3a, v=f(P,s)
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //   s : float
  //     Specific entropy, [kJ/kgK]
  //
  //   Returns
  //   -------
  //   v : float
  //     Specific volume, [m³/kg]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Supplementary Release on Backward Equations for the
  //   Functions T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS
  //   Industrial Formulation 1997 for the Thermodynamic Properties of Water and
  //   Steam, http://www.iapws.org/relguide/Supp-Tv%28ph,ps%293-2014.pdf, Eq 8
  //
  //   Examples
  //   --------
  //   >>> _Backward3a_v_Ps(20,3.8)
  //   0.001733791463
  //   >>> _Backward3a_v_Ps(100,4)
  //   0.001555893131
  //

  public static object _Backward3a_v_Ps(object P, object s)
  {
    var I = new List<int>
    {
      -12,
      -12,
      -12,
      -10,
      -10,
      -10,
      -10,
      -8,
      -8,
      -8,
      -8,
      -6,
      -5,
      -4,
      -3,
      -3,
      -2,
      -2,
      -1,
      -1,
      0,
      0,
      0,
      1,
      2,
      4,
      5,
      6
    };
    var J = new List<int>
    {
      10,
      12,
      14,
      4,
      8,
      10,
      20,
      5,
      6,
      14,
      16,
      28,
      1,
      5,
      2,
      4,
      3,
      8,
      1,
      2,
      0,
      1,
      3,
      0,
      0,
      2,
      2,
      0
    };
    var n = new List<double>
    {
      79.5544074093975,
      -2382.6124298459,
      17681.3100617787,
      -0.00110524727080379,
      -15.3213833655326,
      297.544599376982,
      -35031520.6871242,
      0.277513761062119,
      -0.523964271036888,
      -148011.182995403,
      1600148.99374266,
      1708023226634.27,
      0.000246866996006494,
      1.6532608479798,
      -0.118008384666987,
      2.537986423559,
      0.965127704669424,
      -28.2172420532826,
      0.203224612353823,
      1.10648186063513,
      0.52612794845128,
      0.277000018736321,
      1.08153340501132,
      -0.0744127885357893,
      0.0164094443541384,
      -0.0680468275301065,
      0.025798857610164,
      -0.000145749861944416
    };
    var Pr = P / 100;
    var sigma = s / 4.4;
    var suma = 0;
    foreach (var _tup_1 in zip(I, J, n))
    {
      var i = _tup_1.Item1;
      var j = _tup_1.Item2;
      var ni = _tup_1.Item3;
      suma += ni * Math.Pow(Pr + 0.187, i) * Math.Pow(sigma - 0.755, j);
    }
    return 0.0028 * suma;
  }

  // Backward equation for region 3b, v=f(P,s)
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //   s : float
  //     Specific entropy, [kJ/kgK]
  //
  //   Returns
  //   -------
  //   v : float
  //     Specific volume, [m³/kg]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Supplementary Release on Backward Equations for the
  //   Functions T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS
  //   Industrial Formulation 1997 for the Thermodynamic Properties of Water and
  //   Steam, http://www.iapws.org/relguide/Supp-Tv%28ph,ps%293-2014.pdf, Eq 9
  //
  //   Examples
  //   --------
  //   >>> _Backward3b_v_Ps(20,5)
  //   0.006262101987
  //   >>> _Backward3b_v_Ps(100,5)
  //   0.002449610757
  //

  public static object _Backward3b_v_Ps(object P, object s)
  {
    var I = new List<int>
    {
      -12,
      -12,
      -12,
      -12,
      -12,
      -12,
      -10,
      -10,
      -10,
      -10,
      -8,
      -5,
      -5,
      -5,
      -4,
      -4,
      -4,
      -4,
      -3,
      -2,
      -2,
      -2,
      -2,
      -2,
      -2,
      0,
      0,
      0,
      1,
      1,
      2
    };
    var J = new List<int>
    {
      0,
      1,
      2,
      3,
      5,
      6,
      0,
      1,
      2,
      4,
      0,
      1,
      2,
      3,
      0,
      1,
      2,
      3,
      1,
      0,
      1,
      2,
      3,
      4,
      12,
      0,
      1,
      2,
      0,
      2,
      2
    };
    var n = new List<double>
    {
      5.91599780322238E-05,
      -0.00185465997137856,
      0.0104190510480013,
      0.0059864730203859,
      -0.771391189901699,
      1.72549765557036,
      -0.000467076079846526,
      0.0134533823384439,
      -0.0808094336805495,
      0.508139374365767,
      0.00128584643361683,
      -1.63899353915435,
      5.86938199318063,
      -2.92466667918613,
      -0.00614076301499537,
      5.76199014049172,
      -12.1613320606788,
      1.67637540957944,
      -7.44135838773463,
      0.0378168091437659,
      4.01432203027688,
      16.0279837479185,
      3.17848779347728,
      -3.58362310304853,
      -1159952.60446827,
      0.199256573577909,
      -0.122270624794624,
      -19.1449143716586,
      -0.0150448002905284,
      14.6407900162154,
      -3.2747778718823
    };
    var Pr = P / 100;
    var sigma = s / 5.3;
    var suma = 0;
    foreach (var _tup_1 in zip(I, J, n))
    {
      var i = _tup_1.Item1;
      var j = _tup_1.Item2;
      var ni = _tup_1.Item3;
      suma += ni * Math.Pow(Pr + 0.298, i) * Math.Pow(sigma - 0.816, j);
    }
    return 0.0088 * suma;
  }

  // Backward equation for region 3, v=f(P,s)
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //   s : float
  //     Specific entropy, [kJ/kgK]
  //
  //   Returns
  //   -------
  //   v : float
  //     Specific volume, [m³/kg]
  //

  public static object _Backward3_v_Ps(object P, object s)
  {
    if(s <= sc) { return _Backward3a_v_Ps(P, s); }
    else { return _Backward3b_v_Ps(P, s); }
  }

  // Backward equation for region 3a, T=f(P,s)
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //   s : float
  //     Specific entropy, [kJ/kgK]
  //
  //   Returns
  //   -------
  //   T : float
  //     Temperature, [K]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Supplementary Release on Backward Equations for the
  //   Functions T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS
  //   Industrial Formulation 1997 for the Thermodynamic Properties of Water and
  //   Steam, http://www.iapws.org/relguide/Supp-Tv%28ph,ps%293-2014.pdf, Eq 6
  //
  //   Examples
  //   --------
  //   >>> _Backward3a_T_Ps(20,3.8)
  //   628.2959869
  //   >>> _Backward3a_T_Ps(100,4)
  //   705.6880237
  //

  public static object _Backward3a_T_Ps(object P, object s)
  {
    var I = new List<int>
    {
      -12,
      -12,
      -10,
      -10,
      -10,
      -10,
      -8,
      -8,
      -8,
      -8,
      -6,
      -6,
      -6,
      -5,
      -5,
      -5,
      -4,
      -4,
      -4,
      -2,
      -2,
      -1,
      -1,
      0,
      0,
      0,
      1,
      2,
      2,
      3,
      8,
      8,
      10
    };
    var J = new List<int>
    {
      28,
      32,
      4,
      10,
      12,
      14,
      5,
      7,
      8,
      28,
      2,
      6,
      32,
      0,
      14,
      32,
      6,
      10,
      36,
      1,
      4,
      1,
      6,
      0,
      1,
      4,
      0,
      0,
      3,
      2,
      0,
      1,
      2
    };
    var n = new List<double>
    {
      1500420082.63875,
      -159397258480.424,
      0.000502181140217975,
      -67.2057767855466,
      1450.58545404456,
      -8238.8953488889,
      -0.154852214233853,
      11.2305046746695,
      -29.7000213482822,
      43856513263.5495,
      0.00137837838635464,
      -2.97478527157462,
      9717779473494.13,
      -5.71527767052398E-05,
      28830.794977842,
      -74442828926270.3,
      12.8017324848921,
      -368.275545889071,
      6.64768904779177E+15,
      0.044935925195888,
      -4.22897836099655,
      -0.240614376434179,
      -4.74341365254924,
      0.72409399912611,
      0.923874349695897,
      3.99043655281015,
      0.0384066651868009,
      -0.00359344365571848,
      -0.735196448821653,
      0.188367048396131,
      0.000141064266818704,
      -0.00257418501496337,
      0.00123220024851555
    };
    var Pr = P / 100;
    var sigma = s / 4.4;
    var suma = 0;
    foreach (var _tup_1 in zip(I, J, n))
    {
      var i = _tup_1.Item1;
      var j = _tup_1.Item2;
      var ni = _tup_1.Item3;
      suma += ni * Math.Pow(Pr + 0.24, i) * Math.Pow(sigma - 0.703, j);
    }
    return 760 * suma;
  }

  // Backward equation for region 3b, T=f(P,s)
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //   s : float
  //     Specific entropy, [kJ/kgK]
  //
  //   Returns
  //   -------
  //   T : float
  //     Temperature, [K]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Supplementary Release on Backward Equations for the
  //   Functions T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS
  //   Industrial Formulation 1997 for the Thermodynamic Properties of Water and
  //   Steam, http://www.iapws.org/relguide/Supp-Tv%28ph,ps%293-2014.pdf, Eq 7
  //
  //   Examples
  //   --------
  //   >>> _Backward3b_T_Ps(20,5)
  //   640.1176443
  //   >>> _Backward3b_T_Ps(100,5)
  //   847.4332825
  //

  public static object _Backward3b_T_Ps(object P, object s)
  {
    var I = new List<int>
    {
      -12,
      -12,
      -12,
      -12,
      -8,
      -8,
      -8,
      -6,
      -6,
      -6,
      -5,
      -5,
      -5,
      -5,
      -5,
      -4,
      -3,
      -3,
      -2,
      0,
      2,
      3,
      4,
      5,
      6,
      8,
      12,
      14
    };
    var J = new List<int>
    {
      1,
      3,
      4,
      7,
      0,
      1,
      3,
      0,
      2,
      4,
      0,
      1,
      2,
      4,
      6,
      12,
      1,
      6,
      2,
      0,
      1,
      1,
      0,
      24,
      0,
      3,
      1,
      2
    };
    var n = new List<double>
    {
      0.52711170160166,
      -40.1317830052742,
      153.020073134484,
      -2247.99398218827,
      -0.193993484669048,
      -1.40467557893768,
      42.6799878114024,
      0.752810643416743,
      22.6657238616417,
      -622.873556909932,
      -0.660823667935396,
      0.841267087271658,
      -25.3717501764397,
      485.708963532948,
      880.531517490555,
      2650155.92794626,
      -0.359287150025783,
      -656.991567673753,
      2.41768149185367,
      0.856873461222588,
      0.655143675313458,
      -0.213535213206406,
      0.00562974957606348,
      -316955725450471.0,
      -0.000699997000152457,
      0.0119845803210767,
      1.93848122022095E-05,
      -2.15095749182309E-05
    };
    var Pr = P / 100;
    var sigma = s / 5.3;
    var suma = 0;
    foreach (var _tup_1 in zip(I, J, n))
    {
      var i = _tup_1.Item1;
      var j = _tup_1.Item2;
      var ni = _tup_1.Item3;
      suma += ni * Math.Pow(Pr + 0.76, i) * Math.Pow(sigma - 0.818, j);
    }
    return 860 * suma;
  }

  // Backward equation for region 3, T=f(P,s)
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //   s : float
  //     Specific entropy, [kJ/kgK]
  //
  //   Returns
  //   -------
  //   T : float
  //     Temperature, [K]
  //

  public static object _Backward3_T_Ps(object P, object s)
  {
    object T;
    var sc = 4.41202148223476;
    if(s <= sc) { T = _Backward3a_T_Ps(P, s); }
    else { T = _Backward3b_T_Ps(P, s); }
    return T;
  }

  // Backward equation for region 3a, P=f(h,s)
  //
  //   Parameters
  //   ----------
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //   s : float
  //     Specific entropy, [kJ/kgK]
  //
  //   Returns
  //   -------
  //   P : float
  //     Pressure, [MPa]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Supplementary Release on Backward Equations p(h,s) for
  //   Region 3, Equations as a Function of h and s for the Region Boundaries, and
  //   an Equation Tsat(h,s) for Region 4 of the IAPWS Industrial Formulation 1997
  //   for the Thermodynamic Properties of Water and Steam,
  //   http://www.iapws.org/relguide/Supp-phs3-2014.pdf. Eq 1
  //
  //   Examples
  //   --------
  //   >>> _Backward3a_P_hs(1700,3.8)
  //   25.55703246
  //   >>> _Backward3a_P_hs(2000,4.2)
  //   45.40873468
  //   >>> _Backward3a_P_hs(2100,4.3)
  //   60.78123340
  //

  public static object _Backward3a_P_hs(object h, object s)
  {
    var I = new List<int>
    {
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
      3,
      3,
      3,
      4,
      4,
      4,
      4,
      5,
      6,
      7,
      8,
      10,
      10,
      14,
      18,
      20,
      22,
      22,
      24,
      28,
      28,
      32,
      32
    };
    var J = new List<int>
    {
      0,
      1,
      5,
      0,
      3,
      4,
      8,
      14,
      6,
      16,
      0,
      2,
      3,
      0,
      1,
      4,
      5,
      28,
      28,
      24,
      1,
      32,
      36,
      22,
      28,
      36,
      16,
      28,
      36,
      16,
      36,
      10,
      28
    };
    var n = new List<double>
    {
      7.70889828326934,
      -26.0835009128688,
      267.416218930389,
      17.2221089496844,
      -293.54233214597,
      614.135601882478,
      -61056.2757725674,
      -65127225.1118219,
      73591.9313521937,
      -11664650591.4191,
      35.5267086434461,
      -596.144543825955,
      -475.842430145708,
      69.6781965359503,
      335.674250377312,
      25052.6809130882,
      146997.380630766,
      5.38069315091534E+19,
      1.43619827291346E+21,
      3.64985866165994E+19,
      -2547.41561156775,
      2.40120197096563E+27,
      -3.93847464679496E+29,
      1.47073407024852E+24,
      -4.26391250432059E+31,
      1.94509340621077E+38,
      6.66212132114896E+23,
      7.06777016552858E+33,
      1.75563621975576E+41,
      1.08408607429124E+28,
      7.30872705175151E+43,
      1.5914584739887E+24,
      3.77121605943324E+40
    };
    var nu = h / 2300;
    var sigma = s / 4.4;
    var suma = 0;
    foreach (var _tup_1 in zip(I, J, n))
    {
      var i = _tup_1.Item1;
      var j = _tup_1.Item2;
      var ni = _tup_1.Item3;
      suma += ni * Math.Pow(nu - 1.01, i) * Math.Pow(sigma - 0.75, j);
    }
    return 99 * suma;
  }

  // Backward equation for region 3b, P=f(h,s)
  //
  //   Parameters
  //   ----------
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //   s : float
  //     Specific entropy, [kJ/kgK]
  //
  //   Returns
  //   -------
  //   P : float
  //     Pressure, [MPa]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Supplementary Release on Backward Equations p(h,s) for
  //   Region 3, Equations as a Function of h and s for the Region Boundaries, and
  //   an Equation Tsat(h,s) for Region 4 of the IAPWS Industrial Formulation 1997
  //   for the Thermodynamic Properties of Water and Steam,
  //   http://www.iapws.org/relguide/Supp-phs3-2014.pdf. Eq 1
  //
  //   Examples
  //   --------
  //   >>> _Backward3b_P_hs(2400,4.7)
  //   63.63924887
  //   >>> _Backward3b_P_hs(2600,5.1)
  //   34.34999263
  //   >>> _Backward3b_P_hs(2700,5.0)
  //   88.39043281
  //

  public static object _Backward3b_P_hs(object h, object s)
  {
    var I = new List<int>
    {
      -12,
      -12,
      -12,
      -12,
      -12,
      -10,
      -10,
      -10,
      -10,
      -8,
      -8,
      -6,
      -6,
      -6,
      -6,
      -5,
      -4,
      -4,
      -4,
      -3,
      -3,
      -3,
      -3,
      -2,
      -2,
      -1,
      0,
      2,
      2,
      5,
      6,
      8,
      10,
      14,
      14
    };
    var J = new List<int>
    {
      2,
      10,
      12,
      14,
      20,
      2,
      10,
      14,
      18,
      2,
      8,
      2,
      6,
      7,
      8,
      10,
      4,
      5,
      8,
      1,
      3,
      5,
      6,
      0,
      1,
      0,
      3,
      0,
      1,
      0,
      1,
      1,
      1,
      3,
      7
    };
    var n = new List<double>
    {
      1.25244360717979E-13,
      -0.0126599322553713,
      5.06878030140626,
      31.7847171154202,
      -391041.161399932,
      -9.75733406392044E-11,
      -18.6312419488279,
      510.973543414101,
      373847.005822362,
      2.99804024666572E-08,
      20.0544393820342,
      -4.98030487662829E-06,
      -10.230180636003,
      55.2819126990325,
      -206.211367510878,
      -7940.12232324823,
      7.82248472028153,
      -58.6544326902468,
      3550.73647696481,
      -0.000115303107290162,
      -1.75092403171802,
      257.98168774816,
      -727.048374179467,
      0.000121644822609198,
      0.0393137871762692,
      0.00704181005909296,
      -82.910820069811,
      -0.26517881813125,
      13.7531682453991,
      -52.2394090753046,
      2405.56298941048,
      -22736.1631268929,
      89074.6343932567,
      -23923456.5822486,
      5687958081.29714
    };
    var nu = h / 2800;
    var sigma = s / 5.3;
    var suma = 0;
    foreach (var _tup_1 in zip(I, J, n))
    {
      var i = _tup_1.Item1;
      var j = _tup_1.Item2;
      var ni = _tup_1.Item3;
      suma += ni * Math.Pow(nu - 0.681, i) * Math.Pow(sigma - 0.792, j);
    }
    return 16.6 / suma;
  }

  // Backward equation for region 3, P=f(h,s)
  //
  //   Parameters
  //   ----------
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //   s : float
  //     Specific entropy, [kJ/kgK]
  //
  //   Returns
  //   -------
  //   P : float
  //     Pressure, [MPa]
  //

  public static object _Backward3_P_hs(object h, object s)
  {
    var sc = 4.41202148223476;
    if(s <= sc) { return _Backward3a_P_hs(h, s); }
    else { return _Backward3b_P_hs(h, s); }
  }

  // Backward equation for region 3 for saturated state, vs=f(P,x)
  //
  //   Parameters
  //   ----------
  //   T : float
  //     Temperature, [K]
  //   P : float
  //     Pressure, [MPa]
  //   x : integer
  //     Vapor quality, [-]
  //
  //   Returns
  //   -------
  //   v : float
  //     Specific volume, [m³/kg]
  //
  //   Notes
  //   -----
  //   The vapor quality (x) can be 0 (saturated liquid) or 1 (saturated vapour)
  //

  public static object _Backward3_sat_v_P(object P, object T, object x)
  {
    object region;
    if(x == 0)
    {
      if(P < 19.00881189) { region = "c"; }
      else if(P < 21.0434) { region = "s"; }
      else if(P < 21.9316) { region = "u"; }
      else { region = "y"; }
    }
    else if(P < 20.5) { region = "t"; }
    else if(P < 21.0434) { region = "r"; }
    else if(P < 21.9009) { region = "x"; }
    else { region = "z"; }
    return _Backward3x_v_PT(T, P, region);
  }

  // Backward equation for region 3, v=f(P,T)
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
  //   v : float
  //     Specific volume, [m³/kg]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Supplementary Release on Backward Equations for Specific
  //   Volume as a Function of Pressure and Temperature v(p,T) for Region 3 of the
  //   IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water
  //   and Steam, http://www.iapws.org/relguide/Supp-VPT3-2016.pdf, Table 2 and 10
  //

  public static object _Backward3_v_PT(object P, object T)
  {
    object region;
    if(P > 40)
    {
      if(T <= _tab_P(P)) { region = "a"; }
      else { region = "b"; }
    }
    else if(25 < P <= 40)
    {
      var tcd = _txx_P(P, "cd");
      var tab = _tab_P(P);
      var tef = _tef_P(P);
      if(T <= tcd) { region = "c"; }
      else if(tcd < T <= tab) { region = "d"; }
      else if(tab < T <= tef) { region = "e"; }
      else { region = "f"; }
    }
    else if(23.5 < P <= 25)
    {
      tcd = _txx_P(P, "cd");
      var tgh = _txx_P(P, "gh");
      tef = _tef_P(P);
      var tij = _txx_P(P, "ij");
      var tjk = _txx_P(P, "jk");
      if(T <= tcd) { region = "c"; }
      else if(tcd < T <= tgh) { region = "g"; }
      else if(tgh < T <= tef) { region = "h"; }
      else if(tef < T <= tij) { region = "i"; }
      else if(tij < T <= tjk) { region = "j"; }
      else { region = "k"; }
    }
    else if(23 < P <= 23.5)
    {
      tcd = _txx_P(P, "cd");
      tgh = _txx_P(P, "gh");
      tef = _tef_P(P);
      tij = _txx_P(P, "ij");
      tjk = _txx_P(P, "jk");
      if(T <= tcd) { region = "c"; }
      else if(tcd < T <= tgh) { region = "l"; }
      else if(tgh < T <= tef) { region = "h"; }
      else if(tef < T <= tij) { region = "i"; }
      else if(tij < T <= tjk) { region = "j"; }
      else { region = "k"; }
    }
    else if(22.5 < P <= 23)
    {
      tcd = _txx_P(P, "cd");
      tgh = _txx_P(P, "gh");
      var tmn = _txx_P(P, "mn");
      tef = _tef_P(P);
      var top = _top_P(P);
      tij = _txx_P(P, "ij");
      tjk = _txx_P(P, "jk");
      if(T <= tcd) { region = "c"; }
      else if(tcd < T <= tgh) { region = "l"; }
      else if(tgh < T <= tmn) { region = "m"; }
      else if(tmn < T <= tef) { region = "n"; }
      else if(tef < T <= top) { region = "o"; }
      else if(top < T <= tij) { region = "p"; }
      else if(tij < T <= tjk) { region = "j"; }
      else { region = "k"; }
    }
    else if(_PSat_T(643.15) < P <= 22.5)
    {
      tcd = _txx_P(P, "cd");
      var tqu = _txx_P(P, "qu");
      var trx = _txx_P(P, "rx");
      tjk = _txx_P(P, "jk");
      if(T <= tcd) { region = "c"; }
      else if(tcd < T <= tqu) { region = "q"; }
      else if(tqu < T <= trx)
      {
        // Table 10
        tef = _tef_P(P);
        var twx = _twx_P(P);
        var tuv = _txx_P(P, "uv");
        if(22.11 < P <= 22.5)
        {
          if(T <= tuv) { region = "u"; }
          else if(tuv <= T <= tef) { region = "v"; }
          else if(tef <= T <= twx) { region = "w"; }
          else { region = "x"; }
        }
        else if(22.064 < P <= 22.11)
        {
          if(T <= tuv) { region = "u"; }
          else if(tuv <= T <= tef) { region = "y"; }
          else if(tef <= T <= twx) { region = "z"; }
          else { region = "x"; }
        }
        else if(T > _TSat_P(P))
        {
          if(_PSat_T(643.15) < P <= 21.90096265) { region = "x"; }
          else if(21.90096265 < P <= 22.064)
          {
            if(T <= twx) { region = "z"; }
            else { region = "x"; }
          }
        }
        else if(T <= _TSat_P(P))
        {
          if(_PSat_T(643.15) < P <= 21.93161551) { region = "u"; }
          else if(21.93161551 < P <= 22.064)
          {
            if(T <= tuv) { region = "u"; }
            else { region = "y"; }
          }
        }
      }
      else if(trx < T <= tjk) { region = "r"; }
      else { region = "k"; }
    }
    else if(20.5 < P <= _PSat_T(643.15))
    {
      tcd = _txx_P(P, "cd");
      var Ts = _TSat_P(P);
      tjk = _txx_P(P, "jk");
      if(T <= tcd) { region = "c"; }
      else if(tcd < T <= Ts) { region = "s"; }
      else if(Ts < T <= tjk) { region = "r"; }
      else { region = "k"; }
    }
    else if(19.0088118917393 < P <= 20.5)
    {
      tcd = _txx_P(P, "cd");
      Ts = _TSat_P(P);
      if(T <= tcd) { region = "c"; }
      else if(tcd < T <= Ts) { region = "s"; }
      else { region = "t"; }
    }
    else if(Ps_623 < P <= 19.0088118917393)
    {
      Ts = _TSat_P(P);
      if(T <= Ts) { region = "c"; }
      else { region = "t"; }
    }
    return _Backward3x_v_PT(T, P, region);
  }

  // Backward equation for region 3x, v=f(P,T)
  //
  //   Parameters
  //   ----------
  //   T : float
  //     Temperature, [K]
  //   P : float
  //     Pressure, [MPa]
  //   x : char
  //     Region 3 subregion code
  //
  //   Returns
  //   -------
  //   v : float
  //     Specific volume, [m³/kg]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Supplementary Release on Backward Equations for Specific
  //   Volume as a Function of Pressure and Temperature v(p,T) for Region 3 of the
  //   IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water
  //   and Steam, http://www.iapws.org/relguide/Supp-VPT3-2016.pdf, Eq. 4-5
  //
  //   Examples
  //   --------
  //   >>> _Backward3x_v_PT(630,50,"a")
  //   0.001470853100
  //   >>> _Backward3x_v_PT(670,80,"a")
  //   0.001503831359
  //   >>> _Backward3x_v_PT(710,50,"b")
  //   0.002204728587
  //   >>> _Backward3x_v_PT(750,80,"b")
  //   0.001973692940
  //   >>> _Backward3x_v_PT(630,20,"c")
  //   0.001761696406
  //   >>> _Backward3x_v_PT(650,30,"c")
  //   0.001819560617
  //   >>> _Backward3x_v_PT(656,26,"d")
  //   0.002245587720
  //   >>> _Backward3x_v_PT(670,30,"d")
  //   0.002506897702
  //   >>> _Backward3x_v_PT(661,26,"e")
  //   0.002970225962
  //   >>> _Backward3x_v_PT(675,30,"e")
  //   0.003004627086
  //   >>> _Backward3x_v_PT(671,26,"f")
  //   0.005019029401
  //   >>> _Backward3x_v_PT(690,30,"f")
  //   0.004656470142
  //   >>> _Backward3x_v_PT(649,23.6,"g")
  //   0.002163198378
  //   >>> _Backward3x_v_PT(650,24,"g")
  //   0.002166044161
  //   >>> _Backward3x_v_PT(652,23.6,"h")
  //   0.002651081407
  //   >>> _Backward3x_v_PT(654,24,"h")
  //   0.002967802335
  //   >>> _Backward3x_v_PT(653,23.6,"i")
  //   0.003273916816
  //   >>> _Backward3x_v_PT(655,24,"i")
  //   0.003550329864
  //   >>> _Backward3x_v_PT(655,23.5,"j")
  //   0.004545001142
  //   >>> _Backward3x_v_PT(660,24,"j")
  //   0.005100267704
  //   >>> _Backward3x_v_PT(660,23,"k")
  //   0.006109525997
  //   >>> _Backward3x_v_PT(670,24,"k")
  //   0.006427325645
  //   >>> _Backward3x_v_PT(646,22.6,"l")
  //   0.002117860851
  //   >>> _Backward3x_v_PT(646,23,"l")
  //   0.002062374674
  //   >>> _Backward3x_v_PT(648.6,22.6,"m")
  //   0.002533063780
  //   >>> _Backward3x_v_PT(649.3,22.8,"m")
  //   0.002572971781
  //   >>> _Backward3x_v_PT(649,22.6,"n")
  //   0.002923432711
  //   >>> _Backward3x_v_PT(649.7,22.8,"n")
  //   0.002913311494
  //   >>> _Backward3x_v_PT(649.1,22.6,"o")
  //   0.003131208996
  //   >>> _Backward3x_v_PT(649.9,22.8,"o")
  //   0.003221160278
  //   >>> _Backward3x_v_PT(649.4,22.6,"p")
  //   0.003715596186
  //   >>> _Backward3x_v_PT(650.2,22.8,"p")
  //   0.003664754790
  //   >>> _Backward3x_v_PT(640,21.1,"q")
  //   0.001970999272
  //   >>> _Backward3x_v_PT(643,21.8,"q")
  //   0.002043919161
  //   >>> _Backward3x_v_PT(644,21.1,"r")
  //   0.005251009921
  //   >>> _Backward3x_v_PT(648,21.8,"r")
  //   0.005256844741
  //   >>> _Backward3x_v_PT(635,19.1,"s")
  //   0.001932829079
  //   >>> _Backward3x_v_PT(638,20,"s")
  //   0.001985387227
  //   >>> _Backward3x_v_PT(626,17,"t")
  //   0.008483262001
  //   >>> _Backward3x_v_PT(640,20,"t")
  //   0.006227528101
  //   >>> _Backward3x_v_PT(644.6,21.5,"u")
  //   0.002268366647
  //   >>> _Backward3x_v_PT(646.1,22,"u")
  //   0.002296350553
  //   >>> _Backward3x_v_PT(648.6,22.5,"v")
  //   0.002832373260
  //   >>> _Backward3x_v_PT(647.9,22.3,"v")
  //   0.002811424405
  //   >>> _Backward3x_v_PT(647.5,22.15,"w")
  //   0.003694032281
  //   >>> _Backward3x_v_PT(648.1,22.3,"w")
  //   0.003622226305
  //   >>> _Backward3x_v_PT(648,22.11,"x")
  //   0.004528072649
  //   >>> _Backward3x_v_PT(649,22.3,"x")
  //   0.004556905799
  //   >>> _Backward3x_v_PT(646.84,22,"y")
  //   0.002698354719
  //   >>> _Backward3x_v_PT(647.05,22.064,"y")
  //   0.002717655648
  //   >>> _Backward3x_v_PT(646.89,22,"z")
  //   0.003798732962
  //   >>> _Backward3x_v_PT(647.15,22.064,"z")
  //   0.003701940009
  //

  public static object _Backward3x_v_PT(object T, object P, object x)
  {
    object ni;
    object j;
    object i;
    var par = new Dictionary<object, object>
    {
      { "a", new List<double> { 0.0024, 100, 760, 0.085, 0.817, 1, 1, 1 }},
      { "b", new List<double> { 0.0041, 100, 860, 0.28, 0.779, 1, 1, 1 }},
      { "c", new List<double> { 0.0022, 40, 690, 0.259, 0.903, 1, 1, 1 }},
      { "d", new List<double> { 0.0029, 40, 690, 0.559, 0.939, 1, 1, 4 }},
      { "e", new List<double> { 0.0032, 40, 710, 0.587, 0.918, 1, 1, 1 }},
      { "f", new List<double> { 0.0064, 40, 730, 0.587, 0.891, 0.5, 1, 4 }},
      { "g", new List<double> { 0.0027, 25, 660, 0.872, 0.971, 1, 1, 4 }},
      { "h", new List<double> { 0.0032, 25, 660, 0.898, 0.983, 1, 1, 4 }},
      { "i", new List<double> { 0.0041, 25, 660, 0.91, 0.984, 0.5, 1, 4 }},
      { "j", new List<double> { 0.0054, 25, 670, 0.875, 0.964, 0.5, 1, 4 }},
      { "k", new List<double> { 0.0077, 25, 680, 0.802, 0.935, 1, 1, 1 }},
      { "l", new List<double> { 0.0026, 24, 650, 0.908, 0.989, 1, 1, 4 }},
      { "m", new List<double> { 0.0028, 23, 650, 1.0, 0.997, 1, 0.25, 1 }},
      { "n", new List<double> { 0.0031, 23, 650, 0.976, 0.997, null, null, null }},
      { "o", new List<double> { 0.0034, 23, 650, 0.974, 0.996, 0.5, 1, 1 }},
      { "p", new List<double> { 0.0041, 23, 650, 0.972, 0.997, 0.5, 1, 1 }},
      { "q", new List<double> { 0.0022, 23, 650, 0.848, 0.983, 1, 1, 4 }},
      { "r", new List<double> { 0.0054, 23, 650, 0.874, 0.982, 1, 1, 1 }},
      { "s", new List<double> { 0.0022, 21, 640, 0.886, 0.99, 1, 1, 4 }},
      { "t", new List<double> { 0.0088, 20, 650, 0.803, 1.02, 1, 1, 1 }},
      { "u", new List<double> { 0.0026, 23, 650, 0.902, 0.988, 1, 1, 1 }},
      { "v", new List<double> { 0.0031, 23, 650, 0.96, 0.995, 1, 1, 1 }},
      { "w", new List<double> { 0.0039, 23, 650, 0.959, 0.995, 1, 1, 4 }},
      { "x", new List<double> { 0.0049, 23, 650, 0.91, 0.988, 1, 1, 1 }},
      { "y", new List<double> { 0.0031, 22, 650, 0.996, 0.994, 1, 1, 4 }},
      { "z", new List<double> { 0.0038, 22, 650, 0.993, 0.994, 1, 1, 4 }}
    };
    var I = new Dictionary<object, object>
    {
      {
        "a",
        new List<int>
        {
          -12,
          -12,
          -12,
          -10,
          -10,
          -10,
          -8,
          -8,
          -8,
          -6,
          -5,
          -5,
          -5,
          -4,
          -3,
          -3,
          -3,
          -3,
          -2,
          -2,
          -2,
          -1,
          -1,
          -1,
          0,
          0,
          1,
          1,
          2,
          2
        }},
      {
        "b",
        new List<int>
        {
          -12,
          -12,
          -10,
          -10,
          -8,
          -6,
          -6,
          -6,
          -5,
          -5,
          -5,
          -4,
          -4,
          -4,
          -3,
          -3,
          -3,
          -3,
          -3,
          -2,
          -2,
          -2,
          -1,
          -1,
          0,
          0,
          1,
          1,
          2,
          3,
          4,
          4
        }},
      {
        "c",
        new List<int>
        {
          -12,
          -12,
          -12,
          -10,
          -10,
          -10,
          -8,
          -8,
          -8,
          -6,
          -5,
          -5,
          -5,
          -4,
          -4,
          -3,
          -3,
          -2,
          -2,
          -2,
          -1,
          -1,
          -1,
          0,
          0,
          0,
          1,
          1,
          2,
          2,
          2,
          2,
          3,
          3,
          8
        }},
      {
        "d",
        new List<int>
        {
          -12,
          -12,
          -12,
          -12,
          -12,
          -12,
          -10,
          -10,
          -10,
          -10,
          -10,
          -10,
          -10,
          -8,
          -8,
          -8,
          -8,
          -6,
          -6,
          -5,
          -5,
          -5,
          -5,
          -4,
          -4,
          -4,
          -3,
          -3,
          -2,
          -2,
          -1,
          -1,
          -1,
          0,
          0,
          1,
          1,
          3
        }},
      {
        "e",
        new List<int>
        {
          -12,
          -12,
          -10,
          -10,
          -10,
          -10,
          -10,
          -8,
          -8,
          -8,
          -6,
          -5,
          -4,
          -4,
          -3,
          -3,
          -3,
          -2,
          -2,
          -2,
          -2,
          -1,
          0,
          0,
          1,
          1,
          1,
          2,
          2
        }},
      {
        "f",
        new List<int>
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
          2,
          2,
          3,
          3,
          3,
          4,
          5,
          5,
          6,
          7,
          7,
          10,
          12,
          12,
          12,
          14,
          14,
          14,
          14,
          14,
          16,
          16,
          18,
          18,
          20,
          20,
          20,
          22,
          24,
          24,
          28,
          32
        }},
      {
        "g",
        new List<int>
        {
          -12,
          -12,
          -12,
          -12,
          -12,
          -12,
          -10,
          -10,
          -10,
          -8,
          -8,
          -8,
          -8,
          -6,
          -6,
          -5,
          -5,
          -4,
          -3,
          -2,
          -2,
          -2,
          -2,
          -1,
          -1,
          -1,
          0,
          0,
          0,
          1,
          1,
          1,
          3,
          5,
          6,
          8,
          10,
          10
        }},
      {
        "h",
        new List<int>
        {
          -12,
          -12,
          -10,
          -10,
          -10,
          -10,
          -10,
          -10,
          -8,
          -8,
          -8,
          -8,
          -8,
          -6,
          -6,
          -6,
          -5,
          -5,
          -5,
          -4,
          -4,
          -3,
          -3,
          -2,
          -1,
          -1,
          0,
          1,
          1
        }},
      {
        "i",
        new List<int>
        {
          0,
          0,
          0,
          1,
          1,
          1,
          1,
          2,
          3,
          3,
          4,
          4,
          4,
          5,
          5,
          5,
          7,
          7,
          8,
          8,
          10,
          12,
          12,
          12,
          14,
          14,
          14,
          14,
          18,
          18,
          18,
          18,
          18,
          20,
          20,
          22,
          24,
          24,
          32,
          32,
          36,
          36
        }},
      {
        "j",
        new List<int>
        {
          0,
          0,
          0,
          1,
          1,
          1,
          2,
          2,
          3,
          4,
          4,
          5,
          5,
          5,
          6,
          10,
          12,
          12,
          14,
          14,
          14,
          16,
          18,
          20,
          20,
          24,
          24,
          28,
          28
        }},
      {
        "k",
        new List<int>
        {
          -2,
          -2,
          -1,
          -1,
          0,
          -0,
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
          2,
          2,
          2,
          2,
          2,
          2,
          5,
          5,
          5,
          6,
          6,
          6,
          6,
          8,
          10,
          12
        }},
      {
        "l",
        new List<int>
        {
          -12,
          -12,
          -12,
          -12,
          -12,
          -10,
          -10,
          -8,
          -8,
          -8,
          -8,
          -8,
          -8,
          -8,
          -6,
          -5,
          -5,
          -4,
          -4,
          -3,
          -3,
          -3,
          -3,
          -2,
          -2,
          -2,
          -1,
          -1,
          -1,
          0,
          0,
          0,
          0,
          1,
          1,
          2,
          4,
          5,
          5,
          6,
          10,
          10,
          14
        }},
      {
        "m",
        new List<int>
        {
          0,
          3,
          8,
          20,
          1,
          3,
          4,
          5,
          1,
          6,
          2,
          4,
          14,
          2,
          5,
          3,
          0,
          1,
          1,
          1,
          28,
          2,
          16,
          0,
          5,
          0,
          3,
          4,
          12,
          16,
          1,
          8,
          14,
          0,
          2,
          3,
          4,
          8,
          14,
          24
        }},
      {
        "n",
        new List<int>
        {
          0,
          3,
          4,
          6,
          7,
          10,
          12,
          14,
          18,
          0,
          3,
          5,
          6,
          8,
          12,
          0,
          3,
          7,
          12,
          2,
          3,
          4,
          2,
          4,
          7,
          4,
          3,
          5,
          6,
          0,
          0,
          3,
          1,
          0,
          1,
          0,
          1,
          0,
          1
        }},
      {
        "o",
        new List<int>
        {
          0,
          0,
          0,
          2,
          3,
          4,
          4,
          4,
          4,
          4,
          5,
          5,
          6,
          7,
          8,
          8,
          8,
          10,
          10,
          14,
          14,
          20,
          20,
          24
        }},
      {
        "p",
        new List<int>
        {
          0,
          0,
          0,
          0,
          1,
          2,
          3,
          3,
          4,
          6,
          7,
          7,
          8,
          10,
          12,
          12,
          12,
          14,
          14,
          14,
          16,
          18,
          20,
          22,
          24,
          24,
          36
        }},
      {
        "q",
        new List<int>
        {
          -12,
          -12,
          -10,
          -10,
          -10,
          -10,
          -8,
          -6,
          -5,
          -5,
          -4,
          -4,
          -3,
          -2,
          -2,
          -2,
          -2,
          -1,
          -1,
          -1,
          0,
          1,
          1,
          1
        }},
      {
        "r",
        new List<int>
        {
          -8,
          -8,
          -3,
          -3,
          -3,
          -3,
          -3,
          0,
          0,
          0,
          0,
          3,
          3,
          8,
          8,
          8,
          8,
          10,
          10,
          10,
          10,
          10,
          10,
          10,
          10,
          12,
          14
        }},
      {
        "s",
        new List<int>
        {
          -12,
          -12,
          -10,
          -8,
          -6,
          -5,
          -5,
          -4,
          -4,
          -3,
          -3,
          -2,
          -1,
          -1,
          -1,
          0,
          0,
          0,
          0,
          1,
          1,
          3,
          3,
          3,
          4,
          4,
          4,
          5,
          14
        }},
      {
        "t",
        new List<int>
        {
          0,
          0,
          0,
          0,
          1,
          1,
          2,
          2,
          2,
          3,
          3,
          4,
          4,
          7,
          7,
          7,
          7,
          7,
          10,
          10,
          10,
          10,
          10,
          18,
          20,
          22,
          22,
          24,
          28,
          32,
          32,
          32,
          36
        }},
      {
        "u",
        new List<int>
        {
          -12,
          -10,
          -10,
          -10,
          -8,
          -8,
          -8,
          -6,
          -6,
          -5,
          -5,
          -5,
          -3,
          -1,
          -1,
          -1,
          -1,
          0,
          0,
          1,
          2,
          2,
          3,
          5,
          5,
          5,
          6,
          6,
          8,
          8,
          10,
          12,
          12,
          12,
          14,
          14,
          14,
          14
        }},
      {
        "v",
        new List<int>
        {
          -10,
          -8,
          -6,
          -6,
          -6,
          -6,
          -6,
          -6,
          -5,
          -5,
          -5,
          -5,
          -5,
          -5,
          -4,
          -4,
          -4,
          -4,
          -3,
          -3,
          -3,
          -2,
          -2,
          -1,
          -1,
          0,
          0,
          0,
          1,
          1,
          3,
          4,
          4,
          4,
          5,
          8,
          10,
          12,
          14
        }},
      {
        "w",
        new List<int>
        {
          -12,
          -12,
          -10,
          -10,
          -8,
          -8,
          -8,
          -6,
          -6,
          -6,
          -6,
          -5,
          -4,
          -4,
          -3,
          -3,
          -2,
          -2,
          -1,
          -1,
          -1,
          0,
          0,
          1,
          2,
          2,
          3,
          3,
          5,
          5,
          5,
          8,
          8,
          10,
          10
        }},
      {
        "x",
        new List<int>
        {
          -8,
          -6,
          -5,
          -4,
          -4,
          -4,
          -3,
          -3,
          -1,
          0,
          0,
          0,
          1,
          1,
          2,
          3,
          3,
          3,
          4,
          5,
          5,
          5,
          6,
          8,
          8,
          8,
          8,
          10,
          12,
          12,
          12,
          12,
          14,
          14,
          14,
          14
        }},
      {
        "y",
        new List<int>
        {
          0,
          0,
          0,
          0,
          1,
          2,
          2,
          2,
          2,
          3,
          3,
          3,
          4,
          4,
          5,
          5,
          8,
          8,
          10,
          12
        }},
      {
        "z",
        new List<int>
        {
          -8,
          -6,
          -5,
          -5,
          -4,
          -4,
          -4,
          -3,
          -3,
          -3,
          -2,
          -1,
          0,
          1,
          2,
          3,
          3,
          6,
          6,
          6,
          6,
          8,
          8
        }}};
    var J = new Dictionary<object, object>
    {
      {
        "a",
        new List<int>
        {
          5,
          10,
          12,
          5,
          10,
          12,
          5,
          8,
          10,
          1,
          1,
          5,
          10,
          8,
          0,
          1,
          3,
          6,
          0,
          2,
          3,
          0,
          1,
          2,
          0,
          1,
          0,
          2,
          0,
          2
        }},
      {
        "b",
        new List<int>
        {
          10,
          12,
          8,
          14,
          8,
          5,
          6,
          8,
          5,
          8,
          10,
          2,
          4,
          5,
          0,
          1,
          2,
          3,
          5,
          0,
          2,
          5,
          0,
          2,
          0,
          1,
          0,
          2,
          0,
          2,
          0,
          1
        }},
      {
        "c",
        new List<int>
        {
          6,
          8,
          10,
          6,
          8,
          10,
          5,
          6,
          7,
          8,
          1,
          4,
          7,
          2,
          8,
          0,
          3,
          0,
          4,
          5,
          0,
          1,
          2,
          0,
          1,
          2,
          0,
          2,
          0,
          1,
          3,
          7,
          0,
          7,
          1
        }},
      {
        "d",
        new List<int>
        {
          4,
          6,
          7,
          10,
          12,
          16,
          0,
          2,
          4,
          6,
          8,
          10,
          14,
          3,
          7,
          8,
          10,
          6,
          8,
          1,
          2,
          5,
          7,
          0,
          1,
          7,
          2,
          4,
          0,
          1,
          0,
          1,
          5,
          0,
          2,
          0,
          6,
          0
        }},
      {
        "e",
        new List<int>
        {
          14,
          16,
          3,
          6,
          10,
          14,
          16,
          7,
          8,
          10,
          6,
          6,
          2,
          4,
          2,
          6,
          7,
          0,
          1,
          3,
          4,
          0,
          0,
          1,
          0,
          4,
          6,
          0,
          2
        }},
      {
        "f",
        new List<int>
        {
          -3,
          -2,
          -1,
          0,
          1,
          2,
          -1,
          1,
          2,
          3,
          0,
          1,
          -5,
          -2,
          0,
          -3,
          -8,
          1,
          -6,
          -4,
          1,
          -6,
          -10,
          -8,
          -4,
          -12,
          -10,
          -8,
          -6,
          -4,
          -10,
          -8,
          -12,
          -10,
          -12,
          -10,
          -6,
          -12,
          -12,
          -4,
          -12,
          -12
        }},
      {
        "g",
        new List<int>
        {
          7,
          12,
          14,
          18,
          22,
          24,
          14,
          20,
          24,
          7,
          8,
          10,
          12,
          8,
          22,
          7,
          20,
          22,
          7,
          3,
          5,
          14,
          24,
          2,
          8,
          18,
          0,
          1,
          2,
          0,
          1,
          3,
          24,
          22,
          12,
          3,
          0,
          6
        }},
      {
        "h",
        new List<int>
        {
          8,
          12,
          4,
          6,
          8,
          10,
          14,
          16,
          0,
          1,
          6,
          7,
          8,
          4,
          6,
          8,
          2,
          3,
          4,
          2,
          4,
          1,
          2,
          0,
          0,
          2,
          0,
          0,
          2
        }},
      {
        "i",
        new List<int>
        {
          0,
          1,
          10,
          -4,
          -2,
          -1,
          0,
          0,
          -5,
          0,
          -3,
          -2,
          -1,
          -6,
          -1,
          12,
          -4,
          -3,
          -6,
          10,
          -8,
          -12,
          -6,
          -4,
          -10,
          -8,
          -4,
          5,
          -12,
          -10,
          -8,
          -6,
          2,
          -12,
          -10,
          -12,
          -12,
          -8,
          -10,
          -5,
          -10,
          -8
        }},
      {
        "j",
        new List<int>
        {
          -1,
          0,
          1,
          -2,
          -1,
          1,
          -1,
          1,
          -2,
          -2,
          2,
          -3,
          -2,
          0,
          3,
          -6,
          -8,
          -3,
          -10,
          -8,
          -5,
          -10,
          -12,
          -12,
          -10,
          -12,
          -6,
          -12,
          -5
        }},
      {
        "k",
        new List<int>
        {
          10,
          12,
          -5,
          6,
          -12,
          -6,
          -2,
          -1,
          0,
          1,
          2,
          3,
          14,
          -3,
          -2,
          0,
          1,
          2,
          -8,
          -6,
          -3,
          -2,
          0,
          4,
          -12,
          -6,
          -3,
          -12,
          -10,
          -8,
          -5,
          -12,
          -12,
          -10
        }},
      {
        "l",
        new List<int>
        {
          14,
          16,
          18,
          20,
          22,
          14,
          24,
          6,
          10,
          12,
          14,
          18,
          24,
          36,
          8,
          4,
          5,
          7,
          16,
          1,
          3,
          18,
          20,
          2,
          3,
          10,
          0,
          1,
          3,
          0,
          1,
          2,
          12,
          0,
          16,
          1,
          0,
          0,
          1,
          14,
          4,
          12,
          10
        }},
      {
        "m",
        new List<int>
        {
          0,
          0,
          0,
          2,
          5,
          5,
          5,
          5,
          6,
          6,
          7,
          8,
          8,
          10,
          10,
          12,
          14,
          14,
          18,
          20,
          20,
          22,
          22,
          24,
          24,
          28,
          28,
          28,
          28,
          28,
          32,
          32,
          32,
          36,
          36,
          36,
          36,
          36,
          36,
          36
        }},
      {
        "n",
        new List<int>
        {
          -12,
          -12,
          -12,
          -12,
          -12,
          -12,
          -12,
          -12,
          -12,
          -10,
          -10,
          -10,
          -10,
          -10,
          -10,
          -8,
          -8,
          -8,
          -8,
          -6,
          -6,
          -6,
          -5,
          -5,
          -5,
          -4,
          -3,
          -3,
          -3,
          -2,
          -1,
          -1,
          0,
          1,
          1,
          2,
          4,
          5,
          6
        }},
      {
        "o",
        new List<int>
        {
          -12,
          -4,
          -1,
          -1,
          -10,
          -12,
          -8,
          -5,
          -4,
          -1,
          -4,
          -3,
          -8,
          -12,
          -10,
          -8,
          -4,
          -12,
          -8,
          -12,
          -8,
          -12,
          -10,
          -12
        }},
      {
        "p",
        new List<int>
        {
          -1,
          0,
          1,
          2,
          1,
          -1,
          -3,
          0,
          -2,
          -2,
          -5,
          -4,
          -2,
          -3,
          -12,
          -6,
          -5,
          -10,
          -8,
          -3,
          -8,
          -8,
          -10,
          -10,
          -12,
          -8,
          -12
        }},
      {
        "q",
        new List<int>
        {
          10,
          12,
          6,
          7,
          8,
          10,
          8,
          6,
          2,
          5,
          3,
          4,
          3,
          0,
          1,
          2,
          4,
          0,
          1,
          2,
          0,
          0,
          1,
          3
        }},
      {
        "r",
        new List<int>
        {
          6,
          14,
          -3,
          3,
          4,
          5,
          8,
          -1,
          0,
          1,
          5,
          -6,
          -2,
          -12,
          -10,
          -8,
          -5,
          -12,
          -10,
          -8,
          -6,
          -5,
          -4,
          -3,
          -2,
          -12,
          -12
        }},
      {
        "s",
        new List<int>
        {
          20,
          24,
          22,
          14,
          36,
          8,
          16,
          6,
          32,
          3,
          8,
          4,
          1,
          2,
          3,
          0,
          1,
          4,
          28,
          0,
          32,
          0,
          1,
          2,
          3,
          18,
          24,
          4,
          24
        }},
      {
        "t",
        new List<int>
        {
          0,
          1,
          4,
          12,
          0,
          10,
          0,
          6,
          14,
          3,
          8,
          0,
          10,
          3,
          4,
          7,
          20,
          36,
          10,
          12,
          14,
          16,
          22,
          18,
          32,
          22,
          36,
          24,
          28,
          22,
          32,
          36,
          36
        }},
      {
        "u",
        new List<int>
        {
          14,
          10,
          12,
          14,
          10,
          12,
          14,
          8,
          12,
          4,
          8,
          12,
          2,
          -1,
          1,
          12,
          14,
          -3,
          1,
          -2,
          5,
          10,
          -5,
          -4,
          2,
          3,
          -5,
          2,
          -8,
          8,
          -4,
          -12,
          -4,
          4,
          -12,
          -10,
          -6,
          6
        }},
      {
        "v",
        new List<int>
        {
          -8,
          -12,
          -12,
          -3,
          5,
          6,
          8,
          10,
          1,
          2,
          6,
          8,
          10,
          14,
          -12,
          -10,
          -6,
          10,
          -3,
          10,
          12,
          2,
          4,
          -2,
          0,
          -2,
          6,
          10,
          -12,
          -10,
          3,
          -6,
          3,
          10,
          2,
          -12,
          -2,
          -3,
          1
        }},
      {
        "w",
        new List<int>
        {
          8,
          14,
          -1,
          8,
          6,
          8,
          14,
          -4,
          -3,
          2,
          8,
          -10,
          -1,
          3,
          -10,
          3,
          1,
          2,
          -8,
          -4,
          1,
          -12,
          1,
          -1,
          -1,
          2,
          -12,
          -5,
          -10,
          -8,
          -6,
          -12,
          -10,
          -12,
          -8
        }},
      {
        "x",
        new List<int>
        {
          14,
          10,
          10,
          1,
          2,
          14,
          -2,
          12,
          5,
          0,
          4,
          10,
          -10,
          -1,
          6,
          -12,
          0,
          8,
          3,
          -6,
          -2,
          1,
          1,
          -6,
          -3,
          1,
          8,
          -8,
          -10,
          -8,
          -5,
          -4,
          -12,
          -10,
          -8,
          -6
        }},
      {
        "y",
        new List<int>
        {
          -3,
          1,
          5,
          8,
          8,
          -4,
          -1,
          4,
          5,
          -8,
          4,
          8,
          -6,
          6,
          -2,
          1,
          -8,
          -2,
          -5,
          -8
        }},
      {
        "z",
        new List<int>
        {
          3,
          6,
          6,
          8,
          5,
          6,
          8,
          -2,
          5,
          6,
          2,
          -6,
          3,
          1,
          6,
          -6,
          -2,
          -6,
          -5,
          -4,
          -1,
          -8,
          -4
        }}};
    var n = new Dictionary<object, object>
    {
      {
        "a",
        new List<double>
        {
          0.00110879558823853,
          572.616740810616,
          -76705.1948380852,
          -0.0253321069529674,
          6280.08049345689,
          234105.654131876,
          0.216867826045856,
          -156.237904341963,
          -26989.3956176613,
          -0.000180407100085505,
          0.00116732227668261,
          26.698704085604,
          28277.6617243286,
          -2424.31520029523,
          0.000435217323022733,
          -0.0122494831387441,
          1.79357604019989,
          44.2729521058314,
          -0.00593223489018342,
          0.453186261685774,
          1.3582570312914,
          0.0408748415856745,
          0.474686397863312,
          1.18646814997915,
          0.546987265727549,
          0.195266770452643,
          -0.0502268790869663,
          -0.369645308193377,
          0.0063382803752842,
          0.0797441793901017
        }},
      {
        "b",
        new List<double>
        {
          -0.0827670470003621,
          41.6887126010565,
          0.0483651982197059,
          -29103.2084950276,
          -111.422582236948,
          -0.0202300083904014,
          294.002509338515,
          140.244997609658,
          -344.384158811459,
          361.182452612149,
          -1406.99677420738,
          -0.00202023902676481,
          171.346792457471,
          -4.25597804058632,
          6.91346085000334E-06,
          0.00151140509678925,
          -0.0416375290166236,
          -41.3754957011042,
          -50.6673295721637,
          -0.000572212965569023,
          6.08817368401785,
          23.9600660256161,
          0.0122261479925384,
          2.16356057692938,
          0.398198903368642,
          -0.116892827834085,
          -0.102845919373532,
          -0.492676637589284,
          0.065554045640679,
          -0.24046253507853,
          -0.0269798180310075,
          0.128369435967012
        }},
      {
        "c",
        new List<double>
        {
          3.1196778876303,
          27671.3458847564,
          32258310.3403269,
          -342.416065095363,
          -899732.529907377,
          -79389204.9821251,
          95.3193003217388,
          2297.84742345072,
          175336.675322499,
          7912143.65222792,
          3.19933345844209E-05,
          -65.9508863555767,
          -833426.563212851,
          0.0645734680583292,
          -3820310.20570813,
          4.06398848470079E-05,
          31.0327498492008,
          -0.000892996718483724,
          234.604891591616,
          3775.15668966951,
          0.0158646812591361,
          0.707906336241843,
          12.601622514657,
          0.736143655772152,
          0.676544268999101,
          -17.8100588189137,
          -0.156531975531713,
          11.7707430048158,
          0.0840143653860447,
          -0.186442467471949,
          -44.0170203949645,
          1232904.23502494,
          -0.0240650039730845,
          -1070777.16660869,
          0.0438319858566475
        }},
      {
        "d",
        new List<double>
        {
          -4.52484847171645E-10,
          3.15210389538801E-05,
          -0.00214991352047545,
          508.058874808345,
          -12712303.6845932,
          1153711331204.97,
          -1.97805728776273E-16,
          2.41554806033972E-11,
          -1.56481703640525E-06,
          0.00277211346836625,
          -20.3578994462286,
          1443694.89909053,
          -41125421794.6539,
          6.23449786243773E-06,
          -22.1774281146038,
          -68931.5087933158,
          -19541952.5060713,
          3163.73510564015,
          2240407.54426988,
          -4.36701347922356E-06,
          -0.000404213852833996,
          -348.153203414663,
          -385294.213555289,
          1.35203700099403E-07,
          0.000134648383271089,
          125031.835351736,
          0.0968123678455841,
          225.660517512438,
          -0.000190102435341872,
          -0.0299628410819229,
          0.00500833915372121,
          0.387842482998411,
          -1385.35367777182,
          0.870745245971773,
          1.71946252068742,
          -0.0326650121426383,
          4980.44171727877,
          0.00551478022765087
        }},
      {
        "e",
        new List<double>
        {
          715815808.404721,
          -114328360753.449,
          3.7653100201572E-12,
          -9.03983668691157E-05,
          665695.908836252,
          5353641749.60127,
          79497740233.5603,
          92.2230563421437,
          -142586.073991215,
          -1117963.81424162,
          8961.2162964076,
          -6699.89239070491,
          0.00451242538486834,
          -33.9731325977713,
          -1.20523111552278,
          47599.2667717124,
          -266627.750390341,
          -0.000153314954386524,
          0.305638404828265,
          123.654999499486,
          -1043.90794213011,
          -0.0157496516174308,
          0.685331118940253,
          1.78373462873903,
          -0.54467412487891,
          2045.29931318843,
          -22834.2359328752,
          0.413197481515899,
          -34.1931835910405
        }},
      {
        "f",
        new List<double>
        {
          -2.51756547792325E-08,
          6.01307193668763E-06,
          -0.00100615977450049,
          0.999969140252192,
          2.14107759236486,
          -16.5175571959086,
          -0.00141987303638727,
          2.69251915156554,
          34.9741815858722,
          -30.0208695771783,
          -1.31546288252539,
          -8.39091277286169,
          1.81545608337015E-10,
          -0.000591099206478909,
          1.52115067087106,
          2.52956470663225E-05,
          1.00726265203786E-15,
          -1.4977453386065,
          -7.93940970562969E-10,
          -0.000150290891264717,
          1.51205531275133,
          4.70942606221652E-06,
          1.95049710391712E-13,
          -9.11627886266077E-09,
          0.000604374640201265,
          -2.25132933900136E-16,
          6.10916973582981E-12,
          -3.03063908043404E-07,
          -1.37796070798409E-05,
          -0.000919296736666106,
          6.39288223132545E-10,
          7.53259479898699E-07,
          -4.00321478682929E-13,
          7.56140294351614E-09,
          -9.12082054034891E-12,
          -2.37612381140539E-08,
          2.69586010591874E-05,
          -7.32828135157839E-11,
          2.4199557830666E-10,
          -0.000405735532730322,
          1.89424143498011E-10,
          -4.86632965074563E-10
        }},
      {
        "g",
        new List<double>
        {
          4.12209020652996E-05,
          -1149872.38280587,
          9481808850.3208,
          -1.95788865718971E+17,
          4.962507048713E+24,
          -1.05549884548496E+28,
          -758642165988.278,
          -9.22172769596101E+22,
          7.25379072059348E+29,
          -61.7718249205859,
          10755.5033344858,
          -37954580.2336487,
          228646846221.831,
          -4997410.93010619,
          -2.80214310054101E+30,
          1049154.06769586,
          6.13754229168619E+27,
          8.02056715528378E+31,
          -29861781.9828065,
          -91.0782540134681,
          135033.227281565,
          -7.12949383408211E+18,
          -1.04578785289542E+36,
          30.4331584444093,
          5932507979.59445,
          -3.64174062110798E+27,
          0.921791403532461,
          -0.337693609657471,
          -72.4644143758508,
          -0.110480239272601,
          5.36516031875059,
          -2914.41872156205,
          6.16338176535305E+39,
          -1.2088917586118E+38,
          8.18396024524612E+22,
          940781944.835829,
          -36727.9669545448,
          -8.37513931798655E+15
        }},
      {
        "h",
        new List<double>
        {
          0.0561379678887577,
          7741354215.87083,
          1.11482975877938E-09,
          -0.00143987128208183,
          1936.9655876492,
          -605971823.585005,
          17195156812433.7,
          -1.85461154985145E+16,
          3.8785116807801E-17,
          -3.95464327846105E-14,
          -170.875935679023,
          -2120.1062070122,
          17768333.7348191,
          11.0177443629575,
          -234396.091693313,
          -6561744.21999594,
          1.56362212977396E-05,
          -2.129462570214,
          13.5249306374858,
          0.177189164145813,
          1394.99167345464,
          -0.00703670932036388,
          -0.152011044389648,
          9.81916922991113E-05,
          0.00147199658618076,
          20.2618487025578,
          0.89934551894424,
          -0.211346402240858,
          24.9971752957491
        }},
      {
        "i",
        new List<double>
        {
          1.06905684359136,
          -1.48620857922333,
          259862256980408.0,
          -4.46352055678749E-12,
          -5.66620757170032E-07,
          -0.00235302885736849,
          -0.269226321968839,
          9.22024992944392,
          3.57633505503772E-12,
          -17.3942565562222,
          7.00681785556229E-06,
          -0.000267050351075768,
          -2.31779669675624,
          -7.53533046979752E-13,
          4.81337131452891,
          -2.23286270422356E+21,
          -1.18746004987383E-05,
          0.00646412934136496,
          -4.10588536330937E-10,
          4.22739537057241E+19,
          3.13698180473812E-13,
          1.6439533434504E-24,
          -3.39823323754373E-06,
          -0.0135268639905021,
          -7.23252514211625E-15,
          1.84386437538366E-09,
          -0.0463959533752385,
          -99226310037675.0,
          6.88169154439335E-17,
          -2.22620998452197E-11,
          -5.40843018624083E-08,
          0.00345570606200257,
          42227580030.4086,
          -1.26974478770487E-15,
          9.27237985153679E-10,
          6.12670812016489E-14,
          -7.22693924063497E-12,
          -0.000383669502636822,
          0.000374684572410204,
          -93197.6897511086,
          -0.0247690616026922,
          65.8110546759474
        }},
      {
        "j",
        new List<double>
        {
          -0.00011137131739554,
          1.00342892423685,
          5.30615581928979,
          1.79058760078792E-06,
          -0.000728541958464774,
          -18.7576133371704,
          0.00199060874071849,
          24.357475537729,
          -0.000177040785499444,
          -0.0025968038522713,
          -198.704578406823,
          7.38627790224287E-05,
          -0.00236264692844138,
          -1.61023121314333,
          6223.22971786473,
          -9.60754116701669E-09,
          -5.10572269720488E-11,
          0.00767373781404211,
          6.63855469485254E-15,
          -7.17590735526745E-10,
          1.46564542926508E-05,
          3.09029474277013E-12,
          -4.64216300971708E-16,
          -3.90499637961161E-14,
          -2.36716126781431E-10,
          4.54652854268717E-12,
          -0.00422271787482497,
          2.83911742354706E-11,
          2.70929002720228
        }},
      {
        "k",
        new List<double>
        {
          -401215699.576099,
          48450147831.8406,
          3.94721471363678E-15,
          37262.9967374147,
          -3.69794374168666E-30,
          -3.80436407012452E-15,
          4.75361629970233E-07,
          -0.000879148916140706,
          0.844317863844331,
          12.24331626566,
          -104.529634830279,
          589.702771277429,
          -29102685116444.4,
          1.7034307284185E-06,
          -0.000277617606975748,
          -3.44709605486686,
          22.1333862447095,
          -194.646110037079,
          8.08354639772825E-16,
          -1.8084520914547E-11,
          -6.96664158132412E-06,
          -0.00181057560300994,
          2.55830298579027,
          3289.13873658481,
          -1.73270241249904E-19,
          -6.61876792558034E-07,
          -0.0039568892342125,
          6.04203299819132E-18,
          -4.00879935920517E-14,
          1.60751107464958E-09,
          3.83719409025556E-05,
          -6.49565446702457E-15,
          -1.49095328506E-12,
          5.41449377329581E-09
        }},
      {
        "l",
        new List<double>
        {
          2607020586.47537,
          -188277213604704.0,
          5.54923870289667E+18,
          -7.58966946387758E+22,
          4.13865186848908E+26,
          -815038000738.06,
          -3.81458260489955E+32,
          -0.0123239564600519,
          22609563.1437174,
          -495017809506.72,
          5.29482996422863E+15,
          -4.44359478746295E+22,
          5.21635864527315E+34,
          -4.87095672740742E+54,
          -714430.209937547,
          0.127868634615495,
          -10.0752127917598,
          7774514.3796099,
          -1.08105480796471E+24,
          -3.57578581169659E-06,
          -2.12857169423484,
          2.70706111085238E+29,
          -6.95953622348829E+32,
          0.11060902747228,
          72.1559163361354,
          -306367307532219.0,
          2.6583961888553E-05,
          0.0253392392889754,
          -214.443041836579,
          0.937846601489667,
          2.231840431017,
          33.8401222509191,
          4.94237237179718E+20,
          -0.198068404154428,
          -1.4141534988114E+30,
          -99.3862421613651,
          125.070534142731,
          -996.473529004439,
          47313.7909872765,
          1.16662121219322E+32,
          -3.15874976271533E+15,
          -4.45703369196945E+32,
          6.42794932373694E+32
        }},
      {
        "m",
        new List<double>
        {
          0.811384363481847,
          -5681.99310990094,
          -17865719817.2556,
          7.95537657613427E+31,
          -81456.8209346872,
          -65977456.7602874,
          -15286114865.9302,
          -560165667510.446,
          458384.828593949,
          -38575400038384.8,
          45373580.0004273,
          939454935735.563,
          2.66572856432938E+27,
          -5475783138.99097,
          200725701112386.0,
          1850072455632.39,
          185135446.828337,
          -170451090076.385,
          157890366037614.0,
          -2.02530509748774E+15,
          3.6819392618357E+59,
          1.70215539458936E+17,
          6.39234909918741E+41,
          -821698160721956.0,
          -7.95260241872306E+23,
          2.3341586947851E+17,
          -6.00079934586803E+22,
          5.94584382273384E+24,
          1.89461279349492E+39,
          -8.10093428842645E+45,
          1.88813911076809E+21,
          1.11052244098768E+35,
          2.91133958602503E+45,
          -3.2942192395146E+21,
          -1.37570282536696E+25,
          1.81508996303902E+27,
          -3.46865122768353E+29,
          -2.1196114877426E+37,
          -1.28617899887675E+48,
          4.79817895699239E+64
        }},
      {
        "n",
        new List<double>
        {
          2.80967799943151E-39,
          6.14869006573609E-31,
          5.82238667048942E-28,
          3.90628369238462E-23,
          8.21445758255119E-21,
          4.02137961842776E-15,
          6.51718171878301E-13,
          -2.11773355803058E-08,
          0.00264953354380072,
          -1.35031446451331E-32,
          -6.07246643970893E-24,
          -4.02352115234494E-19,
          -7.44938506925544E-17,
          1.89917206526237E-13,
          3.64975183508473E-06,
          1.77274872361946E-26,
          -3.34952758812999E-19,
          -4.21537726098389E-09,
          -0.0391048167929649,
          5.41276911564176E-14,
          7.05412100773699E-12,
          2.58585887897486E-09,
          -4.93111362030162E-11,
          -1.58649699894543E-06,
          -0.5250374278861,
          0.00220019901729615,
          -0.00643064132636925,
          62.9154149015048,
          135.147318617061,
          2.40560808321713E-07,
          -0.000890763306701305,
          -4402.09599407714,
          -302.807107747776,
          1591.58748314599,
          232534.272709876,
          -792681.2071326,
          -86987136466.2769,
          354542769185.671,
          400849240129329.0
        }},
      {
        "o",
        new List<double>
        {
          1.28746023979718E-35,
          -7.35234770382342E-12,
          0.0028907869214915,
          0.244482731907223,
          1.41733492030985E-24,
          -3.54533853059476E-29,
          -5.94539202901431E-18,
          -5.85188401782779E-09,
          2.01377325411803E-06,
          1.38647388209306,
          -1.73959365084772E-05,
          0.00137680878349369,
          8.14897605805513E-15,
          4.25596631351839E-26,
          -3.87449113787755E-18,
          1.3981474793024E-13,
          -0.00171849638951521,
          6.41890529513296E-22,
          1.18960578072018E-11,
          -1.55282762571611E-18,
          2.33907907347507E-08,
          -1.74093247766213E-13,
          3.77682649089149E-09,
          -5.16720236575302E-11
        }},
      {
        "p",
        new List<double>
        {
          -9.82825342010366E-05,
          1.05145700850612,
          116.033094095084,
          3246.64750281543,
          -1235.92348610137,
          -0.0561403450013495,
          8.56677401640869E-08,
          236.313425393924,
          0.00972503292350109,
          -1.03001994531927,
          -1.49653706199162E-09,
          -2.15743778861592E-05,
          -8.34452198291445,
          0.586602660564988,
          3.43480022104968E-26,
          8.16256095947021E-06,
          0.00294985697916798,
          7.11730466276584E-17,
          4.00954763806941E-10,
          10.7766027032853,
          -4.09449599138182E-07,
          -7.29121307758902E-06,
          6.77107970938909E-09,
          6.02745973022975E-08,
          -3.82323011855257E-11,
          0.00179946628317437,
          -0.000345042834640005
        }},
      {
        "q",
        new List<double>
        {
          -82043.384325995,
          47327151846.1586,
          -0.0805950021005413,
          32.860002543598,
          -3566.1702998249,
          -1729857814.33335,
          35176923.2729192,
          -775489.259985144,
          7.10346691966018E-05,
          99349.9883820274,
          -0.64209417190457,
          -6128.42816820083,
          232.808472983776,
          -1.42808220416837E-05,
          -0.00643596060678456,
          -4.28577227475614,
          2256.89939161918,
          0.0010035565172151,
          0.333491455143516,
          1.09697576888873,
          0.961917379376452,
          -0.0838165632204598,
          2.47795908411492,
          -3191.14969006533
        }},
      {
        "r",
        new List<double>
        {
          0.00144165955660863,
          -7014385996282.58,
          -8.30946716459219E-17,
          0.261975135368109,
          393.097214706245,
          -10433.4030654021,
          490112654.154211,
          -0.000147104222772069,
          1.03602748043408,
          3.05308890065089,
          -3997452.76971264,
          5.6923371959375E-12,
          -0.0464923504407778,
          -5.35400396512906E-18,
          3.99988795693162E-13,
          -5.36479560201811E-07,
          0.0159536722411202,
          2.70303248860217E-15,
          2.44247453858506E-08,
          -9.83430636716454E-06,
          0.0663513144224454,
          -9.93456957845006,
          546.491323528491,
          -14336.5406393758,
          150764.974125511,
          -3.37209709340105E-10,
          3.77501980025469E-09
        }},
      {
        "s",
        new List<double>
        {
          -5.32466612140254E+22,
          1.00415480000824E+31,
          -1.91540001821367E+29,
          1.05618377808847E+16,
          2.02281884477061E+58,
          88458547.2596134,
          1.66540181638363E+22,
          -313563.197669111,
          -1.85662327545324E+53,
          -0.0624942093918942,
          -5041607241.3259,
          18751.4491833092,
          0.00121399979993217,
          1.88317043049455,
          -1670.7350396206,
          0.965961650599775,
          2.94885696802488,
          -65391.5627346115,
          6.04012200163444E+49,
          -0.198339358557937,
          -1.75984090163501E+57,
          3.56314881403987,
          -575.991255144384,
          45621.3415338071,
          -10917404.4987829,
          4.37796099975134E+33,
          -6.16552611135792E+45,
          1935687689.17797,
          9.50898170425042E+53
        }},
      {
        "t",
        new List<double>
        {
          1.55287249586268,
          6.64235115009031,
          -2893.6623672721,
          -3859232023098.48,
          -2.91002915783761,
          -829088246858.083,
          1.76814899675218,
          -534686695.713469,
          1.60464608687834E+17,
          196435.366560186,
          1566374275417.29,
          -1.78154560260006,
          -2.29746237623692E+15,
          38565900.1648006,
          1105544467.90543,
          -67707383068734.9,
          -3.27910592086523E+30,
          -3.41552040860644E+50,
          -5.27251339709047E+20,
          2.45375640937055E+23,
          -1.68776617209269E+26,
          3.58958955867578E+28,
          -6.56475280339411E+35,
          3.55286045512301E+38,
          5.6902145441327E+57,
          -7.00584546433113E+47,
          -7.05772623326374E+64,
          1.66861176200148E+52,
          -3.00475129680486E+60,
          -6.68481295196808E+50,
          4.28432338620678E+68,
          -4.44227367758304E+71,
          -2.81396013562745E+76
        }},
      {
        "u",
        new List<double>
        {
          1.22088349258355E+17,
          1042164686.08488,
          -8.82666931564652E+15,
          2.59929510849499E+19,
          222612779142211.0,
          -8.78473585050085E+17,
          -3.14432577551552E+21,
          -2169349169962.85,
          1.59079648196849E+20,
          -339.567617303423,
          8843876513378.36,
          -8.43405926846418E+20,
          11.4178193518022,
          -0.000122708229235641,
          -106.201671767107,
          9.03443213959313E+24,
          -6.93996270370852E+27,
          6.48916718965575E-09,
          7189.57567127851,
          0.00105581745346187,
          -651903203602581.0,
          -1.60116813274676E+24,
          -5.10254294237837E-09,
          -0.152355388953402,
          677143292290.144,
          276378438378930.0,
          0.0116862983141686,
          -30142694798017.1,
          1.6971981388484E-08,
          1.04674840020929E+26,
          -10801.690456014,
          -9.90623601934295E-13,
          5361164.83602738,
          2.26145963747881E+21,
          -4.8873156577621E-10,
          1.5100154888067E-05,
          -22770.046464392,
          -7.81754507698846E+27
        }},
      {
        "v",
        new List<double>
        {
          -4.15652812061591E-55,
          1.77441742924043E-61,
          -3.57078668203377E-55,
          3.59252213604114E-26,
          -25.9123736380269,
          59461.976619346,
          -62418400710.3158,
          3.13080299915944E+16,
          1.05006446192036E-09,
          -1.92824336984852E-06,
          654144.373749937,
          5131174628650.44,
          -6.97595750347391E+18,
          -1.03977184454767E+28,
          1.19563135540666E-48,
          -4.36677034051655E-42,
          9.26990036530639E-30,
          5.87793105620748E+20,
          2.80375725094731E-18,
          -1.92359972440634E+22,
          7.42705723302738E+26,
          -51.7429682450605,
          8206120.48645469,
          -1.88214882341448E-09,
          0.0184587261114837,
          -1.35830407782663E-06,
          -7.23681885626348E+16,
          -2.23449194054124E+26,
          -1.11526741826431E-35,
          2.76032601145151E-29,
          134856491567853.0,
          6.5244029334586E-10,
          5.1065511977436E+16,
          -4.68138358908732E+31,
          -7.60667491183279E+15,
          -4.17247986986821E-19,
          31254567775610.4,
          -100375333864186.0,
          2.47761392329058E+26
        }},
      {
        "w",
        new List<double>
        {
          -5.86219133817016E-08,
          -89446035500.5526,
          5.31168037519774E-31,
          0.109892402329239,
          -0.0575368389425212,
          22827.6853990249,
          -1.58548609655002E+18,
          3.29865748576503E-28,
          -6.34987981190669E-25,
          6.15762068640611E-09,
          -96110924.0985747,
          -4.06274286652625E-45,
          -4.71103725498077E-13,
          0.725937724828145,
          1.87768525763682E-39,
          -1033.08436323771,
          -0.0662552816342168,
          579.51404176571,
          2.37416732616644E-27,
          2.71700235739893E-15,
          -90.78862134836,
          -1.71242509570207E-37,
          156.792067854621,
          0.92326135790147,
          -5.97865988422577,
          3219887.67636389,
          -3.99441390042203E-30,
          4.93429086046981E-08,
          8.12036983370565E-20,
          -2.07610284654137E-12,
          -3.40821291419719E-07,
          5.42000573372233E-18,
          -8.56711586510214E-13,
          2.66170454405981E-14,
          8.58133791857099E-06
        }},
      {
        "x",
        new List<double>
        {
          3.77373741298151E+18,
          -5071008837229.13,
          -1.0336322559886E+15,
          1.84790814320773E-06,
          -0.000924729378390945,
          -4.25999562292738E+23,
          -4.62307771873973E-13,
          1.07319065855767E+21,
          64866249228.0682,
          2.44200600688281,
          -8515357334.84258,
          1.69894481433592E+21,
          2.1578022250902E-27,
          -0.320850551367334,
          -3.8264244845861E+16,
          -2.75386077674421E-29,
          -563199.253391666,
          -3.26068646279314E+20,
          39794900155318.4,
          1.00824008584757E-07,
          16223.4569738433,
          -43235522531.9745,
          -592874245598.61,
          1.33061647281106,
          1573381.97797544,
          25818961427085.3,
          2.62413209706358E+24,
          -0.0920011937431142,
          0.00220213765905426,
          -11.0433759109547,
          8470048.70612087,
          -592910695.762536,
          -1.8302717326966E-05,
          0.181339603516302,
          -1192.28759669889,
          4308676.58061468
        }},
      {
        "y",
        new List<double>
        {
          -5.25597995024633E-10,
          5834.41305228407,
          -1.34778968457925E+16,
          1.18973500934212E+25,
          -1.59096490904708E+26,
          -3.15839902302021E-07,
          496.212197158239,
          3.27777227273171E+18,
          -5.27114657850696E+21,
          2.10017506281863E-17,
          7.05106224399834E+20,
          -2.66713136106469E+30,
          -1.45370512554562E-08,
          1.4933391705313E+27,
          -14979562.0287641,
          -3.818819062711E+15,
          7.24660165585797E-05,
          -93780816955019.3,
          5144114683.76383,
          -82819.8594040141
        }},
      {
        "z",
        new List<double>
        {
          2.4400789229065E-11,
          -4630574.30331242,
          7288032747.77712,
          3.27776302858856E+15,
          -1105981701.18409,
          -3238999157299.57,
          9.23814007023245E+15,
          8.42250080413712E-13,
          663221436245.506,
          -167170186672139.0,
          2537.49358701391,
          -8.19731559610523E-21,
          328380587890.663,
          -62500479.1171543,
          8.03197957462023E+20,
          -2.04397011338353E-11,
          -3783.91047055938,
          0.0097287654593862,
          15.4355721681459,
          -3739.62862928643,
          -68285901137.4572,
          -0.000248488015614543,
          3945360.49497068
        }}};
    I = I[x];
    J = J[x];
    n = n[x];
    var _tup_1 = par[x];
    var v_ = _tup_1.Item1;
    var P_ = _tup_1.Item2;
    var T_ = _tup_1.Item3;
    var a = _tup_1.Item4;
    var b = _tup_1.Item5;
    var c = _tup_1.Item6;
    var d = _tup_1.Item7;
    var e = _tup_1.Item8;
    var Pr = P / P_;
    var Tr = T / T_;
    var suma = 0;
    if(x == "n")
    {
      foreach (var _tup_2 in zip(I, J, n))
      {
        i = _tup_2.Item1;
        j = _tup_2.Item2;
        ni = _tup_2.Item3;
        suma += ni * Math.Pow(Pr - a, i) * Math.Pow(Tr - b, j);
      }
      return v_ * Math.Exp(suma);
    }
    else
    {
      foreach (var _tup_3 in zip(I, J, n))
      {
        i = _tup_3.Item1;
        j = _tup_3.Item2;
        ni = _tup_3.Item3;
        suma += ni * Math.Pow(Pr - a, c * i) * Math.Pow(Tr - b, j * d);
      }
      return v_ * Math.Pow(suma, e);
    }
  }

  // Region 4
  // Basic equation for region 4
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //   x : float
  //     Vapor quality, [-]
  //
  //   Returns
  //   -------
  //   prop : dict
  //     Dict with calculated properties. The available properties are:
  //
  //       * T: Saturated temperature, [K]
  //       * P: Saturated pressure, [MPa]
  //       * x: Vapor quality, [-]
  //       * v: Specific volume, [m³/kg]
  //       * h: Specific enthalpy, [kJ/kg]
  //       * s: Specific entropy, [kJ/kgK]
  //

  public static Dictionary<string, double> _Region4(double P, double x)
  {
    object P2;
    object P1;
    var T = _TSat_P(P);
    if(T > 623.15)
    {
      var rhol = 1.0 / _Backward3_sat_v_P(P, T, 0);
      P1 = _Region3(rhol, T);
      var rhov = 1.0 / _Backward3_sat_v_P(P, T, 1);
      P2 = _Region3(rhov, T);
    }
    else
    {
      P1 = _Region1(T, P);
      P2 = _Region2(T, P);
    }

    var propiedades = new Dictionary<string, double> { };
    propiedades["T"] = T;
    propiedades["P"] = P;
    propiedades["v"] = P1["v"] + x * (P2["v"] - P1["v"]);
    propiedades["h"] = P1["h"] + x * (P2["h"] - P1["h"]);
    propiedades["s"] = P1["s"] + x * (P2["s"] - P1["s"]);
    propiedades["cp"] = null;
    propiedades["cv"] = null;
    propiedades["w"] = null;
    propiedades["alfav"] = null;
    propiedades["kt"] = null;
    propiedades["region"] = 4;
    propiedades["x"] = x;
    return propiedades;
  }

  // Backward equation for region 4, T=f(h,s)
  //
  //   Parameters
  //   ----------
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //   s : float
  //     Specific entropy, [kJ/kgK]
  //
  //   Returns
  //   -------
  //   T : float
  //     Temperature, [K]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Supplementary Release on Backward Equations p(h,s) for
  //   Region 3, Equations as a Function of h and s for the Region Boundaries, and
  //   an Equation Tsat(h,s) for Region 4 of the IAPWS Industrial Formulation 1997
  //   for the Thermodynamic Properties of Water and Steam,
  //   http://www.iapws.org/relguide/Supp-phs3-2014.pdf. Eq 9
  //
  //   Examples
  //   --------
  //   >>> _Backward4_T_hs(1800,5.3)
  //   346.8475498
  //   >>> _Backward4_T_hs(2400,6.0)
  //   425.1373305
  //   >>> _Backward4_T_hs(2500,5.5)
  //   522.5579013
  //

  public static object _Backward4_T_hs(object h, object s)
  {
    var I = new List<int>
    {
      0,
      0,
      0,
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
      3,
      4,
      4,
      5,
      5,
      5,
      5,
      6,
      6,
      6,
      8,
      10,
      10,
      12,
      14,
      14,
      16,
      16,
      18,
      18,
      18,
      20,
      28
    };
    var J = new List<int>
    {
      0,
      3,
      12,
      0,
      1,
      2,
      5,
      0,
      5,
      8,
      0,
      2,
      3,
      4,
      0,
      1,
      1,
      2,
      4,
      16,
      6,
      8,
      22,
      1,
      20,
      36,
      24,
      1,
      28,
      12,
      32,
      14,
      22,
      36,
      24,
      36
    };
    var n = new List<double>
    {
      0.179882673606601,
      -0.267507455199603,
      1.162767226126,
      0.147545428713616,
      -0.512871635973248,
      0.421333567697984,
      0.56374952218987,
      0.429274443819153,
      -3.3570455214214,
      10.8890916499278,
      -0.248483390456012,
      0.30415322190639,
      -0.494819763939905,
      1.07551674933261,
      0.0733888415457688,
      0.0140170545411085,
      -0.106110975998808,
      0.0168324361811875,
      1.25028363714877,
      1013.16840309509,
      -1.51791558000712,
      52.4277865990866,
      23049.5545563912,
      0.0249459806365456,
      2107964.67412137,
      366836848.613065,
      -144814105.365163,
      -0.0017927637300359,
      4899556021.00459,
      471.262212070518,
      -82929439019.8652,
      -1715.45662263191,
      3557776.82973575,
      586062760258.436,
      -12988763.5078195,
      31724744937.1057
    };
    var nu = h / 2800;
    var sigma = s / 9.2;
    var suma = 0;
    foreach (var _tup_1 in zip(I, J, n))
    {
      var i = _tup_1.Item1;
      var j = _tup_1.Item2;
      var ni = _tup_1.Item3;
      suma += ni * Math.Pow(nu - 0.119, i) * Math.Pow(sigma - 1.07, j);
    }
    return 550 * suma;
  }

  // Region 5
  // Basic equation for region 5
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
  //     Dict with calculated properties. The available properties are:
  //
  //       * v: Specific volume, [m³/kg]
  //       * h: Specific enthalpy, [kJ/kg]
  //       * s: Specific entropy, [kJ/kgK]
  //       * cp: Specific isobaric heat capacity, [kJ/kgK]
  //       * cv: Specific isocoric heat capacity, [kJ/kgK]
  //       * w: Speed of sound, [m/s]
  //       * alfav: Cubic expansion coefficient, [1/K]
  //       * kt: Isothermal compressibility, [1/MPa]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  //   Thermodynamic Properties of Water and Steam August 2007,
  //   http://www.iapws.org/relguide/IF97-Rev.html, Eq 32-34
  //
  //   Examples
  //   --------
  //   >>> _Region5(1500,0.5)["v"]
  //   1.38455090
  //   >>> _Region5(1500,0.5)["h"]
  //   5219.76855
  //   >>> _Region5(1500,0.5)["h"]-500*_Region5(1500,0.5)["v"]
  //   4527.49310
  //   >>> _Region5(1500,30)["s"]
  //   7.72970133
  //   >>> _Region5(1500,30)["cp"]
  //   2.72724317
  //   >>> _Region5(1500,30)["cv"]
  //   2.19274829
  //   >>> _Region5(2000,30)["w"]
  //   1067.36948
  //   >>> _Region5(2000,30)["alfav"]
  //   0.000508830641
  //   >>> _Region5(2000,30)["kt"]
  //   0.0329193892
  //

  public static Dictionary<string, double> _Region5(double T, double P)
  {
    if(P < 0)
    {
      P = Pmin;
    }
    var Tr = 1000 / T;
    var Pr = P / 1;
    var _tup_1 = Region5_cp0(Tr, Pr);
    var go = _tup_1.Item1;
    var gop = _tup_1.Item2;
    var gopp = _tup_1.Item3;
    var got = _tup_1.Item4;
    var gott = _tup_1.Item5;
    var gopt = _tup_1.Item6;
    var Ir = new List<int>
    {
      1,
      1,
      1,
      2,
      2,
      3
    };
    var Jr = new List<int>
    {
      1,
      2,
      3,
      3,
      9,
      7
    };
    var nr = new List<double>
    {
      0.0015736404855259,
      0.00090153761673944,
      -0.0050270077677648,
      2.2440037409485E-06,
      -4.1163275453471E-06,
      3.7919454822955E-08
    };
    var gr = 0;
    foreach (var _tup_2 in zip(Ir, Jr, nr))
    {
      var i = _tup_2.Item1;
      var j = _tup_2.Item2;
      var ni = _tup_2.Item3;
      gr += ni * Math.Pow(Pr, i) * Math.Pow(Tr, j);
      grp += ni * i * Math.Pow(Pr, i - 1) * Math.Pow(Tr, j);
      grpp += ni * i * (i - 1) * Math.Pow(Pr, i - 2) * Math.Pow(Tr, j);
      grt += ni * j * Math.Pow(Pr, i) * Math.Pow(Tr, j - 1);
      grtt += ni * j * (j - 1) * Math.Pow(Pr, i) * Math.Pow(Tr, j - 2);
      grpt += ni * i * j * Math.Pow(Pr, i - 1) * Math.Pow(Tr, j - 1);
    }

    var propiedades = new Dictionary<string, double> { };
    propiedades["T"] = T;
    propiedades["P"] = P;
    propiedades["v"] = Pr * (gop + grp) * R * T / P / 1000;
    propiedades["h"] = Tr * (got + grt) * R * T;
    propiedades["s"] = R * (Tr * (got + grt) - (go + gr));
    propiedades["cp"] = -R * Math.Pow(Tr, 2) * (gott + grtt);
    propiedades["cv"] = R * (-Math.Pow(Tr, 2) * (gott + grtt) + Math.Pow(gop + grp - Tr * (gopt + grpt), 2) / (gopp + grpp));
    propiedades["w"] = Math.Pow(R * T * 1000 * (1 + 2 * Pr * grp + Math.Pow(Pr, 2) * Math.Pow(grp, 2)) / (1 - Math.Pow(Pr, 2) * grpp + Math.Pow(1 + Pr * grp - Tr * Pr * grpt, 2) / Math.Pow(Tr, 2) / (gott + grtt)), 0.5);
    propiedades["alfav"] = (1 + Pr * grp - Tr * Pr * grpt) / (1 + Pr * grp) / T;
    propiedades["kt"] = (1 - Math.Pow(Pr, 2) * grpp) / (1 + Pr * grp) / P;
    propiedades["region"] = 5;
    propiedades["x"] = 1;
    return propiedades;
  }

  // Ideal properties for Region 5
  //
  //   Parameters
  //   ----------
  //   Tr : float
  //     Reduced temperature, [-]
  //   Pr : float
  //     Reduced pressure, [-]
  //
  //   Returns
  //   -------
  //   prop : array
  //     Array with ideal Gibbs energy partial derivatives:
  //
  //       * g: Ideal Specific Gibbs energy, [kJ/kg]
  //       * gp: [∂g/∂P]T
  //       * gpp: [∂²g/∂P²]T
  //       * gt: [∂g/∂T]P
  //       * gtt: [∂²g/∂T²]P
  //       * gpt: [∂²g/∂T∂P]
  //
  //   References
  //   ----------
  //   IAPWS, Revised Release on the IAPWS Industrial Formulation 1997 for the
  //   Thermodynamic Properties of Water and Steam August 2007,
  //   http://www.iapws.org/relguide/IF97-Rev.html, Eq 33
  //

  public static object Region5_cp0(object Tr, object Pr)
  {
    var Jo = new List<int> { 0, 1, -3, -2, -1, 2 };
    var no = new List<double> { -13.179983674201, 6.8540841634434, -0.024805148933466, 0.36901534980333, -3.1161318213925, -0.32961626538917 };
    var go = Math.Log(Pr);
    var gop = Math.Pow(Pr, -1);
    var gopp = -Math.Pow(Pr, -2);
    var got = 0;
    foreach (var _tup_1 in zip(Jo, no))
    {
      var j = _tup_1.Item1;
      var ni = _tup_1.Item2;
      go += ni * Math.Pow(Tr, j);
      got += ni * j * Math.Pow(Tr, j - 1);
      gott += ni * j * (j - 1) * Math.Pow(Tr, j - 2);
    }
    return Tuple.Create(go, gop, gopp, got, gott, gopt);
  }

  // Region definitions
  // Region definition for input T and P
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
  //   region : float
  //     IAPWS-97 region code
  //
  //   References
  //   ----------
  //   Wagner, W; Kretzschmar, H-J: International Steam Tables: Properties of
  //   Water and Steam Based on the Industrial Formulation IAPWS-IF97; Springer,
  //   2008; doi: 10.1007/978-3-540-74234-0. Fig. 2.3
  //

  public static object _Bound_TP(object T, object P)
  {
    object region = null;
    if(1073.15 < T <= 2273.15 && Pmin <= P <= 50) { region = 5; }
    else if(Pmin <= P <= Ps_623)
    {
      var Tsat = _TSat_P(P);
      if(273.15 <= T <= Tsat) { region = 1; }
      else if(Tsat < T <= 1073.15) { region = 2; }
    }
    else if(Ps_623 < P <= 100)
    {
      var T_b23 = _t_P(P);
      if(273.15 <= T <= 623.15) { region = 1; }
      else if(623.15 < T < T_b23) { region = 3; }
      else if(T_b23 <= T <= 1073.15) { region = 2; }
    }
    return region;
  }

  // Region definition for input P y h
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //
  //   Returns
  //   -------
  //   region : float
  //     IAPWS-97 region code
  //
  //   References
  //   ----------
  //   Wagner, W; Kretzschmar, H-J: International Steam Tables: Properties of
  //   Water and Steam Based on the Industrial Formulation IAPWS-IF97; Springer,
  //   2008; doi: 10.1007/978-3-540-74234-0. Fig. 2.5
  //

  public static object _Bound_Ph(object P, object h)
  {
    object p34;
    object h32;
    object h13;
    object region = null;
    if(Pmin <= P <= Ps_623)
    {
      var h14 = _Region1(_TSat_P(P), P)["h"];
      var h24 = _Region2(_TSat_P(P), P)["h"];
      var h25 = _Region2(1073.15, P)["h"];
      var hmin = _Region1(273.15, P)["h"];
      var hmax = _Region5(2273.15, P)["h"];
      if(hmin <= h <= h14) { region = 1; }
      else if(h14 < h < h24) { region = 4; }
      else if(h24 <= h <= h25) { region = 2; }
      else if(h25 < h <= hmax) { region = 5; }
    }
    else if(Ps_623 < P < Pc)
    {
      hmin = _Region1(273.15, P)["h"];
      h13 = _Region1(623.15, P)["h"];
      h32 = _Region2(_t_P(P), P)["h"];
      h25 = _Region2(1073.15, P)["h"];
      hmax = _Region5(2273.15, P)["h"];
      if(hmin <= h <= h13) { region = 1; }
      else if(h13 < h < h32)
      {
        try
        {
          p34 = _PSat_h(h);
        } catch (NotImplementedError)
        {
          p34 = Pc;
        }
        if(P < p34) { region = 4; }
        else { region = 3; }
      }
      else if(h32 <= h <= h25) { region = 2; }
      else if(h25 < h <= hmax) { region = 5; }
    }
    else if(Pc <= P <= 100)
    {
      hmin = _Region1(273.15, P)["h"];
      h13 = _Region1(623.15, P)["h"];
      h32 = _Region2(_t_P(P), P)["h"];
      h25 = _Region2(1073.15, P)["h"];
      hmax = _Region5(2273.15, P)["h"];
      if(hmin <= h <= h13) { region = 1; }
      else if(h13 < h < h32) { region = 3; }
      else if(h32 <= h <= h25) { region = 2; }
      else if(P <= 50 && h25 <= h <= hmax) { region = 5; }
    }
    return region;
  }

  // Region definition for input P and s
  //
  //   Parameters
  //   ----------
  //   P : float
  //     Pressure, [MPa]
  //   s : float
  //     Specific entropy, [kJ/kgK]
  //
  //   Returns
  //   -------
  //   region : float
  //     IAPWS-97 region code
  //
  //   References
  //   ----------
  //   Wagner, W; Kretzschmar, H-J: International Steam Tables: Properties of
  //   Water and Steam Based on the Industrial Formulation IAPWS-IF97; Springer,
  //   2008; doi: 10.1007/978-3-540-74234-0. Fig. 2.9
  //

  public static object _Bound_Ps(object P, object s)
  {
    object p34;
    object s32;
    object s13;
    object region = null;
    if(Pmin <= P <= Ps_623)
    {
      var smin = _Region1(273.15, P)["s"];
      var s14 = _Region1(_TSat_P(P), P)["s"];
      var s24 = _Region2(_TSat_P(P), P)["s"];
      var s25 = _Region2(1073.15, P)["s"];
      var smax = _Region5(2273.15, P)["s"];
      if(smin <= s <= s14) { region = 1; }
      else if(s14 < s < s24) { region = 4; }
      else if(s24 <= s <= s25) { region = 2; }
      else if(s25 < s <= smax) { region = 5; }
    }
    else if(Ps_623 < P < Pc)
    {
      smin = _Region1(273.15, P)["s"];
      s13 = _Region1(623.15, P)["s"];
      s32 = _Region2(_t_P(P), P)["s"];
      s25 = _Region2(1073.15, P)["s"];
      smax = _Region5(2273.15, P)["s"];
      if(smin <= s <= s13) { region = 1; }
      else if(s13 < s < s32)
      {
        try
        {
          p34 = _PSat_s(s);
        } catch (NotImplementedError)
        {
          p34 = Pc;
        }
        if(P < p34) { region = 4; }
        else { region = 3; }
      }
      else if(s32 <= s <= s25) { region = 2; }
      else if(s25 < s <= smax) { region = 5; }
    }
    else if(Pc <= P <= 100)
    {
      smin = _Region1(273.15, P)["s"];
      s13 = _Region1(623.15, P)["s"];
      s32 = _Region2(_t_P(P), P)["s"];
      s25 = _Region2(1073.15, P)["s"];
      smax = _Region5(2273.15, P)["s"];
      if(smin <= s <= s13) { region = 1; }
      else if(s13 < s < s32) { region = 3; }
      else if(s32 <= s <= s25) { region = 2; }
      else if(P <= 50 && s25 <= s <= smax) { region = 5; }
    }
    return region;
  }

  // Region definition for input h and s
  //
  //   Parameters
  //   ----------
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //   s : float
  //     Specific entropy, [kJ/kgK]
  //
  //   Returns
  //   -------
  //   region : float
  //     IAPWS-97 region code
  //
  //   References
  //   ----------
  //   Wagner, W; Kretzschmar, H-J: International Steam Tables: Properties of
  //   Water and Steam Based on the Industrial Formulation IAPWS-IF97; Springer,
  //   2008; doi: 10.1007/978-3-540-74234-0. Fig. 2.14
  //

  public static object _Bound_hs(object h, object s)
  {
    object P;
    object T;
    object region = null;
    var s13 = _Region1(623.15, 100)["s"];
    var s13s = _Region1(623.15, Ps_623)["s"];
    var sTPmax = _Region2(1073.15, 100)["s"];
    var s2ab = _Region2(1073.15, 4)["s"];
    // Left point in h-s plot
    var smin = _Region1(273.15, 100)["s"];
    var hmin = _Region1(273.15, Pmin)["h"];
    // Right point in h-s plot
    var _Pmax = _Region2(1073.15, Pmin);
    var hmax = _Pmax["h"];
    var smax = _Pmax["s"];
    // Region 4 left and right point
    var _sL = _Region1(273.15, Pmin);
    var h4l = _sL["h"];
    var s4l = _sL["s"];
    var _sV = _Region2(273.15, Pmin);
    var h4v = _sV["h"];
    var s4v = _sV["s"];
    if(smin <= s <= s13)
    {
      hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
      var hs = _h1_s(s);
      T = _Backward1_T_Ps(100, s) - 0.0218;
      hmax = _Region1(T, 100)["h"];
      if(hmin <= h < hs) { region = 4; }
      else if(hs <= h <= hmax) { region = 1; }
    }
    else if(s13 < s <= s13s)
    {
      hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
      hs = _h1_s(s);
      var h13 = _h13_s(s);
      var v = _Backward3_v_Ps(100, s) * (1 + 9.6E-05);
      T = _Backward3_T_Ps(100, s) - 0.0248;
      hmax = _Region3(1 / v, T)["h"];
      if(hmin <= h < hs) { region = 4; }
      else if(hs <= h < h13) { region = 1; }
      else if(h13 <= h <= hmax) { region = 3; }
    }
    else if(s13s < s <= sc)
    {
      hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
      hs = _h3a_s(s);
      v = _Backward3_v_Ps(100, s) * (1 + 9.6E-05);
      T = _Backward3_T_Ps(100, s) - 0.0248;
      hmax = _Region3(1 / v, T)["h"];
      if(hmin <= h < hs) { region = 4; }
      else if(hs <= h <= hmax) { region = 3; }
    }
    else if(sc < s < 5.049096828)
    {
      hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
      hs = _h2c3b_s(s);
      v = _Backward3_v_Ps(100, s) * (1 + 9.6E-05);
      T = _Backward3_T_Ps(100, s) - 0.0248;
      hmax = _Region3(1 / v, T)["h"];
      if(hmin <= h < hs) { region = 4; }
      else if(hs <= h <= hmax) { region = 3; }
    }
    else if(5.049096828 <= s < 5.260578707)
    {
      // Specific zone with 2-3 boundary in s shape
      hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
      hs = _h2c3b_s(s);
      var h23max = _Region2(863.15, 100)["h"];
      var h23min = _Region2(623.15, Ps_623)["h"];
      T = _Backward2_T_Ps(100, s) - 0.019;
      hmax = _Region2(T, 100)["h"];
      if(hmin <= h < hs)
      {
        region = 4;
      }
      else if(hs <= h < h23min)
      {
        region = 3;
      }
      else if(h23min <= h < h23max)
      {
        if(_Backward2c_P_hs(h, s) <= _P23_T(_t_hs(h, s))) { region = 2; }
        else { region = 3; }
      }
      else if(h23max <= h <= hmax) { region = 2; }
    }
    else if(5.260578707 <= s < 5.85)
    {
      hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
      hs = _h2c3b_s(s);
      T = _Backward2_T_Ps(100, s) - 0.019;
      hmax = _Region2(T, 100)["h"];
      if(hmin <= h < hs) { region = 4; }
      else if(hs <= h <= hmax) { region = 2; }
    }
    else if(5.85 <= s < sTPmax)
    {
      hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
      hs = _h2ab_s(s);
      T = _Backward2_T_Ps(100, s) - 0.019;
      hmax = _Region2(T, 100)["h"];
      if(hmin <= h < hs) { region = 4; }
      else if(hs <= h <= hmax) { region = 2; }
    }
    else if(sTPmax <= s < s2ab)
    {
      hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
      hs = _h2ab_s(s);
      P = _Backward2_P_hs(h, s);
      hmax = _Region2(1073.15, P)["h"];
      if(hmin <= h < hs) { region = 4; }
      else if(hs <= h <= hmax) { region = 2; }
    }
    else if(s2ab <= s < s4v)
    {
      hmin = h4l + (s - s4l) / (s4v - s4l) * (h4v - h4l);
      hs = _h2ab_s(s);
      P = _Backward2_P_hs(h, s);
      hmax = _Region2(1073.15, P)["h"];
      if(hmin <= h < hs) { region = 4; }
      else if(hs <= h <= hmax) { region = 2; }
    }
    else if(s4v <= s <= smax)
    {
      hmin = _Region2(273.15, Pmin)["h"];
      P = _Backward2a_P_hs(h, s);
      hmax = _Region2(1073.15, P)["h"];
      if(Pmin <= P <= 100 && hmin <= h <= hmax) { region = 2; }
    }
    // Check region 5
    if(!region && _Region5(1073.15, 50)["s"] < s <= _Region5(2273.15, Pmin)["s"] && _Region5(1073.15, 50)["h"] < h <= _Region5(2273.15, Pmin)["h"])
    {
      var _tup_1 = fsolve(funcion, new List<int>
      {
        1400,
        1
      });
      T = _tup_1.Item1;
      P = _tup_1.Item2;
      if(1073.15 < T <= 2273.15 && Pmin <= P <= 50) { region = 5; }
    }
    Func<object, object> funcion = par =>
    {
      return (_Region5(par[0], par[1])["h"] - h, _Region5(par[0], par[1])["s"] - s);
    };
    return region;
  }

  // Ideal gas properties
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
  //     Dict with calculated properties. The available properties are:
  //
  //       * v: Specific volume, [m³/kg]
  //       * h: Specific enthalpy, [kJ/kg]
  //       * s: Specific entropy, [kJ/kgK]
  //       * cp: Specific isobaric heat capacity, [kJ/kgK]
  //       * cv: Specific isocoric heat capacity, [kJ/kgK]
  //       * w: Speed of sound, [m/s]
  //       * alfav: Cubic expansion coefficient, [1/K]
  //       * kt: Isothermal compressibility, [1/MPa]
  //

  public static object prop0(object T, object P)
  {
    object gopt;
    object gott;
    object got;
    object gopp;
    object gop;
    object go;
    object Pr;
    object Tr;
    if(T <= 1073.15)
    {
      Tr = 540 / T;
      Pr = P / 1.0;
      var _tup_1 = Region2_cp0(Tr, Pr);
      go = _tup_1.Item1;
      gop = _tup_1.Item2;
      gopp = _tup_1.Item3;
      got = _tup_1.Item4;
      gott = _tup_1.Item5;
      gopt = _tup_1.Item6;
    }
    else
    {
      Tr = 1000 / T;
      Pr = P / 1.0;
      var _tup_2 = Region5_cp0(Tr, Pr);
      go = _tup_2.Item1;
      gop = _tup_2.Item2;
      gopp = _tup_2.Item3;
      got = _tup_2.Item4;
      gott = _tup_2.Item5;
      gopt = _tup_2.Item6;
    }
    var prop0 = new Dictionary<object, object>
    {
    };
    prop0["v"] = Pr * gop * R * T / P / 1000;
    prop0["h"] = Tr * got * R * T;
    prop0["s"] = R * (Tr * got - go);
    prop0["cp"] = -R * Math.Pow(Tr, 2) * gott;
    prop0["cv"] = R * (-Math.Pow(Tr, 2) * gott - 1);
    prop0["w"] = Math.Pow(R * T * 1000 / (1 + 1 / Math.Pow(Tr, 2) / gott), 0.5);
    prop0["alfav"] = 1 / T;
    prop0["xkappa"] = 1 / P;
    return prop0;
  }

  // Class to model a state of liquid water or steam with the IAPWS-IF97
  //
  //   Parameters
  //   ----------
  //   T : float
  //     Temperature, [K]
  //   P : float
  //     Pressure, [MPa]
  //   h : float
  //     Specific enthalpy, [kJ/kg]
  //   s : float
  //     Specific entropy, [kJ/kgK]
  //   x : float
  //     Vapor quality, [-]
  //   l : float, optional
  //     Wavelength of light, for refractive index, [nm]
  //
  //   Notes
  //   -----
  //   Definitions options:
  //
  //     * T, P: Not valid for two-phases region
  //     * P, h
  //     * P, s
  //     * h, s
  //     * T, x: Only for two-phases region
  //     * P, x: Only for two-phases region
  //
  //   Returns
  //   -------
  //   prop : dict
  //     The calculated instance has the following properties:
  //
  //       * P: Pressure, [MPa]
  //       * T: Temperature, [K]
  //       * g: Specific Gibbs free energy, [kJ/kg]
  //       * a: Specific Helmholtz free energy, [kJ/kg]
  //       * v: Specific volume, [m³/kg]
  //       * rho: Density, [kg/m³]
  //       * h: Specific enthalpy, [kJ/kg]
  //       * u: Specific internal energy, [kJ/kg]
  //       * s: Specific entropy, [kJ/kg·K]
  //       * cp: Specific isobaric heat capacity, [kJ/kg·K]
  //       * cv: Specific isochoric heat capacity, [kJ/kg·K]
  //       * Z: Compression factor, [-]
  //       * fi: Fugacity coefficient, [-]
  //       * f: Fugacity, [MPa]
  //
  //       * gamma: Isoentropic exponent, [-]
  //       * alfav: Isobaric cubic expansion coefficient, [1/K]
  //       * xkappa: Isothermal compressibility, [1/MPa]
  //       * kappas: Adiabatic compresibility, [1/MPa]
  //       * alfap: Relative pressure coefficient, [1/K]
  //       * betap: Isothermal stress coefficient, [kg/m³]
  //       * joule: Joule-Thomson coefficient, [K/MPa]
  //       * deltat: Isothermal throttling coefficient, [kJ/kg·MPa]
  //       * region: Region
  //
  //       * v0: Ideal specific volume, [m³/kg]
  //       * u0: Ideal specific internal energy, [kJ/kg]
  //       * h0: Ideal specific enthalpy, [kJ/kg]
  //       * s0: Ideal specific entropy, [kJ/kg·K]
  //       * a0: Ideal specific Helmholtz free energy, [kJ/kg]
  //       * g0: Ideal specific Gibbs free energy, [kJ/kg]
  //       * cp0: Ideal specific isobaric heat capacity, [kJ/kg·K]
  //       * cv0: Ideal specific isochoric heat capacity [kJ/kg·K]
  //       * w0: Ideal speed of sound, [m/s]
  //       * gamma0: Ideal isoentropic exponent, [-]
  //
  //       * w: Speed of sound, [m/s]
  //       * mu: Dynamic viscosity, [Pa·s]
  //       * nu: Kinematic viscosity, [m²/s]
  //       * k: Thermal conductivity, [W/m·K]
  //       * alfa: Thermal diffusivity, [m²/s]
  //       * sigma: Surface tension, [N/m]
  //       * epsilon: Dielectric constant, [-]
  //       * n: Refractive index, [-]
  //       * Prandt: Prandtl number, [-]
  //       * Pr: Reduced Pressure, [-]
  //       * Tr: Reduced Temperature, [-]
  //       * Hvap: Vaporization heat, [kJ/kg]
  //       * Svap: Vaporization entropy, [kJ/kg·K]
  //
  //   Examples
  //   --------
  //   >>> water=IAPWS97(T=170+273.15, x=0.5)
  //   >>> water.Liquid.cp, water.Vapor.cp, water.Liquid.w, water.Vapor.w
  //   4.3695 2.5985 1418.3 498.78
  //
  //   >>> water=IAPWS97(T=325+273.15, x=0.5)
  //   >>> water.P, water.Liquid.v, water.Vapor.v, water.Liquid.h, water.Vapor.h
  //   12.0505 0.00152830 0.0141887 1493.37 2684.48
  //
  //   >>> water=IAPWS97(T=50+273.15, P=0.0006112127)
  //   >>> water.cp0, water.cv0, water.h0, water.s0, water.w0
  //   1.8714 1.4098 2594.66 9.471 444.93
  //

  public class IAPWS97 : object
  {
    public string _thermo;
    public int a;
    public double a0;
    public string CAS;
    public double cp0;
    public double cp0_cv;
    public double cv0;
    public object dipole;
    public object f_accent;
    public object g;
    public int g0;
    public double gamma0;
    public object h;
    public double h0;
    public object Hvap;
    public object Liquid;
    public double M;
    public string name;
    public object P;
    public object Pc;
    public object phase;
    public object Pr;
    public object region;
    public int rho;
    public object rhoc;
    public object s;
    public double s0;
    public object sigma;
    public object Svap;
    public string synonim;
    public object T;
    public object Tb;
    public object Tc;
    public object Tr;
    public object Tt;
    public int u;
    public double u0;
    public object v;
    public double v0;
    public object Vapor;
    public double w0;
    public object x;

    public Dictionary<string, double> kwargs = new Dictionary<object, object>
    {
      { "T", 0.0},
      { "P", 0.0},
      { "x", null},
      { "h", null},
      { "s", null},
      { "v", 0.0},
      { "l", 0.5893}
    };

    public int status = 0;
    public string msg = "Unknown variables";

    public IAPWS97(Hashtable kwargs)
    {
      this.kwargs = IAPWS97.kwargs.copy();
      this.@__call__(kwargs);
    }

    public virtual object @__call__(Hashtable kwargs)
    {
      this.kwargs.update(kwargs);
      if(this.calculable)
      {
        this.status = 1;
        this.calculo();
        this.msg = "Solved";
      }
    }

    // Check if class is calculable by its kwargs
    public object calculable
    {
      get
      {
        this._thermo = "";
        if(this.kwargs["T"] && this.kwargs["P"]) { this._thermo = "TP"; }
        else if(this.kwargs["P"] && this.kwargs["h"] != null) { this._thermo = "Ph"; }
        else if(this.kwargs["P"] && this.kwargs["s"] != null) { this._thermo = "Ps"; }
        else if(this.kwargs["h"] != null && this.kwargs["s"] != null)
        {
          // TODO: Add other pairs definitions options
          // elif self.kwargs["P"] and self.kwargs["v"]:
          // self._thermo = "Pv"
          // elif self.kwargs["T"] and self.kwargs["s"] is not None:
          // self._thermo = "Ts"
          this._thermo = "hs";
        }
        else if(this.kwargs["T"] && this.kwargs["x"] != null) { this._thermo = "Tx"; }
        else if(this.kwargs["P"] && this.kwargs["x"] != null) { this._thermo = "Px"; }
        return this._thermo;
      }
    }

    public virtual object calculo()
    {
      object Po;
      object s;
      object x;
      object To;
      object h;
      object vo;
      object rho;
      object region;
      object P;
      object T;
      Dictionary<string, double> propiedades = null;

      var args = (this.kwargs[this._thermo[0]], this.kwargs[this._thermo[1]]);
      if(this._thermo == "TP")
      {
        var _tup_1 = args;
        T = _tup_1.Item1;
        P = _tup_1.Item2;
        region = _Bound_TP(T, P);
        if(region == 1) { propiedades = _Region1(T, P); }
        else if(region == 2) { propiedades = _Region2(T, P); }
        else if(region == 3)
        {
          if(T == Tc && P == Pc) { rho = rhoc; }
          else
          {
            vo = _Backward3_v_PT(P, T);
            rho = newton(funcion, 1 / vo);
          }
          propiedades = _Region3(rho, T);
        }
        else if(region == 5)
        {
          propiedades = _Region5(T, P);
        }
        else
        {
          throw new NotImplementedException("Incoming out of bound");
        }
      }
      else if(this._thermo == "Ph")
      {
        var _tup_2 = args;
        P = _tup_2.Item1;
        h = _tup_2.Item2;
        region = _Bound_Ph(P, h);
        if(region == 1)
        {
          To = _Backward1_T_Ph(P, h);
          T = newton(T => _Region1(T, P)["h"] - h, To);
          propiedades = _Region1(T, P);
        }
        else if(region == 2)
        {
          To = _Backward2_T_Ph(P, h);
          T = newton(T => _Region2(T, P)["h"] - h, To);
          propiedades = _Region2(T, P);
        }
        else if(region == 3)
        {
          vo = _Backward3_v_Ph(P, h);
          To = _Backward3_T_Ph(P, h);
          var _tup_3 = fsolve(funcion, new List<int>
          {
            1 / vo,
            To
          });
          rho = _tup_3.Item1;
          T = _tup_3.Item2;
          propiedades = _Region3(rho, T);
        }
        else if(region == 4)
        {
          T = _TSat_P(P);
          if(T <= 623.15)
          {
            var h1 = _Region1(T, P)["h"];
            var h2 = _Region2(T, P)["h"];
            x = (h - h1) / (h2 - h1);
            propiedades = _Region4(P, x);
          }
          else
          {
            vo = _Backward3_v_Ph(P, h);
            To = _Backward3_T_Ph(P, h);
            var _tup_4 = fsolve(funcion, new List<int>
            {
              1 / vo,
              To
            });
            rho = _tup_4.Item1;
            T = _tup_4.Item2;
            propiedades = _Region3(rho, T);
          }
        }
        else if(region == 5)
        {
          T = newton(T => _Region5(T, P)["h"] - h, 1500);
          propiedades = _Region5(T, P);
        }
        else
        {
          throw new NotImplementedException("Incoming out of bound");
        }
      }
      else if(this._thermo == "Ps")
      {
        var _tup_5 = args;
        P = _tup_5.Item1;
        s = _tup_5.Item2;
        region = _Bound_Ps(P, s);
        if(region == 1)
        {
          To = _Backward1_T_Ps(P, s);
          T = newton(T => _Region1(T, P)["s"] - s, To);
          propiedades = _Region1(T, P);
        }
        else if(region == 2)
        {
          To = _Backward2_T_Ps(P, s);
          T = newton(T => _Region2(T, P)["s"] - s, To);
          propiedades = _Region2(T, P);
        }
        else if(region == 3)
        {
          vo = _Backward3_v_Ps(P, s);
          To = _Backward3_T_Ps(P, s);
          var _tup_6 = fsolve(funcion, new List<int>
          {
            1 / vo,
            To
          });
          rho = _tup_6.Item1;
          T = _tup_6.Item2;
          propiedades = _Region3(rho, T);
        }
        else if(region == 4)
        {
          T = _TSat_P(P);
          if(T <= 623.15)
          {
            var s1 = _Region1(T, P)["s"];
            var s2 = _Region2(T, P)["s"];
            x = (s - s1) / (s2 - s1);
            propiedades = _Region4(P, x);
          }
          else
          {
            vo = _Backward3_v_Ps(P, s);
            To = _Backward3_T_Ps(P, s);
            var _tup_7 = fsolve(funcion, new List<int>
            {
              1 / vo,
              To
            });
            rho = _tup_7.Item1;
            T = _tup_7.Item2;
            propiedades = _Region3(rho, T);
          }
        }
        else if(region == 5)
        {
          T = newton(T => _Region5(T, P)["s"] - s, 1500);
          propiedades = _Region5(T, P);
        }
        else
        {
          throw new NotImplementedException("Incoming out of bound");
        }
      }
      else if(this._thermo == "hs")
      {
        var _tup_8 = args;
        h = _tup_8.Item1;
        s = _tup_8.Item2;
        region = _Bound_hs(h, s);
        if(region == 1)
        {
          Po = _Backward1_P_hs(h, s);
          To = _Backward1_T_Ph(Po, h);
          var _tup_9 = fsolve(funcion, new List<object>
          {
            To,
            Po
          });
          T = _tup_9.Item1;
          P = _tup_9.Item2;
          propiedades = _Region1(T, P);
        }
        else if(region == 2)
        {
          Po = _Backward2_P_hs(h, s);
          To = _Backward2_T_Ph(Po, h);
          var _tup_10 = fsolve(funcion, new List<object>
          {
            To,
            Po
          });
          T = _tup_10.Item1;
          P = _tup_10.Item2;
          propiedades = _Region2(T, P);
        }
        else if(region == 3)
        {
          P = _Backward3_P_hs(h, s);
          vo = _Backward3_v_Ph(P, h);
          To = _Backward3_T_Ph(P, h);
          var _tup_11 = fsolve(funcion, new List<int>
          {
            1 / vo,
            To
          });
          rho = _tup_11.Item1;
          T = _tup_11.Item2;
          propiedades = _Region3(rho, T);
        }
        else if(region == 4)
        {
          if(round(s - sc, 6) == 0 && round(h - hc, 6) == 0)
          {
            propiedades = _Region3(rhoc, Tc);
          }
          else
          {
            To = _Backward4_T_hs(h, s);
            if(To < 273.15 || To > Tc)
            {
              To = 300;
            }
            var _tup_12 = fsolve(funcion, new List<double>
            {
              To,
              0.5
            });
            T = _tup_12.Item1;
            x = _tup_12.Item2;
            P = _PSat_T(T);
            if(Pt <= P < Pc && 0 < x < 1)
            {
              propiedades = _Region4(P, x);
            }
            else if(Pt <= P <= Ps_623 && x == 0)
            {
              propiedades = _Region1(T, P);
            }
          }
        }
        else if(region == 5)
        {
          var _tup_13 = fsolve(funcion, new List<int>
          {
            1400,
            1
          });
          T = _tup_13.Item1;
          P = _tup_13.Item2;
          propiedades = _Region5(T, P);
        }
        else
        {
          throw new NotImplementedException("Incoming out of bound");
        }
      }
      else if(this._thermo == "Px")
      {
        var _tup_14 = args;
        P = _tup_14.Item1;
        x = _tup_14.Item2;
        T = _TSat_P(P);
        if(Pt <= P < Pc && 0 < x < 1)
        {
          propiedades = _Region4(P, x);
        }
        else if(Pt <= P <= Ps_623 && x == 0)
        {
          propiedades = _Region1(T, P);
        }
        else if(Pt <= P <= Ps_623 && x == 1)
        {
          propiedades = _Region2(T, P);
        }
        else if(Ps_623 < P < Pc && (0, 1).Contains(x))
        {
          var rhoo = 1.0 / _Backward3_sat_v_P(P, T, x);
          rho = fsolve(funcion, rhoo)[0];
          propiedades = _Region3(rho, T);
        }
        else if(P == Pc && 0 <= x <= 1)
        {
          propiedades = _Region3(rhoc, Tc);
        }
        else
        {
          throw new NotImplementedException("Incoming out of bound");
        }
        this.sigma = _iapws._Tension(T);
        propiedades["x"] = x;
      }
      else if(this._thermo == "Tx")
      {
        var _tup_15 = args;
        T = _tup_15.Item1;
        x = _tup_15.Item2;
        P = _PSat_T(T);
        if(273.15 <= T < Tc && 0 < x < 1) { propiedades = _Region4(P, x); }
        else if(273.15 <= T <= 623.15 && x == 0) { propiedades = _Region1(T, P); }
        else if(273.15 <= T <= 623.15 && x == 1) { propiedades = _Region2(T, P); }
        else if(623.15 < T < Tc && (0, 1).Contains(x))
        {
          rho = 1.0 / _Backward3_sat_v_P(P, T, x);
          propiedades = _Region3(rho, T);
        }
        else if(T == Tc && 0 <= x <= 1) { propiedades = _Region3(rhoc, Tc); }
        else
        {
          throw new NotImplementedException("Incoming out of bound");
        }
        this.sigma = _iapws._Tension(T);
        propiedades["x"] = x;
      }
      Func<object, object> funcion = rho =>
      {
        return _Region3(rho, this.kwargs["T"])["P"] - P;
      };
      Func<object, object> funcion = par => { return (_Region3(par[0], par[1])["h"] - h, _Region3(par[0], par[1])["P"] - P); };
      Func<object, object> funcion = par => { return (_Region3(par[0], par[1])["h"] - h, _Region3(par[0], par[1])["P"] - P); };
      Func<object, object> funcion = par => { return (_Region3(par[0], par[1])["s"] - s, _Region3(par[0], par[1])["P"] - P); };
      Func<object, object> funcion = par => { return (_Region3(par[0], par[1])["s"] - s, _Region3(par[0], par[1])["P"] - P); };
      Func<object, object> funcion = par => { return (_Region1(par[0], par[1])["h"] - h, _Region1(par[0], par[1])["s"] - s); };
      Func<object, object> funcion = par => { return (_Region2(par[0], par[1])["h"] - h, _Region2(par[0], par[1])["s"] - s); };
      Func<object, object> funcion = par => { return (_Region3(par[0], par[1])["h"] - h, _Region3(par[0], par[1])["s"] - s); };
      Func<object, object> funcion = par =>
      {
        if(par[1] < 0) { par[1] = 0; }
        else if(par[1] > 1) { par[1] = 1; }
        if(par[0] < 273.15) { par[0] = 273.15; }
        else if(par[0] > Tc) { par[0] = Tc; }
        var Po = _PSat_T(par[0]);
        var liquid = _Region1(par[0], Po);
        var vapor = _Region2(par[0], Po);
        var hl = liquid["h"];
        var sl = liquid["s"];
        var hv = vapor["h"];
        var sv = vapor["s"];
        return (hv * par[1] + hl * (1 - par[1]) - h, sv * par[1] + sl * (1 - par[1]) - s);
      };
      Func<object, object> funcion = par => { return (_Region5(par[0], par[1])["h"] - h, _Region5(par[0], par[1])["s"] - s); };
      Func<object, object> funcion = rho => { return _Region3(rho, T)["P"] - P; };
      this.M = 18.015257;
      this.Pc = Pc;
      this.Tc = Tc;
      this.rhoc = rhoc;
      this.Tt = Tt;
      this.Tb = Tb;
      this.f_accent = _iapws.f_acent;
      this.dipole = _iapws.Dipole;
      this.x = propiedades["x"];
      this.region = propiedades["region"];
      this.name = "water";
      this.synonim = "R-718";
      this.CAS = "7732-18-5";
      this.T = propiedades["T"];
      this.P = propiedades["P"];
      this.v = propiedades["v"];
      this.rho = 1 / this.v;
      this.phase = _utils.getphase(this.Tc, this.Pc, this.T, this.P, this.x, this.region);
      this.Tr = this.T / this.Tc;
      this.Pr = this.P / this.Pc;
      // Ideal properties
      if(new List<int>
      {
        2,
        5
      }.Contains(this.region))
      {
        var cp0 = prop0(this.T, this.P);
        this.v0 = cp0["v"];
        this.h0 = cp0["h"];
        this.u0 = this.h0 - this.P * 1000 * this.v0;
        this.s0 = cp0["s"];
        this.a0 = this.u0 - this.T * this.s0;
        this.g0 = this.h0 - this.T * this.s0;
        this.cp0 = cp0["cp"];
        this.cv0 = cp0["cv"];
        this.cp0_cv = this.cp0 / this.cv0;
        this.w0 = cp0["w"];
        this.gamma0 = this.cp0_cv;
      }
      else
      {
        this.v0 = null;
        this.h0 = null;
        this.u0 = null;
        this.s0 = null;
        this.a0 = null;
        this.g0 = 0;
        this.cp0 = null;
        this.cv0 = null;
        this.cp0_cv = null;
        this.w0 = null;
        this.gamma0 = null;
      }
      this.Liquid = _utils._fase();
      this.Vapor = _utils._fase();
      if(this.x == 0)
      {
        // only liquid phase
        this.fill(this, propiedades);
        this.fill(this.Liquid, propiedades);
        this.sigma = _iapws._Tension(this.T);
      }
      else if(this.x == 1)
      {
        // only vapor phase
        this.fill(this, propiedades);
        this.fill(this.Vapor, propiedades);
      }
      else
      {
        // two phases
        var liquido = _Region1(this.T, this.P);
        this.fill(this.Liquid, liquido);
        var vapor = _Region2(this.T, this.P);
        this.fill(this.Vapor, vapor);
        this.h = propiedades["h"];
        this.u = this.h - this.P * 1000 * this.v;
        this.s = propiedades["s"];
        this.a = this.u - this.T * this.s;
        this.g = this.h - this.T * this.s;
        this.sigma = _iapws._Tension(this.T);
        this.Hvap = vapor["h"] - liquido["h"];
        this.Svap = vapor["s"] - liquido["s"];
      }
    }

    public virtual void fill(_utils._fase fase, Dictionary<string, double> estado)
    {
      fase.v = estado["v"];
      fase.rho = 1 / fase.v;
      fase.h = estado["h"];
      fase.s = estado["s"];
      fase.u = fase.h - this.P * 1000 * fase.v;
      fase.a = fase.u - this.T * fase.s;
      fase.g = fase.h - this.T * fase.s;
      fase.cv = estado["cv"];
      fase.cp = estado["cp"];
      fase.cp_cv = fase.cp / fase.cv;
      fase.w = estado["w"];
      fase.Z = this.P * fase.v / R * 1000 / this.T;
      fase.alfav = estado["alfav"];
      fase.xkappa = estado["kt"];
      fase.kappas = -1 / fase.v * this.derivative("v", "P", "s", fase);
      fase.joule = this.derivative("T", "P", "h", fase);
      fase.deltat = this.derivative("h", "P", "T", fase);
      fase.gamma = -fase.v / this.P * this.derivative("P", "v", "s", fase);
      fase.alfap = fase.alfav / this.P / fase.xkappa;
      fase.betap = -1 / this.P * this.derivative("P", "v", "T", fase);
      fase.fi = Math.Exp((fase.g - this.g0) / R / this.T);
      fase.f = this.P * fase.fi;
      fase.mu = _iapws._Viscosity(fase.rho, this.T);
      // Use industrial formulation for critical enhancement in thermal
      // conductivity calculation
      fase.drhodP_T = this.derivative("rho", "P", "T", fase);
      fase.k = _iapws._ThCond(fase.rho, this.T, fase);
      fase.nu = fase.mu / fase.rho;
      fase.alfa = fase.k / 1000 / fase.rho / fase.cp;

      try { fase.epsilon = _iapws._Dielectric(fase.rho, this.T); }
      catch (NotImplementedError) { fase.epsilon = null; }

      fase.Prandt = fase.mu * fase.cp * 1000 / fase.k;
      try { fase.n = _iapws._Refractive(fase.rho, this.T, this.kwargs["l"]); }
      catch (NotImplementedError) { fase.n = null; }
    }

    // Wrapper derivative for custom derived properties
    //     where x, y, z can be: P, T, v, u, h, s, g, a
    public virtual object derivative(object z, object x, object y, object fase)
    {
      return _utils.deriv_G(this, z, x, y, fase);
    }

  }

  public class IAPWS97_PT : IAPWS97 { public IAPWS97_PT(object P, object T) : base(T: T, P: P) { } } // Derivated class for direct P and T input
  public class IAPWS97_Ph : IAPWS97 { public IAPWS97_Ph(object P, object h) : base(P: P, h: h) { } } // Derivated class for direct P and h input
  public class IAPWS97_Ps : IAPWS97 { public IAPWS97_Ps(object P, object s) : base(P: P, s: s) { } } // Derivated class for direct P and s input
  public class IAPWS97_Px : IAPWS97 { public IAPWS97_Px(object P, object x) : base(P: P, x: x) { } } // Derivated class for direct P and x input
  public class IAPWS97_Tx : IAPWS97 { public IAPWS97_Tx(object T, object x) : base(T: T, x: x) { } } // Derivated class for direct T and x input

}

