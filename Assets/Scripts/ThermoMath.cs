﻿/*
DOCUMENTATION- phil, 12/16/19
The intent of this class is to (statelessly) abstract implementation details across potentially multiple different "libraries".
The general idea is that IAPWS is a specification which has been implemented in part many different times, with many different APIs.
Some of the APIs offer a subset of functionality only allowing one to query x given y and z, or y given x and z (etc...), and often do so with differing expected units.
The purpose of this class is make compile the total implementation and make universal in form and unit all queries regarding equations of state for water.

This class should be "pure math". The only state should be for caching/performance reasons.

Current implementations ported/used:
IF97 (IAPWS97.cs) https://github.com/CoolProp/IF97
IAPWS95 (IAPWS95.cs) https://code.google.com/archive/p/proph2o/downloads


Note- phil, 12/27/19
I have only implemented functions that were determined to be NEEDED to implement ThermoState.cs
Many other derivations are possible, either directly, or via iterative guessing (see "iterate_x_given_x_verify_x" examples)
*/

using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public static class ThermoMath
{
  /*
  //IF97 API
  public static double rhomass_Tp(double T, double p)     // Get the mass density [kg/m^3] as a function of T [K] and p [Pa]
  public static double hmass_Tp(double T, double p)       // Get the mass enthalpy [J/kg] as a function of T [K] and p [Pa]
  public static double smass_Tp(double T, double p)       // Get the mass entropy [J/kg/K] as a function of T [K] and p [Pa]
  public static double umass_Tp(double T, double p)       // Get the mass internal energy [J/kg] as a function of T [K] and p [Pa]
  public static double cpmass_Tp(double T, double p)      // Get the mass constant-pressure specific heat [J/kg/K] as a function of T [K] and p [Pa]
  public static double cvmass_Tp(double T, double p)      // Get the mass constant-volume specific heat [J/kg/K] as a function of T [K] and p [Pa]
  public static double speed_sound_Tp(double T, double p) // Get the speed of sound [m/s] as a function of T [K] and p [Pa]
  public static double drhodp_Tp(double T, double p)      // Get the [d(rho)/d(p)]T [kg/mï¿½/Pa] as a function of T [K] and p [Pa]
  //IF97 Paper verified units:
  rhomass_Tp(T,p) | 1.0/rhomass_Tp(300, 3) = 0.00100215 | p:MPa, v:M^3/Kg, T:K | expects:K,MPa returns Kg/M^3
  */

  /*
  //IAPWS95 API
  public static double IAPWS95_pressure(double rho, double T);                                   //Input: rho in kg/m3, T in K, Output: Pa
  public static double IAPWS95_internal_energy(double rho, double T);                            //Input: rho in kg/m3, T in K, Output: Pa
  public static double IAPWS95_entropy(double rho, double T);                                    //Input: rho in kg/m3, T in K, Output: kJ/kg-K
  public static double IAPWS95_enthalpy(double rho, double T);                                   //Input: rho in kg/m3, T in K, Output: kJ/kg
  public static double IAPWS95_isochoric_heat_capacity(double rho, double T);                    //Input: rho in kg/m3, T in K, Output: kJ/kg-K
  public static double IAPWS95_isobaric_heat_capacity(double rho, double T);                     //Input: rho in kg/m3, T in K, Output: kJ/kg-K
  public static double IAPWS95_speed_of_sound(double rho, double T);                             //Input: rho in kg/m3, T in K, Output: m/s
  public static double IAPWS95_joule_thompson_coefficient(double rho, double T);                 //Input: rho in kg/m3, T in K
  public static double IAPWS95_isothermal_throttling_coefficient(double rho, double T);          //Input: rho in kg/m3, T in K
  public static double IAPWS95_isentropic_temperature_pressure_coefficent(double rho, double T); //Input: rho in kg/m3, T in K
  //IAPWS95 Paper verified units:
  IAPWS95_pressure(rho, T) | IAPWS95_pressure(999.887406, 275.0)/1000 = 0.0006982125 | p:MPa, v:Kg/M^3, t:K | expects:Kg/M^3,K returns KPa
  */

  /*
  pressure = p
  specificvolume = v
  temperature = t
  internalenergy = u
  entropy = s
  enthalpy = h
  quality = x
  */

  //Pa
  public static double p_min; // 611.213
  public static double p_max; // 100000000
  //Pa
  public static double psat_min; //TODO: comment actual value for quick reference
  public static double psat_max; //TODO: comment actual value for quick reference
  //M^3/kg
  public static double v_min;  // 0.0003
  public static double v_max; // 1000
  //K
  public static double t_min; // 273.15
  public static double t_max; // 1073.15
  //J/kg
  public static double u_min; //TODO: comment actual value for quick reference
  public static double u_max; //TODO: comment actual value for quick reference
  //J/kgK
  public static double s_min; //TODO: comment actual value for quick reference
  public static double s_max; //TODO: comment actual value for quick reference
  //J/kg
  public static double h_min; //TODO: comment actual value for quick reference
  public static double h_max; //TODO: comment actual value for quick reference
  //%
  public static double x_min; //0
  public static double x_max; //1

  public static void Init()
  {
    IF97.initRegions();

    //Pa
    p_min = IF97.get_Pmin()*1000000.0; // 611.213
    p_max = IF97.get_Pmax()*1000000.0; // 100000000
    //Pa
    psat_min = IF97.get_ptrip()*1000000.0; //TODO: comment actual value for quick reference
    psat_max = IF97.get_pcrit()*1000000.0; //TODO: comment actual value for quick reference
    //M^3/kg
    v_min = 1.0/3000;  // 0.0003
    v_max = 1.0/0.001; // 1000
    //K
    t_min = IF97.get_Tmin(); // 273.15
    t_max = IF97.get_Tmax(); // 1073.15
    //J/kg
    u_min = 0; //TODO:find actual min
    u_max = 9999999999; //TODO:find actual max
    //J/kgK
    s_min = IF97.Smin; //TODO: comment actual value for quick reference
    s_max = IF97.Smax; //TODO: comment actual value for quick reference
    //J/kg
    h_min = IF97.Hmin(s_min); //TODO: unsure if this is correct calculation!
    h_max = IF97.Hmax(s_max); //TODO: unsure if this is correct calculation!

    //DEBUG INFO:
    //IF97.print_tables();
    //IAPWS95.print_tables();
    //compare_impls();
  }

  static double sample_lbase = 1.6f;
  static double sample(double t) { return Math.Pow(t,sample_lbase); }
  public static void compare_impls()
  {
    int samples = 350; //match value in ThermoState to match sample points

    for(int y = 0; y < samples; y++)
    {
      double pt = ((double)y/(samples-1));
      for(int z = 0; z < samples; z++)
      {
        double tt = ((double)z/(samples-1));
        double pst = sample(pt);
        double tst = sample(tt);
        double p = Lerpd(p_min,p_max,pst);
        double t = Lerpd(t_min,t_max,tst);
        double v = v_given_pt(p,t);
        //pvt in Pa, M^3/Kg, K
        double _p = p_given_vt(v,t);

        Debug.LogFormat("error:{0} p:{1}Pa ({2}Pa), v:{3}M^3/Kg, t:{4}K",p-_p,p,_p,v,t);
      }
    }

  }

  static double Lerpd(double a, double b, double t) { return (b-a)*t+a; }

  //region: 0 subcooled liquid, 1 two-phase, 2 superheated vapor
  public static int region_given_pvt(double p, double v, double t)
  {
    IF97.IF97REGIONS r = IF97.RegionDetermination_TP(t, p/1000000.0);
    switch(r)
    {
      case IF97.IF97REGIONS.REGION_1: //liquid
        Debug.Log("1"); //subcooled liquid
        break;
      case IF97.IF97REGIONS.REGION_2: //vapor
        Debug.Log("2"); //superheated vapor
        break;
      case IF97.IF97REGIONS.REGION_3:
        Debug.Log("3"); //superheated vapor
        break;
      case IF97.IF97REGIONS.REGION_4: //two-phase
        Debug.Log("4"); //superheated vapor
        break;
      case IF97.IF97REGIONS.REGION_5:
        Debug.Log("5"); //superheated vapor
        break;
    }

    switch(r)
    {
      case IF97.IF97REGIONS.REGION_1: //liquid
        return 0; //subcooled liquid
      case IF97.IF97REGIONS.REGION_2: //vapor
        return 2; //superheated vapor
      case IF97.IF97REGIONS.REGION_3:
      case IF97.IF97REGIONS.REGION_4: //two-phase
      case IF97.IF97REGIONS.REGION_5:
        return 1; //two-phase
    }
    return -1;
  }

  //helpers to easily generate values "randomly", ensuring we're within the range we care about- NOT PHYSICALLY BASED
  //"percent" = "percent between min value and max value across given dimension"
  public static double p_given_percent(   double t) { return Lerpd(p_min,p_max,t); }
  public static double psat_given_percent(double t) { return Lerpd(psat_min,psat_max,t); }
  public static double v_given_percent(   double t) { return Lerpd(v_min,v_max,t); }
  public static double t_given_percent(   double t) { return Lerpd(t_min,t_max,t); }
  public static double u_given_percent(   double t) { return Lerpd(u_min,u_max,t); }
  public static double s_given_percent(   double t) { return Lerpd(s_min,s_max,t); }
  public static double h_given_percent(   double t) { return Lerpd(h_min,h_max,t); }
  public static double x_given_percent(   double t) { return t; } //x already is a percent

  public static double percent_given_p(   double p) { return (p-p_min)/(p_max-p_min); }
  public static double percent_given_psat(double psat) { return (psat-psat_min)/(psat_max-psat_min); }
  public static double percent_given_v(   double v) { return (v-v_min)/(v_max-v_min); }
  public static double percent_given_t(   double t) { return (t-t_min)/(t_max-t_min); }
  public static double percent_given_u(   double u) { return (u-u_min)/(u_max-u_min); }
  public static double percent_given_s(   double s) { return (s-s_min)/(s_max-s_min); }
  public static double percent_given_h(   double h) { return (h-h_min)/(h_max-h_min); }
  public static double percent_given_x(   double x) { return x; } //x already is a percent

  //rule of naming for consistency: prefer lexical ordering "p < v < t < u < s < h < q", ie "p_given_vt" rather than "p_given_tv"
  //note- where it says "UNIT CONVERSION UNTESTED", that measn I haven't yet verified that it performs the correct translation of units between APIs (ie, one may be expecting KPa, and another just Pa). I have learned to not trust the documentation of either API (IAPWS95 or IF97)

  public static double p_given_vt(double v, double t)
  {
    return IAPWS95.IAPWS95_pressure(1.0/v,t)*1000.0; //expects:Kg/M^3,K returns KPa
  }

  public static double v_given_pt(double p, double t) //DO NOT USE IN VAPOR DOME
  {
    return 1.0/IF97.rhomass_Tp(t,p/1000000.0); //expects:K,MPa returns Kg/M^3
  }

  public static double v_given_ph(double p, double h)
  {
    return 1.0/IF97.rhomass_phmass(p/1000000.0,h/1000.0); //UNIT CONVERSION UNTESTED!
  }

  public static double t_given_ph(double p, double h)
  {
    return IF97.T_phmass(p/1000000.0,h/1000.0); //UNIT CONVERSION UNTESTED!
  }

  public static double tsat_given_p(double p)
  {
    return IF97.Tsat97(p/1000000.0); //UNIT CONVERSION UNTESTED!
  }

  public static double vliq_given_p(double p)
  {
    return 1.0/IF97.rholiq_p(p/1000000.0); //expects:MPa returns Kg/M^3
  }

  public static double vvap_given_p(double p)
  {
    return 1.0/IF97.rhovap_p(p/1000000.0); //expects:MPa returns Kg/M^3
  }

  public static double u_given_pt(double p, double t) //DO NOT USE IN VAPOR DOME
  {
    return IF97.umass_Tp(t, p/1000000.0); //UNIT CONVERSION UNTESTED!
  }

  public static double u_given_vt(double v, double t)
  {
    return IAPWS95.IAPWS95_internal_energy(1f/v,t); //UNIT CONVERSION UNTESTED!
  }

  public static double s_given_vt(double v, double t)
  {
    return IAPWS95.IAPWS95_entropy(1f/v,t)*1000f; //UNIT CONVERSION UNTESTED!
  }

  public static double h_given_vt(double v, double t)
  {
    return IAPWS95.IAPWS95_enthalpy(1f/v,t)*1000f; //UNIT CONVERSION UNTESTED!
  }

  public static double x_given_pv(double p, double v) //ONLY USE IN VAPOR DOME
  {
    //f means saturated liquid,
    //g means saturated gas
    double vf = IF97.cvliq_p(p/1000000.0);
    double vg = IF97.cvvap_p(p/1000000.0);
    return (v-vf)/(vg-vf);
  }

  public static double iterate_t_given_p_verify_u(double t, double p, double u) //t = "first guess"
  {
    int MAX_ITERS = 100;     //TODO: define intentionally
    double MAX_DELTA = 0.01; //TODO: define intentionally
    double step = 0.01;      //TODO: define intentionally
    double guess = t;
    double mark = u;
    double delta = Math.Abs(u_given_pt(p,guess)-mark);
    for(int i = 0; i < MAX_ITERS || delta < MAX_DELTA; i++)
    {
      double delta_a = Math.Abs(u_given_pt(p,guess+ step     )-mark);
      double delta_b = Math.Abs(u_given_pt(p,guess-(step/2.0))-mark);
      if(delta_a < delta_b)
      {
        delta = delta_a;
      }
      else
      {
        delta = delta_b;
        step = step/-2.0;
      }
      guess += step;
    }
    return guess;
  }

  public static double iterate_t_given_v_verify_u(double t, double v, double u)
  {
    int MAX_ITERS = 100;      //TODO: define intentionally
    double MAX_DELTA = 0.01;  //TODO: define intentionally
    double step = 0.01;       //TODO: define intentionally
    double guess = t;
    double mark = u;
    double delta = u_given_vt(v,guess)-mark;
    for(int i = 0; i < MAX_ITERS || delta < MAX_DELTA; i++)
    {
      double delta_a = Math.Abs(u_given_vt(v,guess+ step     )-mark);
      double delta_b = Math.Abs(u_given_vt(v,guess-(step/2.0))-mark);
      if(delta_a < delta_b)
      {
        delta = delta_a;
      }
      else
      {
        delta = delta_b;
        step = step/-2.0;
      }
      guess += step;
    }
    return guess;
  }
}

