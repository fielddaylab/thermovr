/*
DOCUMENTATION- phil, 12/16/19
The intent of this class is to (statelessly) abstract implementation details across potentially multiple different "libraries".
The general idea is that IAPWS is a specification which has been implemented in part many different times, with many different APIs.
Some of the APIs offer a subset of functionality only allowing one to query x given y and z, or y given x and z (etc...), and often do so with differing expected units.
The purpose of this class is make compile the total implementation and make universal in form and unit all queries regarding equations of state for water.

This class should be "pure math". The only state should be for caching/performance reasons.

Current implementations ported/used:
IF97 (IAPWS97.cs) https://github.com/CoolProp/IF97
IAPWS95 (IAPWS95.cs) https://code.google.com/archive/p/proph2o/downloads
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
  internalenergy = i (<- thermodynamics calls this "Q", but that's what we use for "quality". so I'm using "i")
  entropy = s
  enthalpy = h
  quality = q
  */

  //Pa
  public static double p_min; // 611.213
  public static double p_max; // 100000000
  //Pa
  public static double psat_min; //TODO: add actual value as comment for quick reference
  public static double psat_max; //TODO: add actual value as comment for quick reference
  //M^3/kg
  public static double v_min;  // 0.0003
  public static double v_max; // 1000
  //K
  public static double t_min; // 273.15
  public static double t_max; // 1073.15
  //J/kg
  public static double i_min; //0 //TODO: add actual value as comment for quick reference
  public static double i_max; //0 //TODO: add actual value as comment for quick reference
  //J/kgK
  public static double s_min; //0 //TODO: add actual value as comment for quick reference
  public static double s_max; //0 //TODO: add actual value as comment for quick reference
  //J/kg
  public static double h_min; //0 //TODO: add actual value as comment for quick reference
  public static double h_max; //0 //TODO: add actual value as comment for quick reference
  //%
  public static double q_min; //0
  public static double q_max; //1

  public static void Init()
  {
    IF97.initRegions();

    //Pa
    p_min = IF97.get_Pmin()*1000000.0; // 611.213
    p_max = IF97.get_Pmax()*1000000.0; // 100000000
    //Pa
    psat_min = IF97.get_ptrip()*1000000.0; //TODO: add actual value as comment for quick reference
    psat_max = IF97.get_pcrit()*1000000.0; //TODO: add actual value as comment for quick reference
    //M^3/kg
    v_min = 1.0/3000;  // 0.0003
    v_max = 1.0/0.001; // 1000
    //K
    t_min = IF97.get_Tmin(); // 273.15
    t_max = IF97.get_Tmax(); // 1073.15
    //J/kg
    i_min = 0; //TODO:find actual min
    i_max = 1; //TODO:find actual max
    //J/kgK
    s_min = IF97.Smin; //TODO: add actual value as comment for quick reference
    s_max = IF97.Smax; //TODO: add actual value as comment for quick reference
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

  //helpers to easily generate values "randomly", ensuring we're within the range we care about- NOT PHYSICALLY BASED
  //"percent" = "percent between min value and max value across given dimension"
  public static double p_given_percent(   double t) { return Lerpd(p_min,p_max,t); }
  public static double psat_given_percent(double t) { return Lerpd(psat_min,psat_max,t); }
  public static double v_given_percent(   double t) { return Lerpd(v_min,v_max,t); }
  public static double t_given_percent(   double t) { return Lerpd(t_min,t_max,t); }
  public static double i_given_percent(   double t) { return Lerpd(i_min,i_max,t); }
  public static double s_given_percent(   double t) { return Lerpd(s_min,s_max,t); }
  public static double h_given_percent(   double t) { return Lerpd(h_min,h_max,t); }
  public static double q_given_percent(   double t) { return t; } //q already is a percent

  public static double percent_given_p(   double p) { return p-p_min/(p_max-p_min); }
  public static double percent_given_psat(double psat) { return psat-psat_min/(psat_max-psat_min); }
  public static double percent_given_v(   double v) { return v-v_min/(v_max-v_min); }
  public static double percent_given_t(   double t) { return t-t_min/(t_max-t_min); }
  public static double percent_given_i(   double i) { return i-i_min/(i_max-i_min); }
  public static double percent_given_s(   double s) { return s-s_min/(s_max-s_min); }
  public static double percent_given_h(   double h) { return h-h_min/(h_max-h_min); }
  public static double percent_given_q(   double q) { return q; } //q already is a percent

  //rule of naming for consistency: prefer lexical ordering "p < v < t", ie "p_given_vt" rather than "p_given_tv"

  public static double p_given_vt(double v, double t)
  {
    return IAPWS95.IAPWS95_pressure(1.0/v,t)*1000.0; //expects:Kg/M^3,K returns KPa
  }

  public static double v_given_pt(double p, double t)
  {
    return 1.0/IF97.rhomass_Tp(t,p/1000000.0); //expects:K,MPa returns Kg/M^3
  }

  public static double t_given_pv(double p, double v)
  {
    return t_given_percent(0.5); //bogus implementation!
  }

  public static double tsat_given_p(double p)
  {
    return IF97.Tsat97(p/1000000.0);
  }

  public static double vliq_given_p(double p)
  {
    return 1.0/IF97.rholiq_p(p/1000000.0); //expects:MPa returns Kg/M^3
  }

  public static double vvap_given_p(double p)
  {
    return 1.0/IF97.rhovap_p(p/1000000.0); //expects:MPa returns Kg/M^3
  }

  public static double i_given_pt(double p, double t)
  {
    return i_given_percent(0.5); //TODO:
  }

  public static double s_given_pt(double p, double t)
  {
    return s_given_percent(0.5); //TODO:
  }

  public static double h_given_pt(double p, double t)
  {
    return h_given_percent(0.5); //TODO:
  }

  public static double q_given_pt(double p, double t)
  {
    return q_given_percent(0.5); //TODO:
  }


}

