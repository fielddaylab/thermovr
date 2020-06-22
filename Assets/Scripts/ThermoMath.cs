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


Note- phil, 12/27/19
I have only implemented functions that were determined to be NEEDED to implement ThermoState.cs
Many other derivations are possible, either directly, or via iterative guessing (see "iterate_x_given_x_verify_x" examples)
*/

using System;
using System.Collections;
using System.Collections.Generic;
using System.Threading.Tasks;
using UnityEngine;

public static class ThermoMath
{
  private const int MAX_MILLISECOND = 50;
  /*
  pressure = p
  specificvolume = v
  temperature = t
  internalenergy = u
  entropy = s
  enthalpy = h
  quality = x
  */

  /*
  max-min = the total range expected by the simulation
  neutral = "room temperature, 1 atm, 1kg water"
  smallstep = a "small step" in the given range (useful for cache invalidation threshhold)
  */

  //Pa
  public static double p_min;
  public static double p_max;
  public static double p_neutral;
  public static double p_smallstep;
  //Pa
  public static double psat_min;
  public static double psat_max;
  //M³/kg
  public static double v_min;
  public static double v_max;
  public static double v_neutral;
  public static double v_smallstep;
  //K
  public static double t_min;
  public static double t_max;
  public static double t_neutral;
  public static double t_smallstep;
  //J/kg
  public static double u_min;
  public static double u_max;
  public static double u_neutral;
  public static double u_smallstep;
  //J/kgK
  public static double s_min;
  public static double s_max;
  public static double s_neutral;
  public static double s_smallstep;
  //J/kg
  public static double h_min;
  public static double h_max;
  public static double h_neutral;
  public static double h_smallstep;
  //%
  public static double x_min;
  public static double x_max;
  public static double x_neutral;
  public static double x_smallstep;

  // A sort of notifier for users of the class.
  // Whenever an exception occurs, we set this to true.
  // User of the class can then check and react appropriately, setting false after handling the error.
  public static bool got_error;

  public static void Init()
  {
    got_error = false;
    IF97.initRegions();

    //Pa
    p_min = IF97.get_Pmin()*1000000.0; // 611.213
    p_max = IF97.get_Pmax()*1000000.0; // 100000000
    p_neutral = 101325.0;
    p_smallstep = 1.0;

    //Pa
    psat_min = IF97.get_ptrip()*1000000.0; //TODO: comment actual value for quick reference
    psat_max = IF97.get_pcrit()*1000000.0; //TODO: comment actual value for quick reference

    //M³/kg
    v_min = 1.0/3000;  // 0.0003
    v_max = 1.0/0.001; // 1000
    v_neutral = 0.001;
    v_smallstep = 0.00001;

    //K
    t_min = IF97.get_Tmin(); // 273.15
    t_max = IF97.get_Tmax(); // 1073.15
    t_neutral = 293.0;
    t_smallstep = 0.001;

    //J/kg
    u_min = 0.0; //TODO:find actual min
    u_max =  9999999999; //TODO:find actual max
    u_neutral = 83.28;
    u_smallstep = 0.0; //TODO: find

    //J/kgK
    //s_min = IF97.Smin; //TODO: comment actual value for quick reference //I don't think this is correct
    //s_max = IF97.Smax; //11.9210548250511 //I don't think this is correct...
    s_min = 0.0; //TODO: actually find something coherent
    s_max =  999999999999.0; //TODO: actually find something coherent
    s_neutral = 294.322;
    s_smallstep = 0.0; //TODO: find

    //J/kg
    //h_? = ; //experimentally derived- room temp water 
    //h_min = IF97.Hmin(s_min); //TODO: I don't think this is correct...
    //h_max = IF97.Hmax(s_max); //4171.65498424024 given s_max 11.9... //TODO: I don't think this is correct...
    h_min = 0.0; //TODO: actually come up with something coherent
    h_max =  9999999.0; //TODO: actually come up with something coherent
    h_neutral = 83377.0;
    h_smallstep = 0.0; //TODO: find

    //%
    x_min = 0.0;
    x_max = 1.0;
    x_neutral = 0;
    x_smallstep = 0.00001;

    //DEBUG INFO:
    //IF97.print_tables();
    //IAPWS95.print_tables();
    //compare_impls();
  }

  static double SAMPLE_LBASE = 1.6f;
  static double sample(double t) { return Math.Pow(t,SAMPLE_LBASE); }
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
        //pvt in Pa, M³/Kg, K
        double _p = p_given_vt(v,t);

        Debug.LogFormat("error:{0} p:{1}Pa ({2}Pa), v:{3}M³/Kg, t:{4}K",p-_p,p,_p,v,t);
      }
    }

  }

  static double Lerpd(double a, double b, double t) { return (b-a)*t+a; }

  //a bunch of options for getting region here- still need to figure out most reliable
  //region: 0 subcooled liquid, 1 two-phase, 2 superheated vapor
  public static int region_given_pvt(double p, double v, double t)
  {
    int liq = 0;
    int two = 1;
    int vapor = 2;

    //broad check w/ t
    if(t-t_smallstep > IF97.Tsat97(p/1000000.0)) return vapor;
    if(t+t_smallstep < IF97.Tsat97(p/1000000.0)) return liq;

    //broad check w/ p - unneeded
    //if(p-p_smallstep > IF97.psat97(t)) return liq;
    //if(p+p_smallstep < IF97.psat97(t)) return vapor;

    //fine check w/ v
    //f means saturated liquid,
    //g means saturated gas
    double vf = 1.0/IF97.rholiq_p(p/1000000.0);
    if(v <= vf) return liq;
    double vg = 1.0/IF97.rhovap_p(p/1000000.0);
    if(v >= vg) return vapor;
    return two;
  }
  /*
  public static int region_given_pvt(double p, double v, double t)
  {
    IF97.IF97REGIONS r = IF97.RegionDetermination_TP(t, p/1000000.0);
    switch(r)
    {
      case 1: //liquid
        return 0; //subcooled liquid
      case 2: //vapor
        return 2; //superheated vapor
      case 3:
      case 4: //two-phase
      case 5:
        return 1; //two-phase
    }
    return -1;
  }
  */
  public static int region_given_ph(double p, double h)
  {
    int r = IF97.Region_ph(p/1000000.0, h/1000.0); 
    switch(r)
    {
      case 1: //liquid
        return 0; //subcooled liquid
      case 2: //vapor
        return 2; //superheated vapor
      case 3:
      case 4: //two-phase
      case 5:
        return 1; //two-phase
    }
    return -1;
  }
  public static int region_given_ps(double p, double s)
  {
    int r = IF97.Region_ps(p/1000000.0, s/1000.0);
    switch(r)
    {
      case 1: //liquid
        return 0; //subcooled liquid
      case 2: //vapor
        return 2; //superheated vapor
      case 3:
      case 4: //two-phase
      case 5:
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
    double ret_val;
    try
    {
      ret_val = IAPWS95.IAPWS95_pressure(1.0 / v, t) * 1000.0; //expects:Kg/M³,K returns KPa
      return ret_val;
    }
    catch (Exception ex)
    {
      Debug.Log( String.Format("Got an exception: {0}\nReturning {1}", ex.Message, p_neutral) );
      got_error = true;
      return p_neutral;
    }
    //return IAPWS95.IAPWS95_pressure(1.0/v,t)*1000.0; //expects:Kg/M³,K returns KPa
  }

  public static double v_given_pt(double p, double t) //DO NOT USE IN VAPOR DOME
  {
    double ret_val;
    try
    {
      ret_val = 1.0/IF97.rhomass_Tp(t,p/1000000.0); //expects:K,MPa returns Kg/M³
      return ret_val;
    }
    catch (Exception ex)
    {
      Debug.Log( String.Format("Got an exception: {0}\nReturning {1}", ex.Message, v_neutral) );
      got_error = true;
      return v_neutral;
    }
    //return 1.0/IF97.rhomass_Tp(t,p/1000000.0); //expects:K,MPa returns Kg/M³
  }

  public static double v_given_ph(double p, double h)
  {
    double ret_val;
    try
    {
      ret_val = 1.0/IF97.rhomass_phmass(p/1000000.0,h/1000.0); //UNIT CONVERSION UNTESTED!
      return ret_val;
    }
    catch (Exception ex)
    {
      Debug.Log( String.Format("Got an exception: {0}\nReturning {1}", ex.Message, v_neutral) );
      got_error = true;
      return v_neutral;
    }
    //return 1.0/IF97.rhomass_phmass(p/1000000.0,h/1000.0); //UNIT CONVERSION UNTESTED!
  }

  public static double v_given_px(double p,double x) //ONLY USE IN VAPOR DOME
  {
    double ret_val;
    try
    {
      ret_val = 1.0/IF97.rhomass_pQ(p/1000000.0,x); //UNIT CONVERSION UNTESTED!
      return ret_val;
    }
    catch (Exception ex)
    {
      Debug.Log( String.Format("Got an exception: {0}\nReturning {1}", ex.Message, v_neutral) );
      got_error = true;
      return v_neutral;
    }
    //return 1.0/IF97.rhomass_pQ(p/1000000.0,x); //UNIT CONVERSION UNTESTED!
  }

  public static double t_given_ph(double p, double h)
  {
    double ret_val;
    try
    {
      ret_val = IF97.T_phmass(p/1000000.0,h/1000.0); //UNIT CONVERSION UNTESTED!
      return ret_val;
    }
    catch (Exception ex)
    {
      Debug.Log( String.Format("Got an exception: {0}\nReturning {1}", ex.Message, t_neutral) );
      got_error = true;
      return t_neutral;
    }
    //return IF97.T_phmass(p/1000000.0,h/1000.0); //UNIT CONVERSION UNTESTED!
  }

  public static double tsat_given_p(double p)
  {
    double ret_val;
    try
    {
      ret_val = IF97.Tsat97(p/1000000.0); //UNIT CONVERSION UNTESTED!
      return ret_val;
    }
    catch (Exception ex)
    {
      Debug.Log( String.Format("Got an exception: {0}\nReturning {1}", ex.Message, t_neutral) );
      got_error = true;
      return t_neutral;
    }
    //return IF97.Tsat97(p/1000000.0); //UNIT CONVERSION UNTESTED!
  }

  public static double vliq_given_p(double p)
  {
    double ret_val;
    try
    {
      ret_val = 1.0/IF97.rholiq_p(p/1000000.0); //expects:MPa returns Kg/M³
      return ret_val;
    }
    catch (Exception ex)
    {
      Debug.Log( String.Format("Got an exception: {0}\nReturning {1}", ex.Message, v_neutral) );
      got_error = true;
      return v_neutral;
    }
    //return 1.0/IF97.rholiq_p(p/1000000.0); //expects:MPa returns Kg/M³
  }

  public static double vvap_given_p(double p)
  {
    double ret_val;
    try
    {
      ret_val = 1.0/IF97.rhovap_p(p/1000000.0); //expects:MPa returns Kg/M³
      return ret_val;
    }
    catch (Exception ex)
    {
      Debug.Log( String.Format("Got an exception: {0}\nReturning {1}", ex.Message, v_neutral) );
      got_error = true;
      return v_neutral;
    }
    //return 1.0/IF97.rhovap_p(p/1000000.0); //expects:MPa returns Kg/M³
  }

  public static double u_given_pt(double p, double t) //DO NOT USE IN VAPOR DOME
  {
    double ret_val;
    try
    {
      ret_val = IF97.umass_Tp(t, p/1000000.0)*1000f; //UNIT CONVERSION UNTESTED!
      return ret_val;
    }
    catch (Exception ex)
    {
      Debug.Log( String.Format("Got an exception: {0}\nReturning {1}", ex.Message, u_neutral) );
      got_error = true;
      return u_neutral;
    }
    //return IF97.umass_Tp(t, p/1000000.0)*1000f; //UNIT CONVERSION UNTESTED!
  }

  public static double u_given_vt(double v, double t)
  {
    double ret_val;
    try
    {
      ret_val = IAPWS95.IAPWS95_internal_energy(1f/v,t)*1000f; //UNIT CONVERSION UNTESTED!
      return ret_val;
    }
    catch (Exception ex)
    {
      Debug.Log( String.Format("Got an exception: {0}\nReturning {1}", ex.Message, u_neutral) );
      got_error = true;
      return u_neutral;
    }
    //return IAPWS95.IAPWS95_internal_energy(1f/v,t)*1000f; //UNIT CONVERSION UNTESTED!
  }

  public static double u_given_px(double p, double x)
  {
    double ret_val;
    try
    {
      ret_val = IF97.umass_pQ(p/1000000.0,x)*1000f; //UNIT CONVERSION UNTESTED!
      return ret_val;
    }
    catch (Exception ex)
    {
      Debug.Log( String.Format("Got an exception: {0}\nReturning {1}", ex.Message, u_neutral) );
      got_error = true;
      return u_neutral;
    }
    //return IF97.umass_pQ(p/1000000.0,x)*1000f; //UNIT CONVERSION UNTESTED!
  }

  public static double s_given_vt(double v, double t)
  {
    double ret_val;
    try
    {
      ret_val = IAPWS95.IAPWS95_entropy(1f/v,t)*1000f; //UNIT CONVERSION UNTESTED!
      return ret_val;
    }
    catch (Exception ex)
    {
      Debug.Log( String.Format("Got an exception: {0}\nReturning {1}", ex.Message, s_neutral) );
      got_error = true;
      return s_neutral;
    }
    //return IAPWS95.IAPWS95_entropy(1f/v,t)*1000f; //UNIT CONVERSION UNTESTED!
  }

  public static double s_given_px(double p, double x)
  {
    double ret_val;
    try
    {
      ret_val = IF97.smass_pQ(p/1000000.0,x)*1000f; //UNIT CONVERSION UNTESTED!
      return ret_val;
    }
    catch (Exception ex)
    {
      Debug.Log( String.Format("Got an exception: {0}\nReturning {1}", ex.Message, s_neutral) );
      got_error = true;
      return s_neutral;
    }
    //return IF97.smass_pQ(p/1000000.0,x)*1000f; //UNIT CONVERSION UNTESTED!
  }

  public static double h_given_vt(double v, double t)
  {
    double ret_val;
    try
    {
      ret_val = IAPWS95.IAPWS95_enthalpy(1f/v,t)*1000f; //UNIT CONVERSION UNTESTED!
      return ret_val;
    }
    catch (Exception ex)
    {
      Debug.Log( String.Format("Got an exception: {0}\nReturning {1}", ex.Message, h_neutral) );
      got_error = true;
      return h_neutral;
    }
    //return IAPWS95.IAPWS95_enthalpy(1f/v,t)*1000f; //UNIT CONVERSION UNTESTED!
  }

  public static double x_given_pv(double p, double v) //ONLY USE IN VAPOR DOME
  {
    double ret_val;
    try
    {
      //f means saturated liquid,
      //g means saturated gas
      double vf = 1.0/IF97.rholiq_p(p/1000000.0);
      double vg = 1.0/IF97.rhovap_p(p/1000000.0);
      ret_val = (v-vf)/(vg-vf); //UNIT CONVERSION UNTESTED!
      return ret_val;
    }
    catch (Exception ex)
    {
      Debug.Log( String.Format("Got an exception: {0}\nReturning {1}", ex.Message, x_neutral) );
      got_error = true;
      return x_neutral;
    }
    //return (v-vf)/(vg-vf); //UNIT CONVERSION UNTESTED!
  }

  public static double x_given_ph(double p, double h) //ONLY USE IN VAPOR DOME
  {
    double ret_val;
    try
    {
      ret_val = IF97.Q_phmass(p/1000000.0,h/1000.0); //UNIT CONVERSION UNTESTED!
      return ret_val;
    }
    catch (Exception ex)
    {
      Debug.Log( String.Format("Got an exception: {0}\nReturning {1}", ex.Message, x_neutral) );
      got_error = true;
      return x_neutral;
    }
    //return IF97.Q_phmass(p/1000000.0,h/1000.0); //UNIT CONVERSION UNTESTED!
  }

  public static double iterate_t_given_p_verify_u(double t, double p, double u) //t = first guess
  {
    try
    {
      var task = Task.Run(() =>
      {
        int MAX_ITERS = 100;     //max # of iterations before giving up //TODO: define intentionally
        double MAX_DELTA = 0.01; //acceptible solution error //TODO: define intentionally
        double step = 0.01;      //size of first step (shrinks every time it overshoots) //TODO: define intentionally
        double guess = t;
        double mark = u;
        double delta = Math.Abs(u_given_pt(p,guess)-mark);
        for(int i = 0; i < MAX_ITERS && delta < MAX_DELTA; i++)
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
      });
      if (!task.Wait(MAX_MILLISECOND))
      {
        Debug.Log( String.Format("Did not complete iterate_t_given_p_verify_u in {0} msec, returning initial guess {1}", MAX_MILLISECOND, t) );
        return t;
      }
      else
      {
        //Debug.Log("Successfully completed iterate_t_given_p_verify_u");
        return task.Result;
      }
    }
    catch (Exception ex)
    {
      Debug.Log( String.Format("Got an exception: {0}\nReturning {1}", ex.Message, t_neutral) );
      got_error = true;
      return t_neutral;
    }
  }

  public static double iterate_t_given_v_verify_u(double t, double v, double u) //t = first guess
  {
    try
    {
      var task = Task.Run(() =>
      {
        int MAX_ITERS = 100;     //max # of iterations before giving up //TODO: define intentionally
        double MAX_DELTA = 0.01; //acceptible solution error //TODO: define intentionally
        double step = 0.01;      //size of first step (shrinks every time it overshoots) //TODO: define intentionally
        double guess = t;
        double mark = u;
        double delta = u_given_vt(v,guess)-mark;
        for(int i = 0; i < MAX_ITERS && delta < MAX_DELTA; i++)
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
      });
      if (!task.Wait(50))
      {
        Debug.Log( String.Format("Did not complete iterate_t_given_v_verify_u in {0} sec, returning initial guess {1}", MAX_MILLISECOND, t) );
        return t;
      }
      else
      {
        //Debug.Log("Successfully completed iterate_t_given_v_verify_u");
        return task.Result;
      }
    }
    catch (Exception ex)
    {
      Debug.Log( String.Format("Got an exception: {0}\nReturning {1}", ex.Message, t_neutral) );
      got_error = true;
      return t_neutral;
    }
  }
}

