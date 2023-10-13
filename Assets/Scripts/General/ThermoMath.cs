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
using ThermoVR.State;
using UnityEngine;

public static class ThermoMath
{
    private const int MAX_MILLISECOND = 50;
    private const int DEFAULT_ITERS = 50; // The lower the number, the better the framerate (but speed at which an accurate answer is reached is slower)
    private const double DEFAULT_STEP = 10.0;

    //public double surfacearea = Math.Pow(3.141592*radius,2.0); //M^2 //hardcoded answer below
    public static double surfacearea = 0.024674011; //M^2 //hardcoded answer to eqn above
    public static double surfacearea_insqr = 38.2447935395871; //in^2 //hardcoded conversion from m^2 to in^2

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

    //region: 0 subcooled liquid, 1 two-phase, 2 superheated vapor
    public const int region_liquid = 0;
    public const int region_twophase = 1;
    public const int region_vapor = 2;

    //Pa
    public static double p_min;
    public static double p_max;
    public static double[] p_neutral;
    public static double p_smallstep;
    //Pa
    public static double psat_min;
    public static double psat_max;
    //M³/kg
    public static double v_min;
    public static double v_max;
    public static double[] v_neutral;
    public static double v_smallstep;
    //K
    public static double t_min;
    public static double t_max;
    public static double[] t_neutral;
    public static double t_smallstep;
    public static double t_crit;

    //J/kg
    public static double u_min;
    public static double u_max;
    public static double[] u_neutral;
    public static double u_smallstep;
    //J/kgK
    public static double s_min;
    public static double s_max;
    public static double[] s_neutral;
    public static double s_smallstep;
    //J/kg
    public static double h_min;
    public static double h_max;
    public static double[] h_neutral;
    public static double h_smallstep;
    //%
    public static double x_min;
    public static double x_max;
    public static double[] x_neutral;
    public static double x_smallstep;

    // A sort of notifier for users of the class.
    // Whenever an exception occurs, we set this to true.
    // User of the class can then check and react appropriately, setting false after handling the error.
    public static bool got_error;

    public static void Init() {
        got_error = false;
        IF97.initRegions();

        //Pa
        p_min = IF97.get_Pmin() * 1000000.0; // 611.213
        p_max = IF97.get_Pmax() * 1000000.0; // 100000000
        p_neutral = new double[] { 101325.0, 3142, 3142 }; 
        p_smallstep = 1.0;

        //Pa
        psat_min = IF97.get_ptrip() * 1000000.0; // 611.656
        psat_max = IF97.get_pcrit() * 1000000.0; // 22064000

        //M³/kg
        v_min = 1.0 / 3000;  // 0.0003
        v_max = 1.0 / 0.001; // 1000
        v_neutral = new double[] { 0.001, 21.85, 0 }; //note: v_neutral[2] calculated/set below
        v_smallstep = 0.00001;

        //K
        t_min = IF97.get_Tmin(); // 273.15
        t_max = IF97.get_Tmax(); // 1073.15
        t_neutral = new double[] { 293.0, 298.0, 400.0 };
        t_smallstep = 0.01;
        t_crit = IF97.get_Tcrit();

        //J/kg
        u_min = 123.8;
        u_max = 3700000.0;
        u_neutral = new double[] { 83280.0, 8328.0, 8328.0 }; //note: all get calculated/set below
        u_smallstep = 1;

        //J/kgK
        s_min = IF97.Smin * 1000; //0.0
        s_max = IF97.Smax * 1000; //11921.0548250511
        s_neutral = new double[] { 294.322, 294.322, 294.322 };
        s_smallstep = 0.0001;

        //J/kg
        h_min = IF97.Hmin(s_min); //0.0006286
        h_max = IF97.Hmax(s_max / 1000.0) * 1000.0; //4171654.98424024
        h_neutral = new double[] { 83377.0, 83377.0, 83377.0 };
        h_smallstep = 0.001;

        //%
        x_min = 0.0;
        x_max = 1.0;
        x_neutral = new double[] { 0, 0.5, 1 };
        x_smallstep = 0.00001;

        //fill in missing/sloppy neutrals (note: if any of these fallback on "fallback_region", defeats the purpose. do not let that happen)
        v_neutral[2] = v_given_pt(p_neutral[2], t_neutral[2], 2);
        for (int i = 0; i < 3; i++) u_neutral[i] = u_given_vt(v_neutral[i], t_neutral[i], i);
        for (int i = 0; i < 3; i++) s_neutral[i] = s_given_vt(v_neutral[i], t_neutral[i], i);
        for (int i = 0; i < 3; i++) h_neutral[i] = h_given_vt(v_neutral[i], t_neutral[i], i);
        x_neutral[1] = x_given_pv(p_neutral[1], v_neutral[1], 1);

        //DEBUG INFO:
        //IF97.print_tables();
        //IAPWS95.print_tables();
        //compare_impls();
    }

    static double SAMPLE_LBASE = 1.6f;
    static double sample(double t) { return Math.Pow(t, SAMPLE_LBASE); }


    /*
      //FOR DEBUGGING
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
            //pvt in Pa, M³/kg, K
            double _p = p_given_vt(v,t);

            Debug.LogFormat("error:{0} p:{1}Pa ({2}Pa), v:{3}M³/kg, t:{4}K",p-_p,p,_p,v,t);
          }
        }

      }
    */

    static double Lerpd(double a, double b, double t) { return (b - a) * t + a; }

    //a bunch of options for getting region here- still need to figure out most reliable
    //region: 0 subcooled liquid, 1 two-phase, 2 superheated vapor
    public static int region_given_pvt(double p, double v, double t) {
        //broad check w/ t
        if (p <= ThermoMath.psat_max) {
            if (t - t_smallstep > IF97.Tsat97(p / 1000000.0)) return region_vapor;
            if (t + t_smallstep < IF97.Tsat97(p / 1000000.0)) return region_liquid;

            //broad check w/ p - unneeded
            //if(p-p_smallstep > IF97.psat97(t)) return liq;
            //if(p+p_smallstep < IF97.psat97(t)) return vapor;

            //fine check w/ v
            //f means saturated liquid,
            //g means saturated gas
            double vf = 1.0 / IF97.rholiq_p(p / 1000000.0);
            if (v <= vf) return region_liquid;
            double vg = 1.0 / IF97.rhovap_p(p / 1000000.0);
            if (v >= vg) return region_vapor;
            return region_twophase;
        }
        else {
            // TODO: find a way to determing region above p_crit line
            throw new Exception("Unable to determine region given p, v, and t");
        }
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

    // Currently unused; 6/30/2020
    /*
      public static int region_given_ph(double p, double h)
      {
        //Be careful - don’t calculate T from these if you are close to the saturation line (within 25 mK).
        int r = IF97.Region_ph(p/1000000.0, h/1000.0); 
        switch(r)
        {
          case 1: //liquid
            return region_liquid; //subcooled liquid
          case 2: //vapor
            return region_vapor; //superheated vapor
          case 3:
          case 4: //two-phase
          case 5:
            return region_twophase; //two-phase
        }
        return -1;
      }
    */
    public static int region_given_ps(double p, double s) {
        int r = IF97.Region_ps(p / 1000000.0, s / 1000.0);
        switch (r) {
            case 1: //liquid
                return region_liquid; //subcooled liquid
            case 2: //vapor
                return region_vapor; //superheated vapor
            case 3:
            case 4: //two-phase
            case 5:
                return region_twophase; //two-phase
        }
        return -1;
    }

    //helpers to easily generate values "randomly", ensuring we're within the range we care about- NOT PHYSICALLY BASED
    //"percent" = "percent between min value and max value across given dimension"
    public static double p_given_percent(double t) { return Lerpd(p_min, p_max, t); }
    public static double psat_given_percent(double t) { return Lerpd(psat_min, psat_max, t); }
    public static double v_given_percent(double t) { return Lerpd(v_min, v_max, t); }
    public static double t_given_percent(double t) { return Lerpd(t_min, t_max, t); }
    public static double u_given_percent(double t) { return Lerpd(u_min, u_max, t); }
    public static double s_given_percent(double t) { return Lerpd(s_min, s_max, t); }
    public static double h_given_percent(double t) { return Lerpd(h_min, h_max, t); }
    public static double x_given_percent(double t) { return t; } //x already is a percent

    public static double percent_given_p(double p) { return (p - p_min) / (p_max - p_min); }
    public static double percent_given_psat(double psat) { return (psat - psat_min) / (psat_max - psat_min); }
    public static double percent_given_v(double v) { return (v - v_min) / (v_max - v_min); }
    public static double percent_given_t(double t) { return (t - t_min) / (t_max - t_min); }
    public static double percent_given_u(double u) { return (u - u_min) / (u_max - u_min); }
    public static double percent_given_s(double s) { return (s - s_min) / (s_max - s_min); }
    public static double percent_given_h(double h) { return (h - h_min) / (h_max - h_min); }
    public static double percent_given_x(double x) { return x; } //x already is a percent

    //rule of naming for consistency: prefer lexical ordering "p < v < t < u < s < h < q", ie "p_given_vt" rather than "p_given_tv"

    public static double p_given_vt(double v, double t, int fallback_region = 0) //experimentally only valid in the superheated vapor region
    {
        try {
            return IAPWS95.IAPWS95_pressure(1.0 / v, t) * 1000.0; //expects:kg/M³,K returns kPa
        }
        catch (Exception ex) {
            Debug.Log(String.Format("Got an exception: {0}\nReturning {1}", ex.Message, p_neutral[fallback_region]));
            Debug.Log("[Error] " + ex.Message);
            got_error = true;
            return p_neutral[fallback_region];
        }
    }
    public static double v_given_pt(double p, double t, int fallback_region = 0, bool projecting = false) //DO NOT USE IN VAPOR DOME
    {
        try {
            return 1.0 / IF97.rhomass_Tp(t, p / 1000000.0); //expects:K,MPa returns kg/M³
        }
        catch (Exception ex) {
            if (projecting) {
                throw ex;
            }
            else {
                Debug.Log(String.Format("Got an exception: {0}\nReturning {1}", ex.Message, v_neutral[fallback_region]));
                Debug.Log("[Error] " + ex.Message);
                got_error = true;
                return v_neutral[fallback_region];
            }
        }
    }

    public static double v_given_ph(double p, double h, int fallback_region = 0, bool projecting = false) {
        try {
            return 1.0 / IF97.rhomass_phmass(p / 1000000.0, h / 1000.0);
        }
        catch (Exception ex) {
            if (projecting) {
                throw ex;
            }
            else {
                Debug.Log(String.Format("Got an exception: {0}\nReturning {1}", ex.Message, v_neutral[fallback_region]));
                Debug.Log("[Error] " + ex.Message);
                got_error = true;
                return v_neutral[fallback_region];
            }
        }
    }

    /// <summary>
    ///  Same as above, but used for projecting purposes. Does not set error state, since it is only a hypothetical.
    /// </summary>
    /// <param name="p"></param>
    /// <param name="h"></param>
    /// <param name="fallback_region"></param>
    /// <returns></returns>
    public static double v_given_ph_projected(double p, double h, int fallback_region = 0) {
        try {
            return 1.0 / IF97.rhomass_phmass(p / 1000000.0, h / 1000.0);
        }
        catch (Exception ex) {
            Debug.Log(String.Format("Got an exception: {0}\nReturning {1}", ex.Message, v_neutral[fallback_region]));
            Debug.Log("[Error] " + ex.Message);
            throw ex;
        }
    }

    public static double v_given_px(double p, double x, int fallback_region = 0, bool projecting = false) //ONLY USE IN VAPOR DOME
    {
        try {
            return 1.0 / IF97.rhomass_pQ(p / 1000000.0, x);
        }
        catch (Exception ex) {
            if (projecting) {
                throw ex;
            }
            else {
                Debug.Log(String.Format("Got an exception: {0}\nReturning {1}", ex.Message, v_neutral[fallback_region]));
                Debug.Log("[Error] " + ex.Message);
                got_error = true;
                return v_neutral[fallback_region];
            }
        }
    }

    public static double t_given_ph(double p, double h, int fallback_region = 0) {
        try {
            return IF97.T_phmass(p / 1000000.0, h / 1000.0);
        }
        catch (Exception ex) {
            Debug.Log(String.Format("Got an exception: {0}\nReturning {1}", ex.Message, t_neutral[fallback_region]));
            Debug.Log("[Error] " + ex.Message);
            got_error = true;
            return t_neutral[fallback_region];
        }
    }

    public static double tsat_given_p(double p, int fallback_region = 0, bool customHandle = false) {
        try {
            return IF97.Tsat97(p / 1000000.0);
        }
        catch (Exception ex) {
            Debug.Log(String.Format("Got an exception: {0}\nReturning {1}", ex.Message, t_neutral[fallback_region]));
            Debug.Log("[Error] " + ex.Message);
            if (customHandle) {
                throw ex;
            }
            else {
                got_error = true;
                return t_neutral[fallback_region];
            }


        }
    }

    public static double vliq_given_p(double p, int fallback_region = 0) {
        try {
            return 1.0 / IF97.rholiq_p(p / 1000000.0); //expects:MPa returns kg/M³
        }
        catch (Exception ex) {
            Debug.Log(String.Format("Got an exception: {0}\nReturning {1}", ex.Message, v_neutral[fallback_region]));
            Debug.Log("[Error] " + ex.Message);
            got_error = true;
            return v_neutral[fallback_region];
        }
    }

    public static double vvap_given_p(double p, int fallback_region = 0) {
        try {
            return 1.0 / IF97.rhovap_p(p / 1000000.0); //expects:MPa returns kg/M³
        }
        catch (Exception ex) {
            Debug.Log(String.Format("Got an exception: {0}\nReturning {1}", ex.Message, v_neutral[fallback_region]));
            Debug.Log("[Error] " + ex.Message);
            got_error = true;
            return v_neutral[fallback_region];
        }
    }

    public static double u_given_pt(double p, double t, int fallback_region = 0) //DO NOT USE IN VAPOR DOME
    {
        try {
            return IF97.umass_Tp(t, p / 1000000.0) * 1000;
        }
        catch (Exception ex) {
            Debug.Log(String.Format("Got an exception: {0}\nReturning {1}", ex.Message, u_neutral[fallback_region]));
            Debug.Log("[Error] " + ex.Message);
            got_error = true;
            return u_neutral[fallback_region];
        }
    }

    public static double u_given_vt(double v, double t, int fallback_region = 0) {
        try {
            return IAPWS95.IAPWS95_internal_energy(1f / v, t) * 1000;
        }
        catch (Exception ex) {
            Debug.Log(String.Format("Got an exception: {0}\nReturning {1}", ex.Message, u_neutral[fallback_region]));
            Debug.Log("[Error] " + ex.Message);
            got_error = true;
            return u_neutral[fallback_region];
        }
    }

    public static double u_given_px(double p, double x, int fallback_region = 0) {
        try {
            return IF97.umass_pQ(p / 1000000.0, x) * 1000;
        }
        catch (Exception ex) {
            Debug.Log(String.Format("Got an exception: {0}\nReturning {1}", ex.Message, u_neutral[fallback_region]));
            Debug.Log("[Error] " + ex.Message);
            got_error = true;
            return u_neutral[fallback_region];
        }
    }

    public static double s_given_vt(double v, double t, int fallback_region = 0) {
        try {
            return IAPWS95.IAPWS95_entropy(1f / v, t) * 1000f;
        }
        catch (Exception ex) {
            Debug.Log(String.Format("Got an exception: {0}\nReturning {1}", ex.Message, s_neutral[fallback_region]));
            Debug.Log("[Error] " + ex.Message);
            got_error = true;
            return s_neutral[fallback_region];
        }
    }

    public static double s_given_px(double p, double x, int fallback_region = 0) {
        try {
            return IF97.smass_pQ(p / 1000000.0, x) * 1000f;
        }
        catch (Exception ex) {
            Debug.Log(String.Format("Got an exception: {0}\nReturning {1}", ex.Message, s_neutral[fallback_region]));
            Debug.Log("[Error] " + ex.Message);
            throw ex;
            // got_error = true;
            // return s_neutral[fallback_region];
        }
    }

    public static double h_given_vt(double v, double t, int fallback_region = 0) { // DOES NOT APPEAR TO WORK IN VAPOR DOME
        try {
            return IAPWS95.IAPWS95_enthalpy(1f / v, t) * 1000f;
        }
        catch (Exception ex) {
            Debug.Log(String.Format("Got an exception: {0}\nReturning {1}", ex.Message, h_neutral[fallback_region]));
            Debug.Log("[Error] " + ex.Message);
            got_error = true;
            return h_neutral[fallback_region];
        }
    }

    public static double h_given_px(double p, double x, int fallback_region = 0) {
        try {
            return IF97.hmass_pQ(p / 1000000.0, x) * 1000f;
        }
        catch (Exception ex) {
            Debug.Log(String.Format("Got an exception: {0}\nReturning {1}", ex.Message, s_neutral[fallback_region]));
            Debug.Log("[Error] " + ex.Message);
            throw ex;
            // got_error = true;
            // return s_neutral[fallback_region];
        }
    }

    public static double x_given_pv(double p, double v, int fallback_region = 0) //ONLY USE IN VAPOR DOME
    {
        try {
            //f means saturated liquid,
            //g means saturated gas
            double vf = 1.0 / IF97.rholiq_p(p / 1000000.0);
            double vg = 1.0 / IF97.rhovap_p(p / 1000000.0);
            return (v - vf) / (vg - vf);
        }
        catch (Exception ex) {
            Debug.Log(String.Format("Got an exception: {0}\nReturning {1}", ex.Message, x_neutral[fallback_region]));
            Debug.Log("[Error] " + ex.Message);
            got_error = true;
            return x_neutral[fallback_region];
        }
    }

    public static double x_given_ph(double p, double h, int fallback_region = 0) //ONLY USE IN VAPOR DOME
    {
        try {
            return IF97.Q_phmass(p / 1000000.0, h / 1000.0);
        }
        catch (Exception ex) {
            Debug.Log(String.Format("Got an exception: {0}\nReturning {1}", ex.Message, x_neutral[fallback_region]));
            Debug.Log("[Error] " + ex.Message);
            got_error = true;
            return x_neutral[fallback_region];
        }
    }

    public static double x_given_pu(double p, double u, int fallback_region = 0) //ONLY USE IN VAPOR DOME
    {
        try {
            return IF97.Q_pumass(p / 1000000.0, u / 1000.0);
        }
        catch (Exception ex) {
            Debug.Log(String.Format("Got an exception: {0}\nReturning {1}", ex.Message, x_neutral[fallback_region]));
            Debug.Log("[Error] " + ex.Message);
            got_error = true;
            return x_neutral[fallback_region];
        }
    }

    public static double iterate_t_given_p_verify_u(double t, double p, double u, int fallback_region = 0) //t = first guess
    {
        // NOTE: Uses step
        try {
            int MAX_ITERS = DEFAULT_ITERS; //max # of iterations before giving up
            double MAX_DELTA = 0.01; //acceptible solution error
            double step = DEFAULT_STEP; //size of first step (shrinks every time it overshoots)
            double guess = t;
            double mark = u;
            double delta = Math.Abs(u_given_pt(p, guess) - mark);
            int i = 0;
            for (i = 0; i < MAX_ITERS && delta > MAX_DELTA; i++) {
                double delta_a = Math.Abs(u_given_pt(p, guess + step) - mark);
                double delta_b = Math.Abs(u_given_pt(p, guess - step) - mark);
                if (delta < delta_a && delta < delta_b) //original guess superior
                    step = step / 2.0;
                else if (delta_a < delta_b) {
                    delta = delta_a;
                    guess += step;
                }
                else {
                    delta = delta_b;
                    guess -= step;
                }
            }
            //Debug.LogFormat("{0} iters, {1} delta, {2} guess", i, delta, guess);
            return guess;
        }
        catch (Exception ex) {
            Debug.Log(String.Format("Got an exception: {0}\nReturning {1}", ex.Message, t_neutral[fallback_region]));
            Debug.Log("[Error] " + ex.Message);
            got_error = true;
            return t_neutral[fallback_region];
        }
    }

    public static double iterate_t_given_v_verify_u(double t, double v, double u, int fallback_region = 0) //t = first guess
    {
        // NOTE: Uses step
        try {
            int MAX_ITERS = DEFAULT_ITERS; //max # of iterations before giving up
            double MAX_DELTA = 0.0001; //acceptible solution error
            double step = DEFAULT_STEP; //size of first step (shrinks every time it overshoots)
            double guess = t;
            double mark = u;
            double delta = Math.Abs(u_given_vt(v, guess, fallback_region) - mark);
            int i = 0;
            for (i = 0; i < MAX_ITERS && delta > MAX_DELTA; i++) {
                double delta_a = Math.Abs(u_given_vt(v, guess + step, fallback_region) - mark);
                double delta_b = Math.Abs(u_given_vt(v, guess - step, fallback_region) - mark);
                if (delta < delta_a && delta < delta_b) //original guess superior
                    step = step / 2.0;
                else if (delta_a < delta_b) {
                    delta = delta_a;
                    guess += step;
                }
                else {
                    delta = delta_b;
                    guess -= step;
                }
            }
            //Debug.LogFormat("{0} iters, {1} delta, {2} guess", i, delta, guess);
            return guess;
        }
        catch (Exception ex) {
            Debug.Log(String.Format("Got an exception: {0}\nReturning {1}", ex.Message, t_neutral[fallback_region]));
            Debug.Log("[Error] " + ex.Message);
            got_error = true;
            return t_neutral[fallback_region];
        }
    }
    public static double iterate_t_given_pv(double t, double p, double v, int fallback_region = 0, bool projecting = false) //t = first guess
    {
        // NOTE: Uses step
        try {
            int MAX_ITERS = DEFAULT_ITERS; //max # of iterations before giving up
            double MAX_DELTA = 0.01; //acceptible solution error
            double step = DEFAULT_STEP * (p / p_max); //size of first step (shrinks every time it overshoots)
            double guess = t;
            double vdelta = MAX_DELTA + 1.0;
            double pdelta = MAX_DELTA + 1.0;
            int i = 0;
            for (i = 0; i < MAX_ITERS && (vdelta > MAX_DELTA || pdelta > MAX_DELTA); i++) {
                //one iteration on v
                if ((p < ThermoMath.psat_max) && region_given_pvt(p, v, guess) != region_twophase) {
                    try {
                        vdelta = Math.Abs(v_given_pt(p, guess, fallback_region, projecting) - v);
                    }
                    catch (Exception ex) {
                        handle_step_error(ex, projecting, fallback_region);
                    }
                    double vdelta_a = vdelta;
                    double vdelta_b = vdelta;
                    try {
                        vdelta_a = Math.Abs(v_given_pt(p, guess + step, fallback_region, projecting) - v);
                    }
                    catch (Exception ex) {
                        handle_step_error(ex, projecting, fallback_region);
                    }
                    try {
                        vdelta_b = Math.Abs(v_given_pt(p, guess - step, fallback_region, projecting) - v);
                    }
                    catch (Exception ex) {
                        handle_step_error(ex, projecting, fallback_region);
                    }
                    if (vdelta < vdelta_a && vdelta < vdelta_b) //unaltered guess is superior
                        step = step / 6.0 * (p / p_max);
                    else if (vdelta_a < vdelta_b) {
                        vdelta = vdelta_a;
                        guess += step;
                        step *= 2;
                    }
                    else {
                        vdelta = vdelta_b;
                        guess -= step;
                        step *= 2;
                    }
                    i++; //force "iteration counter", bc we do two iters per loop
                }

                //another iteration on p
                try {
                    pdelta = Math.Abs(p_given_vt(v, guess, fallback_region) - p);
                }
                catch (Exception ex) {
                    handle_step_error(ex, projecting, fallback_region);
                }
                double pdelta_a = pdelta;
                double pdelta_b = pdelta;
                try {
                    pdelta_a = Math.Abs(p_given_vt(v, guess + step, fallback_region) - p);
                }
                catch (Exception ex) {
                    handle_step_error(ex, projecting, fallback_region);
                }
                try {
                    pdelta_b = Math.Abs(p_given_vt(v, guess - step, fallback_region) - p);
                }
                catch (Exception ex){
                    handle_step_error(ex, projecting, fallback_region);
                }

                if (pdelta < pdelta_a && pdelta < pdelta_b) //unaltered guess is superior
                {
                    step = step / 2.0;
                }
                else if (pdelta_a < pdelta_b) {
                    pdelta = pdelta_a;
                    guess += step;
                }
                else {
                    pdelta = pdelta_b;
                    guess -= step;
                }
            }
            //Debug.LogFormat("{0} iters, {1} vdelta, {2} pdelta, {3} guess", i, vdelta, pdelta, guess);
            try {
                v_given_pt(p, guess, fallback_region, projecting);
                p_given_vt(v, guess, fallback_region);
            }
            catch (Exception ex) {
                handle_step_error(ex, projecting, fallback_region, true);
            }

            return guess;
        }
        catch (Exception ex) {
            return handle_step_error(ex, projecting, fallback_region);
        }
    }

    static private double handle_step_error(Exception ex, bool projecting, int fallback_region, bool final = false) {
        // only throw errors if the best guess has been finalized
        if (!final && !projecting) { return 0.0; }

        // TODO -- handle boundary check for applying pressure at min temp vs boundary check for max pressure
        if (projecting) {
            // propogate error
            throw ex;
        }
        else {
            Debug.Log(String.Format("Got an exception: {0}\nReturning {1}", ex.Message, t_neutral[fallback_region]));
            Debug.Log("[Error] " + ex.Message);
            got_error = true;
            return t_neutral[fallback_region];
        }
    }

    public static double iterate_p_given_vu(double p, double v, double u, int fallback_region = 0) //pq = first guess, step = first step
    {
        // NOTE: Uses step
        try {
            int MAX_ITERS = DEFAULT_ITERS; //max # of iterations before giving up
            double MAX_DELTA = 0.0005; //acceptible solution error
            double step = DEFAULT_STEP; //size of first step (shrinks every time it overshoots)
            double guess = p;
            double delta = Math.Abs(x_given_pu(guess, u, fallback_region) - x_given_pv(guess, v, fallback_region));
            int i = 0;
            for (i = 0; i < MAX_ITERS && delta > MAX_DELTA; i++) {
                double delta_a = Math.Abs(x_given_pu(guess + step, u, fallback_region) - x_given_pv(guess + step, v, fallback_region));
                double delta_b = Math.Abs(x_given_pu(guess - step, u, fallback_region) - x_given_pv(guess - step, v, fallback_region));
                if (delta < delta_a && delta < delta_b) //original guess superior
                    step = step / 2.0;
                else if (delta_a < delta_b) {
                    delta = delta_a;
                    guess += step;
                }
                else {
                    delta = delta_b;
                    guess -= step;
                }
            }
            //Debug.LogFormat("{0} iters, {1} delta, {2} guess", i, delta, guess);
            return guess;
        }
        catch (Exception ex) {
            Debug.Log(String.Format("Got an exception: {0}\nReturning {1}", ex.Message, t_neutral[fallback_region]));
            Debug.Log("[Error] " + ex.Message);
            got_error = true;
            return t_neutral[fallback_region];
        }
    }
}

