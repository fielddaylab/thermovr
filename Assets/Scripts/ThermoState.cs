/*
DOCUMENTATION- phil, 12/16/19
This is the stateful wrapper around ThermoMath.
It represents a "current" state of water. After initialization, the state must always remain consistent.
For this reason, the API consists of applying deltas to some assumed consistent state.
It should be safe to assume that after any method call, the state remains consistent.
*/

using System;
using System.IO;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Collections.Specialized;

public class ThermoState : MonoBehaviour
{
    //state
    //xyz corresponds to vpt (Y = "up")
    public double pressure;       //p //pascals
    public double temperature;    //t //°kelvin
    public double volume;         //v //M³/kg
    public double internalenergy; //u //J/kg
    public double entropy;        //s //J/kgK
    public double enthalpy;       //h //J/kg
    public double quality;        //x //%
    public int region;            //0 subcooled liquid, 1 two-phase, 2 superheated vapor
    
    public double prev_pressure;
    public double prev_temperature;
    public double prev_volume;
    public double prev_internalenergy;
    public double prev_entropy;
    public double prev_enthalpy;
    public double prev_quality;
    public int prev_region;

    //static properties of system
    public double mass = 1; //kg
    public double radius = 0.05; //M
                                 //public double surfacearea = Math.Pow(3.141592*radius,2.0); //M^2 //hardcoded answer below
    public double surfacearea = 0.024674011; //M^2 //hardcoded answer to eqn above
    public double surfacearea_insqr = 38.2447935395871; //in^2 //hardcoded conversion from m^2 to in^2

    public void reset() {
        //ensure consistent state
        pressure = ThermoMath.p_neutral[region];
        temperature = ThermoMath.t_neutral[region];
        //from this point, the rest should be derived!
        quality = ThermoMath.x_neutral[region];
        if (region == ThermoMath.region_twophase)
            volume = ThermoMath.v_given_px(pressure, quality, region);
        else
            volume = ThermoMath.v_given_pt(pressure, temperature, region);
        internalenergy = ThermoMath.u_given_pt(pressure, temperature, region);
        enthalpy = ThermoMath.h_given_vt(volume, temperature, region);
        entropy = ThermoMath.s_given_vt(volume, temperature, region);
        region = ThermoMath.region_given_pvt(pressure, volume, temperature); //should certainly stay the same, as bases were calculated from assumed region

        prev_pressure = -1;
        prev_temperature = -1;
        prev_volume = -1;
        prev_internalenergy = -1;
        prev_entropy = -1;
        prev_enthalpy = -1;
        prev_quality = -1;
        prev_region = -1;
    }

    public void stamp_prev() {
        prev_pressure = pressure;
        prev_temperature = temperature;
        prev_volume = volume;
        prev_internalenergy = internalenergy;
        prev_entropy = entropy;
        prev_enthalpy = enthalpy;
        prev_quality = quality;
    }

    public double Clampd(double v, double min, double max) { if (v < min) return min; if (v > max) return max; return v; }
    public float Clampf(float v, float min, float max) { if (v < min) return min; if (v > max) return max; return v; }
    void clamp_state() {
        if (Double.IsNaN(pressure)) pressure = prev_pressure;
        if (Double.IsNaN(temperature)) temperature = prev_temperature;
        if (Double.IsNaN(volume)) volume = prev_volume;
        if (Double.IsNaN(internalenergy)) internalenergy = prev_internalenergy;
        if (Double.IsNaN(entropy)) entropy = prev_entropy;
        if (Double.IsNaN(enthalpy)) enthalpy = prev_enthalpy;
        if (Double.IsNaN(quality)) quality = prev_quality;
        if (region == -1) region = prev_region;

        //DEBUGGING! [COMMENT OUT IN PROD]
        /*
        double npressure       = Clampd(pressure,       ThermoMath.p_min,ThermoMath.p_max);
        double nvolume         = Clampd(volume,         ThermoMath.v_min,ThermoMath.v_max);
        double ntemperature    = Clampd(temperature,    ThermoMath.t_min,ThermoMath.t_max);
        double ninternalenergy = Clampd(internalenergy, ThermoMath.u_min,ThermoMath.u_max);
        double nentropy        = Clampd(entropy,        ThermoMath.s_min,ThermoMath.s_max);
        double nenthalpy       = Clampd(enthalpy,       ThermoMath.h_min,ThermoMath.h_max);
        double nquality        = Clampd(quality,        ThermoMath.x_min,ThermoMath.x_max);

        if(npressure       != pressure)       Debug.LogFormat("pressure!       {0} clamped to {1}",pressure,npressure);
        if(nvolume         != volume)         Debug.LogFormat("volume!         {0} clamped to {1}",volume,nvolume);
        if(ntemperature    != temperature)    Debug.LogFormat("temperature!    {0} clamped to {1}",temperature,ntemperature);
        if(ninternalenergy != internalenergy) Debug.LogFormat("internalenergy! {0} clamped to {1}",internalenergy,ninternalenergy);
        if(nentropy        != entropy)        Debug.LogFormat("entropy!        {0} clamped to {1}",entropy,nentropy);
        if(nenthalpy       != enthalpy)       Debug.LogFormat("enthalpy!       {0} clamped to {1}",enthalpy,nenthalpy);
        if(nquality        != quality)        Debug.LogFormat("quality!        {0} clamped to {1}",quality,nquality);
        */
        //END DEBUGGING

        /*
            //This isn't worth it: even tiny errors in fixed variables warp results,
            //and this function is devoid of the context to know which warps are appropriate.

            //snap to known region using region-stable perspectives
            if(prev_region != region)
            { //make sure new state snapped to new region!
              switch(region)
              {
                case ThermoMath.region_liquid:  //subcooled liquid
                case ThermoMath.region_vapor:  //superheated vapor
                  //use P,T to fix the rest
                  internalenergy = ThermoMath.u_given_pt(pressure, temperature, region);
                  volume         = ThermoMath.v_given_pt(pressure, temperature, region);
                  enthalpy       = ThermoMath.h_given_vt(volume, temperature, region);
                  entropy        = ThermoMath.s_given_vt(volume, temperature, region);
                  if(region == ThermoMath.region_liquid) quality = 0;
                  if(region == ThermoMath.region_vapor)  quality = 1;
                  break;
                case ThermoMath.region_twophase:  //two-phase region
                  //use P,V to fix the rest
                  //temperature = ThermoMath.t_given_pv();
                  temperature = ThermoMath.tsat_given_p(pressure, region);
                  entropy = ThermoMath.s_given_vt(volume, temperature, region);
                  enthalpy = ThermoMath.h_given_vt(volume, temperature, region);
                  quality = ThermoMath.x_given_pv(pressure, volume, region);
                  internalenergy = ThermoMath.u_given_px(pressure, quality, region);
                  break;
              }
            }
        */

        pressure = Clampd(pressure, ThermoMath.p_min, ThermoMath.p_max);
        volume = Clampd(volume, ThermoMath.v_min, ThermoMath.v_max);
        temperature = Clampd(temperature, ThermoMath.t_min, ThermoMath.t_max);
        internalenergy = Clampd(internalenergy, ThermoMath.u_min, ThermoMath.u_max);
        entropy = Clampd(entropy, ThermoMath.s_min, ThermoMath.s_max);
        enthalpy = Clampd(enthalpy, ThermoMath.h_min, ThermoMath.h_max);
        quality = Clampd(quality, ThermoMath.x_min, ThermoMath.x_max);
    }

    //assume starting/ending point consistent for whole API!

    public void debug_deltas(bool debug_write, StreamWriter debug_file) {
        if (!debug_write) return;
        debug_file.WriteLine("pressure {0} changed to {1} (delta {2})", prev_pressure, pressure, pressure - prev_pressure);
        debug_file.WriteLine("temperature {0} changed to {1} (delta {2})", prev_temperature, temperature, temperature - prev_temperature);
        debug_file.WriteLine("volume {0} changed to {1} (delta {2})", prev_volume, volume, volume - prev_volume);
        debug_file.WriteLine("internalenergy {0} changed to {1} (delta {2})", prev_internalenergy, internalenergy, internalenergy - prev_internalenergy);
        debug_file.WriteLine("entropy {0} changed to {1} (delta {2})", prev_entropy, entropy, entropy - prev_entropy);
        debug_file.WriteLine("enthalpy {0} changed to {1} (delta {2})", prev_enthalpy, enthalpy, enthalpy - prev_enthalpy);
        debug_file.WriteLine("quality {0} changed to {1} (delta {2})", prev_quality, quality, quality - prev_quality);
    }

    public void warp_pv(double p, double v, double t) {
        try {
            pressure = p;
            volume = v;
            temperature = t;

            region = ThermoMath.region_given_pvt(pressure, volume, temperature);
            entropy = ThermoMath.s_given_vt(volume, temperature, region);
            enthalpy = ThermoMath.h_given_vt(volume, temperature, region);

            switch (region) {
                case ThermoMath.region_liquid: quality = 0; internalenergy = ThermoMath.u_given_vt(volume, temperature, region); break;
                case ThermoMath.region_twophase: quality = ThermoMath.x_given_pv(pressure, volume, region); internalenergy = ThermoMath.u_given_px(pressure, quality, region); break;
                case ThermoMath.region_vapor: quality = 1; internalenergy = ThermoMath.u_given_vt(volume, temperature, region); break;
            }
        }
        catch (Exception e) {
            temperature = prev_temperature;
            pressure = prev_pressure;
            volume = prev_volume;
            region = prev_region;
            entropy = prev_entropy;
            enthalpy = prev_enthalpy;
            internalenergy = prev_internalenergy;
        }
        clamp_state();
    }

    public void add_heat_constant_p_per_delta_time(double applied_heat, double insulation_coefficient, double delta_time) {
        try {
            double delta_h = delta_time / mass * (applied_heat * insulation_coefficient);  // time eqtn 6a

            double new_h = enthalpy + delta_h;

            switch (region) {
                case ThermoMath.region_twophase:
                    double new_x = ThermoMath.x_given_ph(pressure, new_h, region);
                    //at this point, we have enough internal state to derive the rest
                    enthalpy = new_h;
                    quality = new_x;
                    volume = ThermoMath.v_given_px(pressure, new_x, region);
                    temperature = ThermoMath.tsat_given_p(pressure, region);
                    entropy = ThermoMath.s_given_px(pressure, new_x, region);
                    internalenergy = ThermoMath.u_given_px(pressure, new_x, region);
                    break;
                case ThermoMath.region_liquid:
                case ThermoMath.region_vapor:
                    //at this point, we have enough internal state to derive the rest
                    enthalpy = new_h;
                    volume = ThermoMath.v_given_ph(pressure, new_h, region);
                    temperature = ThermoMath.t_given_ph(pressure, new_h, region);
                    entropy = ThermoMath.s_given_vt(volume, temperature, region);
                    internalenergy = ThermoMath.u_given_vt(volume, temperature, region);
                    break;
            }

            region = ThermoMath.region_given_pvt(pressure, volume, temperature);
            switch (region) {
                case ThermoMath.region_liquid: quality = 0; break;
                case ThermoMath.region_twophase: quality = ThermoMath.x_given_pv(pressure, volume, region); break;
                case ThermoMath.region_vapor: quality = 1; break;
            }
        }
        catch (Exception e) { }

        clamp_state();
    }

    public void add_heat_constant_v_per_delta_time(double applied_heat, double insulation_coefficient, double delta_time) {
        try {
            double delta_u = delta_time / mass * (applied_heat * insulation_coefficient); // time eqtn 6b

            double new_u = internalenergy + delta_u;
            double new_t = temperature;
            double new_p = pressure;

            if (region != ThermoMath.region_twophase) {
                new_t = ThermoMath.iterate_t_given_v_verify_u(temperature, volume, new_u, region); //try to move t assuming we stay in starting region
                if (region == ThermoMath.region_liquid && new_t > ThermoMath.tsat_given_p(pressure)) //overshot form liquid
                {
                    new_t = ThermoMath.tsat_given_p(pressure);
                    region = ThermoMath.region_twophase;
                }
                else if (region == ThermoMath.region_vapor && new_t < ThermoMath.tsat_given_p(pressure)) //overshot from vapor
                {
                    new_t = ThermoMath.tsat_given_p(pressure);
                    region = ThermoMath.region_twophase;
                }
                else {
                    // remain in liquid or vapor region
                    internalenergy = new_u;
                    temperature = new_t;
                    pressure = ThermoMath.p_given_vt(volume, temperature, region);
                }
            }

            
            // two-phase region
            if (region == ThermoMath.region_twophase) //either newly, or all along
            {
                new_p = ThermoMath.iterate_p_given_vu(pressure, volume, new_u, region); // time eqtn 6
                new_t = ThermoMath.tsat_given_p(new_p);
                internalenergy = new_u;
                pressure = new_p;
                temperature = new_t;
                quality = ThermoMath.x_given_pv(pressure, volume, region);
            }

            enthalpy = ThermoMath.h_given_vt(volume, temperature, region);
            entropy = ThermoMath.s_given_vt(volume, temperature, region);
        }
        catch (Exception e) { }

        clamp_state();
    }

    public void add_pressure_uninsulated_per_delta_time(double p, double delta_time) {
        try {
            double new_p = pressure + p * delta_time;

            //default guess
            double new_u = internalenergy;
            double new_v = volume;
            // Pressure Constrained -> Insulated -> delta pressure (all phases)
            new_v = ThermoMath.v_given_pt(new_p, temperature, region); //ERROR: DO NOT USE IN VAPOR DOME (safe assuming any delta p sufficient to _leave_ vapor dome) TODO: enforce this check
            new_u = ThermoMath.u_given_pt(new_p, temperature, region);
            //at this point, we have enough internal state to derive the rest
            pressure = new_p;
            volume = new_v;
            internalenergy = new_u;
            enthalpy = ThermoMath.h_given_vt(volume, temperature, region);
            entropy = ThermoMath.s_given_vt(volume, temperature, region);
            region = ThermoMath.region_given_pvt(pressure, volume, temperature);
            switch (region) {
                case ThermoMath.region_liquid: quality = 0; break;
                case ThermoMath.region_twophase: quality = ThermoMath.x_given_pv(pressure, volume, region); break;
                case ThermoMath.region_vapor: quality = 1; break;
            }

        }
        catch (Exception e) { }

        clamp_state();
    }

    public void add_pressure_insulated_per_delta_time(double p, double delta_time) {
        try {
            double new_p = pressure + p * delta_time;

            switch (region) {
                case ThermoMath.region_liquid: //subcooled liquid
                case ThermoMath.region_twophase: //two-phase region
                {
                        //AVOID THESE SCENARIOS // Pressure Constrained -> Insulated -> delta pressure (liquid and two-phase)
                    }
                    break;
                case ThermoMath.region_vapor: //superheated vapor
                {
                        //default guess
                        double new_t = temperature;
                        double new_u = internalenergy;
                        double new_v = volume;

                        double k = 1.27;
                        new_v = volume * Math.Pow(pressure / new_p, 1.0 / k);
                        new_u = internalenergy - ((new_p * new_v - pressure * volume) / (1 - k));
                        new_t = ThermoMath.iterate_t_given_p_verify_u(temperature, pressure, new_u, region);
                        //new_t = ThermoMath.iterate_t_given_pv(temperature,pressure,new_v, region);

                        //at this point, we have enough internal state to derive the rest
                        pressure = new_p;
                        volume = new_v;
                        temperature = new_t;
                        internalenergy = new_u;
                        enthalpy = ThermoMath.h_given_vt(volume, temperature, region);
                        entropy = ThermoMath.s_given_vt(volume, temperature, region);
                        region = ThermoMath.region_given_pvt(pressure, volume, temperature);
                    }
                    break;
            }
        }
        catch (Exception e) { }

        clamp_state();
    }

}

