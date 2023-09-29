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
using Oculus.Platform;
using ThermoVR.Tools;
using ThermoVR;

namespace ThermoVR.State
{
    /// <summary>
    /// Unique identifiers for simulation variables.
    /// </summary>
    public enum VarID : byte
    {
        Region,
        Pressure,
        Temperature,
        Volume,
        InternalEnergy,
        Entropy,
        Enthalpy,
        Quality,
        VolumeStop
    }

    /// <summary>
    /// Struct for updating a simulation variable readout
    /// </summary>
    public struct VarUpdate
    {
        public VarID ID;
        public string NewText;

        public VarUpdate(VarID id, string newText) {
            ID = id;
            NewText = newText;
        }
    }

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

        public double v_stop1; // volume stop specified by tool_stop1
        public double v_stop2; // volume stop specified by tool_stop2

        List<VolumeStop> v_stops;
        private static float STOP_BUFFER = 0.005f;
        private static double LIQ_2_DIVISION = 0.002;

        public void reset() {
            prev_region = region;
            region = 0;
            //ensure consistent state
            pressure = ThermoMath.p_neutral[region];
            temperature = ThermoMath.t_neutral[region];
            //from this point, the rest should be derived!
            quality = ThermoMath.x_neutral[region];
            if (region == ThermoMath.region_twophase) {
                volume = ThermoMath.v_given_px(pressure, quality, region);
                try {
                    enthalpy = ThermoMath.h_given_px(pressure, quality, region);
                    entropy = ThermoMath.s_given_px(pressure, quality, region);
                }
                catch (Exception ex) {
                    ThermoMath.got_error = true;
                }
            }
            else {
                volume = ThermoMath.v_given_pt(pressure, temperature, region);
                enthalpy = ThermoMath.h_given_vt(volume, temperature, region);
                entropy = ThermoMath.s_given_vt(volume, temperature, region);
            }
            internalenergy = ThermoMath.u_given_pt(pressure, temperature, region);
            region = ThermoMath.region_given_pvt(pressure, volume, temperature); //should certainly stay the same, as bases were calculated from assumed region

            prev_pressure = -1;
            prev_temperature = -1;
            prev_volume = -1;
            prev_internalenergy = -1;
            prev_entropy = -1;
            prev_enthalpy = -1;
            prev_quality = -1;
            prev_region = region;

            v_stops = new List<VolumeStop>();

            GameMgr.Events.Dispatch(GameEvents.WarpPVT);
        }

        #region Stops (Clamps)

        /// <summary>
        /// Remove clamp volume bounds
        /// </summary>
        public void release_v_stop(Tool source) {
            if (!StopExists(source)) {
                return;
            }

            RemoveStop(source);
        }

        /// <summary>
        /// Set a new clamp stop
        /// </summary>
        public void add_v_stop(double v_stop, Tool source) {
            if (StopExists(source)) {
                return;
            }

            // constrain stop's values to global bounds
            v_stop = Clampd(v_stop, ThermoMath.v_min, ThermoMath.v_max);

            VolumeStop new_stop = new VolumeStop(v_stop, source);
            v_stops.Add(new_stop);
        }

        private bool StopExists(Tool source) {
            for (int i = 0; i < v_stops.Count; i++) {
                if (v_stops[i].Source == source) {
                    return true;
                }
            }
            return false;
        }

        private void RemoveStop(Tool source) {

            for (int i = 0; i < v_stops.Count; i++) {
                if (v_stops[i].Source == source) {
                    v_stops.RemoveAt(i);
                    return;
                }
            }
        }

        public void update_v_stop(double v_stop_val, Tool source) {
            for (int i = 0; i < v_stops.Count; i++) {
                if (v_stops[i].Source == source) {
                    VolumeStop temp_stop = v_stops[i];
                    temp_stop.Volume = v_stop_val;
                    v_stops[i] = temp_stop;
                    return;
                }
            }
        }


        #endregion // Stops (Clamps)

        #region Enforce State

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

            enthalpy = ClampEnthalpy(enthalpy, pressure);

            quality = Clampd(quality, ThermoMath.x_min, ThermoMath.x_max);
        }

        private double ClampEnthalpy(double new_h, double p) {
            enthalpy = Clampd(enthalpy, ThermoMath.h_min, ThermoMath.h_max);

            double h_min_given_p;
            double h_max_given_p;
            try {
                IF97.h_bounds_given_p(p, out h_min_given_p, out h_max_given_p);
                new_h = Clampd(new_h, h_min_given_p, h_max_given_p);
            }
            catch (Exception e) { }

            return new_h;
        }

        #endregion // Enforce State

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

                // iterative_weight = ambient_pressure;

                region = ThermoMath.region_given_pvt(pressure, volume, temperature);
                if (region == ThermoMath.region_twophase) {
                    try {
                        entropy = ThermoMath.s_given_px(pressure, quality, region);
                        enthalpy = ThermoMath.h_given_px(pressure, quality, region);
                    }
                    catch (Exception ex) {
                        ThermoMath.got_error = true;
                    }
                }
                else {
                    entropy = ThermoMath.s_given_vt(volume, temperature, region);
                    enthalpy = ThermoMath.h_given_vt(volume, temperature, region);
                }

                // region = ThermoMath.region_given_ps(pressure, entropy);

                switch (region) {
                    case ThermoMath.region_liquid: {
                            quality = 0;
                            internalenergy = ThermoMath.u_given_vt(volume, temperature, region);
                            break;
                        }
                    case ThermoMath.region_twophase: {
                            /*
                            quality = ThermoMath.x_given_pv(pressure, volume, region);
                            entropy = ThermoMath.s_given_px(pressure, quality, region);
                            internalenergy = ThermoMath.u_given_px(pressure, quality, region);
                            */
                            // We've never quite gotten warp_pv to work in two-phase. So for now, disallow it.
                            // TODO: revisit this with new enthalpy calculation
                            reset();
                            break;
                        }
                    case ThermoMath.region_vapor: {
                            quality = 1;
                            internalenergy = ThermoMath.u_given_vt(volume, temperature, region);
                            break;
                        }
                }

                GameMgr.Events.Dispatch(GameEvents.WarpPVT);
            }
            catch (Exception e) {
                reset();
            }
            clamp_state();
        }

        private void revert_state() {
            temperature = prev_temperature;
            pressure = prev_pressure;
            volume = prev_volume;
            region = prev_region;
            entropy = prev_entropy;
            enthalpy = prev_enthalpy;
            internalenergy = prev_internalenergy;
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


        #region Apply Heat/Pressure

        public void add_heat_per_delta_time(double applied_heat, double insulation_coefficient, double delta_time, double p_outside, bool is_internal, double temperature_gradient) {

            // prevent errors from temperature reaching the boundaries of the simulation
            if (temperature_bounded(applied_heat, insulation_coefficient, delta_time)) {
                return;
            }

            // determine whether currently applied tools force a constant volume
            bool constant_v = treat_as_constant_v_add_heat(applied_heat, insulation_coefficient, delta_time, p_outside);

            if (constant_v) {
                add_heat_constant_v_per_delta_time(applied_heat, insulation_coefficient, delta_time);
            }
            else {
                add_heat_constant_p_per_delta_time(applied_heat, insulation_coefficient, delta_time);

                /*
                if (!is_internal && insulation_coefficient == 1) {
                    // For now, only makes sense for external (internal would be multiplied by 0, rendering calculation pointless)
                    Debug.Log("Constant temperature");

                    // case constant_t, using entropy instead of u or h
                    add_heat_constant_t_per_delta_time(applied_heat, insulation_coefficient, delta_time, temperature_gradient);
                }
                else {
                    Debug.Log("Constant pressure");

                    add_heat_constant_p_per_delta_time(applied_heat, insulation_coefficient, delta_time);
                }
                */
            }
        }

        private void add_heat_constant_p_per_delta_time(double applied_heat, double insulation_coefficient, double delta_time) {  // Pressure Constrained -> Insulated && Uninsulated -> delta energy    // time eqtn 6b
            try {
                double delta_h = delta_time / mass * (applied_heat * insulation_coefficient);  // time eqtn 6a

                double new_h = enthalpy + delta_h;
                new_h = ClampEnthalpy(new_h, pressure);
                double new_v = volume;

                switch (region) {
                    case ThermoMath.region_liquid:
                    case ThermoMath.region_vapor:
                        // check that h is within bounds
                        new_h = ClampEnthalpy(new_h, pressure);
                        clamp_state();
                        //at this point, we have enough internal state to derive the rest
                        try {
                            new_v = ThermoMath.v_given_ph(pressure, enthalpy, region, true);
                        }
                        catch (Exception e) {
                            new_v = ThermoMath.v_given_ph(pressure, new_h, region);
                        }

                        if (region == ThermoMath.region_liquid) {
                            if (pressure < ThermoMath.psat_max) {
                                if (new_v > LIQ_2_DIVISION) {
                                    region = ThermoMath.region_twophase;
                                    break;
                                }
                            }
                            else {
                                if (ThermoMath.t_given_ph(pressure, enthalpy, region) >= ThermoMath.t_crit) {
                                    // Technically supercritical fluid. But we don't have that.
                                    region = ThermoMath.region_vapor;
                                }
                            }
                        }
                        else if (region == ThermoMath.region_vapor) {
                            if (pressure >= ThermoMath.psat_max) {
                                if (ThermoMath.t_given_ph(pressure, enthalpy, region) < ThermoMath.t_crit) {
                                    region = ThermoMath.region_liquid;
                                }
                            }

                        }

                        enthalpy = new_h;
                        temperature = ThermoMath.t_given_ph(pressure, enthalpy, region);
                        entropy = ThermoMath.s_given_vt(new_v, temperature, region);
                        internalenergy = ThermoMath.u_given_vt(new_v, temperature, region);
                        break;
                }

                if (region == ThermoMath.region_twophase) {
                    double new_x = ThermoMath.x_given_ph(pressure, new_h, region);

                    //at this point, we have enough internal state to derive the rest
                    clamp_state();
                    new_v = ThermoMath.v_given_px(pressure, new_x, region);
                    //new_v = v_with_enforced_stops(new_v); // enforce volume stops
                    temperature = ThermoMath.tsat_given_p(pressure, region);
                    entropy = ThermoMath.s_given_px(pressure, new_x, region);
                    internalenergy = ThermoMath.u_given_px(pressure, new_x, region);
                    enthalpy = new_h;
                    quality = new_x;
                }

                volume = new_v;

                int prev_region = region;
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

        private void add_heat_constant_v_per_delta_time(double applied_heat, double insulation_coefficient, double delta_time) { // Volume Constrained -> Insulated && Uninsulated ->  delta energy     // time eqtn 6a
            try {
                double delta_u = delta_time / mass * (applied_heat * insulation_coefficient); // time eqtn 6b

                double new_u = internalenergy + delta_u;
                double new_t = temperature;

                if (region != ThermoMath.region_twophase) {
                    new_t = ThermoMath.iterate_t_given_v_verify_u(temperature, volume, new_u, region); //try to move t assuming we stay in starting region
                    if (region == ThermoMath.region_liquid && pressure < ThermoMath.psat_max && new_t > ThermoMath.tsat_given_p(pressure)) //overshot from liquid
                    {
                        new_t = ThermoMath.tsat_given_p(pressure);
                        region = ThermoMath.region_twophase;
                    }
                    else if (region == ThermoMath.region_vapor && pressure < ThermoMath.psat_max && new_t < ThermoMath.tsat_given_p(pressure)) //overshot from vapor
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
                    double new_p = ThermoMath.iterate_p_given_vu(pressure, volume, new_u, region); // time eqtn 6
                    new_t = ThermoMath.tsat_given_p(new_p);
                    internalenergy = new_u;
                    pressure = new_p;
                    temperature = new_t;
                    quality = ThermoMath.x_given_pv(pressure, volume, region);
                    // TODO: figure out why h_given_vt is so off
                    /*
                    enthalpy = ThermoMath.h_given_vt(volume, temperature, region);
                    entropy = ThermoMath.s_given_vt(volume, temperature, region);
                    */
                    try {
                        enthalpy = ThermoMath.h_given_px(pressure, quality, region);
                        entropy = ThermoMath.s_given_px(pressure, quality, region);
                    }
                    catch (Exception e) {
                        // quality out of range
                        enthalpy = ThermoMath.h_given_vt(volume, temperature, region);
                        entropy = ThermoMath.s_given_vt(volume, temperature, region);
                        if (quality < 0) {
                            region = ThermoMath.region_liquid;
                        }
                        else {
                            region = ThermoMath.region_vapor;
                        }
                    }
                }
                else {
                    enthalpy = ThermoMath.h_given_vt(volume, temperature, region);
                    entropy = ThermoMath.s_given_vt(volume, temperature, region);
                }
            }
            catch (Exception e) { }

            clamp_state();
        }

        private void add_heat_constant_t_per_delta_time(double applied_heat, double insulation_coefficient, double delta_time, double temperature_gradient) {  // Temperature Constrained (0% insulation)   // time eqtn 5
            try {
                double new_s = entropy;
                double new_h = enthalpy;

                double delta_t = temperature_gradient * 1000;

                double delta_s = delta_time / mass * ((delta_t * insulation_coefficient) / temperature); // time eqtn 5
                new_s = entropy + delta_s;

                // TODO: verify below
                double delta_h = delta_time / mass * (delta_t * insulation_coefficient);  // time eqtn 6a
                new_h = enthalpy + delta_h;

                new_h = ClampEnthalpy(new_h, pressure);
                double new_v = volume;

                switch (region) {
                    case ThermoMath.region_liquid:
                    case ThermoMath.region_vapor:
                        // check that h is within bounds
                        new_h = ClampEnthalpy(new_h, pressure);
                        clamp_state();
                        //at this point, we have enough internal state to derive the rest
                        try {
                            new_v = ThermoMath.v_given_ph(pressure, enthalpy, region, true);
                        }
                        catch (Exception e) {
                            new_v = ThermoMath.v_given_ph(pressure, new_h, region);
                        }

                        if (region == ThermoMath.region_liquid) {
                            if (pressure < ThermoMath.psat_max) {
                                if (new_v > LIQ_2_DIVISION) {
                                    region = ThermoMath.region_twophase;
                                    break;
                                }
                            }
                            else {
                                if (ThermoMath.t_given_ph(pressure, enthalpy, region) >= ThermoMath.t_crit) {
                                    // Technically supercritical fluid. But we don't have that.
                                    region = ThermoMath.region_vapor;
                                }
                            }
                        }
                        else if (region == ThermoMath.region_vapor) {
                            if (pressure >= ThermoMath.psat_max) {
                                if (ThermoMath.t_given_ph(pressure, enthalpy, region) < ThermoMath.t_crit) {
                                    region = ThermoMath.region_liquid;
                                }
                            }

                        }

                        enthalpy = new_h;
                        temperature = ThermoMath.t_given_ph(pressure, enthalpy, region);
                        internalenergy = ThermoMath.u_given_vt(new_v, temperature, region);
                        break;
                }

                if (region == ThermoMath.region_twophase) {
                    double new_x = ThermoMath.x_given_ph(pressure, new_h, region);

                    //at this point, we have enough internal state to derive the rest
                    clamp_state();
                    new_v = ThermoMath.v_given_px(pressure, new_x, region);
                    //new_v = v_with_enforced_stops(new_v); // enforce volume stops
                    temperature = ThermoMath.tsat_given_p(pressure, region);
                    internalenergy = ThermoMath.u_given_px(pressure, new_x, region);
                    enthalpy = new_h;
                    quality = new_x;
                }

                volume = new_v;

                int prev_region = region;
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

        public void add_pressure_uninsulated_per_delta_time(double p, double delta_time, double insulation_coefficient, double p_outside, double temperature_gradient) {
            if (blocked_by_stops(p_outside)) {
                return;
            }

            double added_p = p;
            double new_p = pressure + added_p; // * delta_time;

            double p_dif = 0;
            if (new_p < ThermoMath.psat_min) {
                p_dif = new_p - ThermoMath.psat_min;
                new_p = ThermoMath.psat_min;
            }

            double iterative_dif = (added_p - p_dif);

            if (region != ThermoMath.region_twophase && enthalpy_bounded(new_p, enthalpy)) {
                return;
            }
            if (region == ThermoMath.region_liquid /*AND temperature is near lower edge*/) {
                try {
                    double new_t = ThermoMath.iterate_t_given_pv(temperature, new_p, volume, region, true); // new_v is constant in liquid
                }
                catch (ArgumentOutOfRangeException e) {
                    // temperature out of range along lower bound
                    return;
                }
            }
            if (treat_as_constant_v_add_p_uninsulated(new_p, added_p, iterative_dif, delta_time, insulation_coefficient)) {
                return;
            }
            if (new_p < ThermoMath.psat_min) {
                // lower bound
                return;
            }

            try {
                //default guess
                double new_u = internalenergy;
                double new_v = volume;

                double new_s = entropy;
                double new_h = enthalpy;

                double new_x = quality;

                switch (region) {
                    case ThermoMath.region_liquid: //subcooled liquid
                    case ThermoMath.region_vapor: //vapor region
                        if (ThermoMath.region_given_ps(new_p, entropy) != region) { // check for transition into 2-phase from other states
                            if (ThermoMath.region_given_ps(new_p, entropy) == ThermoMath.region_twophase) {
                                region = ThermoMath.region_twophase;
                                break;
                            }
                        }

                        //default guess
                        new_u = internalenergy;
                        new_v = volume;

                        new_s = entropy;
                        new_h = enthalpy;

                        if (region == ThermoMath.region_vapor) {
                            double k = 1.27;
                            new_v = volume * Math.Pow(pressure / new_p, 1.0 / k);
                            update_vapor_vis(pressure - new_p, insulation_coefficient);
                            new_u = internalenergy - ((new_p * new_v - pressure * volume) / (1 - k));
                        }

                        double old_t = temperature;
                        double new_t = ThermoMath.iterate_t_given_pv(temperature, new_p, new_v, region); // new_v is constant in liquid
                        double delta_t = new_t - old_t;
                        new_t = old_t + delta_t * (1 - insulation_coefficient); // when insulation is 0%, T = T_old

                        double delta_h = delta_time / mass * (delta_t * insulation_coefficient);  // time eqtn 6a
                        new_h = enthalpy + delta_h;

                        new_v = ThermoMath.v_given_ph(new_p, new_h); // reconstrain v

                        if (region == ThermoMath.region_liquid) {
                            // new_t = ThermoMath.iterate_t_given_pv(temperature, new_p, new_v, region); // new_v is constant in liquid

                            new_v = ThermoMath.v_given_pt(new_p, new_t, region);
                            new_u = ThermoMath.u_given_pt(new_p, new_t, region);

                            if (new_p < ThermoMath.psat_max) {
                                if (new_v > LIQ_2_DIVISION) {
                                    region = ThermoMath.region_twophase;
                                    break;
                                }
                            }
                            else {
                                if (ThermoMath.t_given_ph(new_p, new_h, region) >= ThermoMath.t_crit) {
                                    // Technically supercritical fluid. But we don't have that.
                                    region = ThermoMath.region_vapor;
                                }
                            }
                        }
                        else if (region == ThermoMath.region_vapor) {
                            bool remain_vapor = true;

                            if (remain_vapor) {
                                // new_v = ThermoMath.v_given_ph(new_p, enthalpy);
                                update_vapor_vis(pressure - new_p, insulation_coefficient);
                                new_u = ThermoMath.u_given_vt(new_v, new_t, region);
                            }
                        }

                        //at this point, we have enough internal state to derive the rest

                        new_s = ThermoMath.s_given_vt(new_v, new_t, region);
                        new_h = ThermoMath.h_given_vt(new_v, new_t, region);

                        pressure = new_p;
                        volume = new_v;
                        internalenergy = new_u;
                        temperature = new_t;

                        entropy = new_s;
                        enthalpy = new_h;

                        if (region == ThermoMath.region_vapor) {
                            region = ThermoMath.region_given_ps(pressure, entropy);
                            // TODO: hangs along two-phase line
                        }
                        else {
                            region = ThermoMath.region_given_pvt(pressure, volume, temperature);
                        }

                        switch (region) {
                            case ThermoMath.region_liquid: { quality = 0; break; }
                            case ThermoMath.region_twophase: quality = ThermoMath.x_given_pv(pressure, volume, region); break;
                            case ThermoMath.region_vapor: quality = 1; break;
                        }

                        break;
                    default:
                        break;
                }

                if (region == ThermoMath.region_twophase) //two-phase region, either newly or all along
                 {
                    /* 
                     * A starting point in the 2-phase region must be in terms of P (or T) and some other parameter, such as v, h, u, or s.
                     * The starting point in the 2 phase region cannot be defined by P and T alone.
                     * 
                     * If this condition (uninsulated) falls in the category that the conductivity of the wall is infinity and therefore the temperature is really constant,
                     * then the new state is defined by the conditions:
                     * 
                     * P = P_new and T = T_old.
                     * 
                     * From a starting position in the two-phase region, an increase of pressure with T=constant will send the state over to the sub cooled liquid range.
                     * The process will move along a T=constant, P=constant line until it hits the saturated liquid line.
                     * 
                     * If the conductivity of the wall is less than infinity (insulation % > 0), this same process will occur, although at a slower rate. 
                     */

                    /*
                     * Any percent other than 100% insulation follows a combination of time eqtns formulas 2 and 5
                     * 
                     * Formula 2 is insulation coefficient * delta_t
                     */

                    double old_t = temperature;
                    // inside - outside
                    double new_t = ThermoMath.tsat_given_p(new_p);
                    double delta_t = temperature_gradient;

                    double delta_s = delta_time / mass * ((delta_t * insulation_coefficient) / temperature); // time eqtn 5
                    new_s = entropy + delta_s;

                    double delta_h = delta_time / mass * (delta_t * insulation_coefficient);  // time eqtn 6a
                    new_h = enthalpy + delta_h;

                    // from this point, we have enough internal state to derive the rest

                    new_x = ThermoMath.x_given_ph(new_p, new_h);

                    if (new_x <= 0.0f) {
                        region = ThermoMath.region_liquid;
                        new_v = ThermoMath.vliq_given_p(new_p, region);
                        new_u = ThermoMath.u_given_px(new_p, new_x, region);
                    }
                    else if (new_x > 1.0f) {
                        // TODO: this sticks to the edge of two-phase and vapor instead of transitioning

                        region = ThermoMath.region_vapor;
                    }
                    else {
                        new_v = ThermoMath.v_given_px(new_p, new_x, region);
                        new_u = ThermoMath.u_given_px(new_p, new_x, region);
                    }

                    pressure = new_p;
                    enthalpy = new_h;
                    // enthalpy = ThermoMath.h_given_vt(new_v, new_t, region);

                    try {
                        enthalpy = ThermoMath.h_given_px(new_p, new_x, region);
                    }
                    catch (Exception ex) {
                        ThermoMath.got_error = true;
                    }

                    entropy = new_s;
                    temperature = new_t;
                    internalenergy = new_u;
                    volume = new_v;
                    quality = new_x;

                    region = ThermoMath.region_given_pvt(pressure, volume, temperature);
                }
            }
            catch (Exception e) { }

            clamp_state();
        }

        private bool blocked_by_stops(double p_outside) {
            if (region == 0) { return false; } // liquid volume doesn't change

            bool within_vstop_buffer = false;

            for (int i = 0; i < v_stops.Count; i++) {
                within_vstop_buffer = false;

                VolumeStop curr_stop = v_stops[i];
                double compare_v = curr_stop.Volume;

                if (volume >= compare_v - STOP_BUFFER && volume <= compare_v + STOP_BUFFER) {
                    within_vstop_buffer = true;
                }

                if (within_vstop_buffer) {
                    if (volume < compare_v && pressure >= p_outside) {
                        return true;
                    }
                    else if (volume > compare_v && pressure <= p_outside) {
                        return true;
                    }
                }
            }

            return false;
        }

        public void add_pressure_insulated_per_delta_time(double p, double delta_time, double p_outside, double temperature_gradient) {
            if (blocked_by_stops(p_outside)) {
                return;
            }

            double new_p = pressure + p; // * delta_time;
                                         // if (Math.Abs(p * delta_time) < World.DELTA_PRESSURE_CUTOFF) { new_p = pressure + p; } // small enough step; finish transition

            if (enthalpy_bounded(new_p, enthalpy) && region != ThermoMath.region_twophase) {
                return;
            }
            if (region == ThermoMath.region_liquid /*AND temperature is near lower edge*/) {
                try {
                    double new_t = ThermoMath.iterate_t_given_pv(temperature, new_p, volume, region, true); // new_v is constant in liquid
                }
                catch (ArgumentOutOfRangeException e) {
                    // temperature out of range along lower bound
                    return;
                }
            }
            if (treat_as_constant_v_add_p_insulated(p, delta_time)) {
                return;
            }
            if (new_p < ThermoMath.psat_min) {
                // lower bound
                return;
            }

            double insulation_coefficient = 1; // 100%

            try {
                double new_h = enthalpy;
                double new_u = internalenergy;
                double new_t = temperature;
                double new_x = quality;
                double new_v = volume;
                // entropy = entropy // s = constant with 100% insulation

                switch (region) {
                    case ThermoMath.region_liquid: //subcooled liquid
                    case ThermoMath.region_twophase: //two-phase region
                    {
                            // Pressure Constrained -> Insulated -> delta pressure (liquid and two-phase)
                            //at this point, we have enough internal state to derive the rest

                            double old_t = temperature;

                            if (region == ThermoMath.region_twophase) {
                                new_t = ThermoMath.tsat_given_p(new_p, region);
                            }
                            else {
                                new_t = ThermoMath.iterate_t_given_pv(temperature, new_p, new_v, region); // new_v is constant in liquid
                            }

                            double delta_t = new_t - old_t;
                            double delta_h = delta_time / mass * (delta_t * insulation_coefficient);  // time eqtn 6a
                            new_h = enthalpy + delta_h;

                            if (region == ThermoMath.region_liquid) {
                                new_v = ThermoMath.v_given_pt(new_p, new_t, region);
                                new_u = ThermoMath.u_given_pt(new_p, new_t, region);

                                if (new_p < ThermoMath.psat_max) {
                                    // Below critical point
                                    if (new_v > LIQ_2_DIVISION || (ThermoMath.region_given_pvt(new_p, new_v, new_t) == ThermoMath.region_vapor && new_p < ThermoMath.psat_max)) {
                                        // Far enough out of liquid to be in two phase
                                        region = ThermoMath.region_twophase;
                                        new_t = ThermoMath.tsat_given_p(new_p, region);
                                    }
                                }
                                else {
                                    // Above critical point
                                    if (ThermoMath.t_given_ph(new_p, new_h, region) >= ThermoMath.t_crit) {
                                        // Far enough out to be vapor. Technically supercritical fluid. But we don't have that.
                                        region = ThermoMath.region_vapor;
                                        break;
                                    }
                                }
                            }

                            // either newly or all along
                            if (region == ThermoMath.region_twophase) {
                                new_x = ThermoMath.x_given_ph(new_p, new_h, region);
                                if (new_x <= 0.0f) {
                                    Debug.Log("[weight] two-phase to liquid");
                                    region = ThermoMath.region_liquid;
                                    new_v = ThermoMath.vliq_given_p(new_p, region);
                                    new_u = ThermoMath.u_given_px(new_p, new_x, region);
                                }
                                else if (new_x > 1.0f) {
                                    // TODO: this sticks to the edge of two-phase and vapor instead of transitioning
                                    Debug.Log("[weight] two-phase to vapor");
                                    region = ThermoMath.region_vapor;
                                }
                                else {
                                    new_v = ThermoMath.v_given_px(new_p, new_x, region);
                                    new_u = ThermoMath.u_given_px(new_p, new_x, region);
                                }
                            }

                            pressure = new_p;
                            if (region == ThermoMath.region_twophase) {
                                enthalpy = new_h;
                            }
                            else {
                                enthalpy = ThermoMath.h_given_vt(new_v, new_t, region);
                            }
                            temperature = new_t;
                            internalenergy = new_u;
                            volume = new_v;
                            quality = new_x;

                            region = ThermoMath.region_given_pvt(pressure, volume, temperature);
                        }
                        break;
                }
                if (region == ThermoMath.region_vapor) {
                    helper_add_p_insulated_vapor(new_t, new_u, new_v, new_p, insulation_coefficient);
                }
            }
            catch (Exception e) { }

            clamp_state();
        }

        #endregion // Apply Heat/PressureM

        #region Helpers

        private void helper_add_p_insulated_vapor(double new_t, double new_u, double new_v, double new_p, double insulation_coefficient) {
            Debug.Log("[weight] calculating vapor");
            //default guess
            new_t = temperature;
            new_u = internalenergy;

            double k = 1.27;
            new_v = volume * Math.Pow(pressure / new_p, 1.0 / k);
            update_vapor_vis(pressure - new_p, insulation_coefficient);
            new_u = internalenergy - ((new_p * new_v - pressure * volume) / (1 - k));
            new_t = ThermoMath.iterate_t_given_pv(temperature, new_p, new_v, region);

            //at this point, we have enough internal state to derive the rest
            pressure = new_p;
            volume = new_v;
            temperature = new_t;
            internalenergy = new_u;
            enthalpy = ThermoMath.h_given_vt(volume, temperature, region);
            // entropy = ThermoMath.s_given_vt(volume, temperature, region); // s = constant with 100% insulation
            region = ThermoMath.region_given_pvt(pressure, volume, temperature);
        }

        #endregion // Helpers

        #region Projections

        /// <summary>
        /// Probes the simulation by projecting change in heat in order to check if treating it as constant volume would be more accurate.
        /// When modifying this equation, it is recommended to finalize the ordering of operations in the add_p_insulated() function,
        /// then replicate that order in this function, terminating when the new volume would be set.
        /// </summary>
        /// <param name="applied_heat"></param>
        /// <param name="insulation_coefficient"></param>
        /// <param name="delta_time"></param>
        /// <returns></returns>
        private bool treat_as_constant_v_add_heat(double applied_heat, double insulation_coefficient, double delta_time, double p_outside) {
            return blocked_by_stops(p_outside);

            /*

            try {
                // project where v will be; if it would overshoot a volume stop, treat it as constant volume
                double delta_h = delta_time / mass * (applied_heat * insulation_coefficient);  // time eqtn 6a

                double new_h = enthalpy + delta_h;
                new_h = ClampEnthalpy(new_h, pressure);

                double raw_v, new_v;
                bool hit_stop;
                switch (region) {
                    case ThermoMath.region_twophase:
                        double new_x = ThermoMath.x_given_ph(pressure, new_h, region);
                        //at this point, we have enough internal state to derive the rest
                        clamp_state();
                        raw_v = ThermoMath.v_given_px(pressure, new_x, region);
                        new_v = v_with_enforced_stops(raw_v, out hit_stop); // enforce volume stops
                        if (hit_stop) {
                            // hit a stop (should actually be treating as constant v); roll back to previous volume
                            // volume = new_v;
                            // calculate new delta_h
                            return true;
                        }
                        break;
                    case ThermoMath.region_liquid:
                    case ThermoMath.region_vapor:
                        // check that h is within bounds
                        // at this point, we have enough internal state to derive the rest
                        clamp_state();
                        raw_v = ThermoMath.v_given_ph(pressure, new_h, region);
                        new_v = v_with_enforced_stops(raw_v, out hit_stop); // enforce volume stops
                        if (hit_stop) {
                            // hit a stop (should actually be treating as constant v); roll back to previous volume
                            // volume = new_v;
                            return true;
                        }
                        break;
                }
                return false;
            }
            catch (Exception e) { }

            return true; // since an error was thrown when projecting volume, treating it as constant p will not work
            */
        }

        /// <summary>
        /// Probes the simulation by projecting change in pressure (uninsulated) in order to check if treating it as constant volume would be more accurate.
        /// When modifying this equation, it is recommended to finalize the ordering of operations in the add_p_uninsulated() function,
        /// then replicate that order in this function, terminating when the new volume would be set.
        /// </summary>
        /// <param name="applied_heat"></param>
        /// <param name="insulation_coefficient"></param>
        /// <param name="delta_time"></param>
        /// <returns></returns>
        private bool treat_as_constant_v_add_p_uninsulated(double new_p, double added_p, double iterative_dif, double delta_time, double insulation_coefficient) {
            try {
                //default guess
                double new_u = internalenergy;
                double new_v = volume;

                double new_s = entropy;
                double new_h = enthalpy;

                double new_x = quality;

                bool hit_stop;
                switch (region) {
                    case ThermoMath.region_liquid: //subcooled liquid
                    case ThermoMath.region_vapor: //vapor region
                                                  // Pressure Constrained -> Insulated -> delta pressure (all phases)
                                                  //default guess
                        new_u = internalenergy;
                        new_v = volume;

                        new_s = entropy;
                        new_h = enthalpy;

                        if (region == ThermoMath.region_vapor) {
                            double k = 1.27;
                            new_v = volume * Math.Pow(pressure / new_p, 1.0 / k);
                            update_vapor_vis(pressure - new_p, insulation_coefficient);
                            new_u = internalenergy - ((new_p * new_v - pressure * volume) / (1 - k));
                        }

                        double old_t = temperature;
                        double new_t = 0;
                        try {
                            new_t = ThermoMath.iterate_t_given_pv(temperature, new_p, new_v, region, true); // new_v is constant in liquid
                        }
                        catch (ArgumentOutOfRangeException e) {
                            // temperature out of range
                            ThermoMath.got_error = false;
                            return true;
                        }
                        double delta_t = new_t - old_t;
                        new_t = old_t + delta_t * (1 - insulation_coefficient); // when insulation is 0%, T = T_old

                        double delta_h = delta_time / mass * (delta_t * insulation_coefficient);  // time eqtn 6a
                        new_h = enthalpy + delta_h;

                        if (region == ThermoMath.region_vapor) {
                            double raw_v = ThermoMath.v_given_ph(new_p, new_h);
                            new_v = v_with_enforced_stops(raw_v, out hit_stop); // enforce volume stops
                            if (hit_stop) {
                                // hit a stop (should actually be treating as constant v); roll back to previous volume
                                volume = new_v;
                                return true;
                            }
                        }

                        if (region == ThermoMath.region_liquid) {
                            // new_t = ThermoMath.iterate_t_given_pv(temperature, new_p, new_v, region); // new_v is constant in liquid

                            double raw_v = ThermoMath.v_given_pt(new_p, new_t, region);
                            new_v = v_with_enforced_stops(raw_v, out hit_stop); // enforce volume stops
                            if (hit_stop) {
                                // hit a stop (should actually be treating as constant v); roll back to previous volume
                                volume = new_v;
                                return true;
                            }
                        }

                        break;
                    default:
                        break;
                }

                if (region == ThermoMath.region_twophase) //two-phase region, either newly or all along
                 {

                    double old_t = temperature;
                    double new_t = ThermoMath.tsat_given_p(new_p);
                    double delta_t = new_t - old_t;

                    double delta_s = delta_time / mass * ((delta_t * insulation_coefficient) / temperature); // time eqtn 5
                    new_s = entropy + delta_s;

                    double delta_h = delta_time / mass * (delta_t * insulation_coefficient);  // time eqtn 6a
                    new_h = enthalpy + delta_h;

                    // from this point, we have enough internal state to derive the rest

                    new_x = ThermoMath.x_given_ph(new_p, new_h);

                    int new_region = region;
                    if (new_x <= 0.0f) {
                        new_region = ThermoMath.region_liquid;

                        double raw_v = ThermoMath.vliq_given_p(new_p, new_region);
                        new_v = v_with_enforced_stops(raw_v, out hit_stop); // enforce volume stops
                        if (hit_stop) {
                            // hit a stop (should actually be treating as constant v); roll back to previous volume
                            volume = new_v;
                            return true;
                        }

                    }
                    else if (new_x > 1.0f) {
                        // TODO: this sticks to the edge of two-phase and vapor instead of transitioning

                        new_region = ThermoMath.region_vapor;
                    }
                    else {
                        double raw_v = ThermoMath.v_given_px(new_p, new_x, region);

                        new_v = v_with_enforced_stops(raw_v, out hit_stop); // enforce volume stops
                        if (hit_stop) {
                            // hit a stop (should actually be treating as constant v); roll back to previous volume
                            volume = new_v;
                            return true;
                        }
                    }
                }
            }
            catch (Exception e) { }
            return false;
        }

        /// <summary>
        /// Probes the simulation by projecting change in pressure (insulated) in order to check if treating it as constant volume would be more accurate.
        /// </summary>
        /// <param name="applied_heat"></param>
        /// <param name="insulation_coefficient"></param>
        /// <param name="delta_time"></param>
        /// <returns></returns>
        private bool treat_as_constant_v_add_p_insulated(double p, double delta_time) {
            try {
                double new_p = pressure + p; //  * delta_time;
                                             // if (Math.Abs(p) < World.DELTA_PRESSURE_CUTOFF) { new_p = pressure + p; } // small enough step; finish transition

                double new_h = enthalpy;
                double new_u = internalenergy;
                double new_t = temperature;

                switch (region) {
                    case ThermoMath.region_liquid: //subcooled liquid
                        {
                            return false;
                        }
                    case ThermoMath.region_twophase: //two-phase region
                    {
                            double insulation_coefficient = 1;

                            double old_t = temperature;
                            new_t = ThermoMath.tsat_given_p(new_p, region);

                            double delta_t = new_t - old_t;
                            double delta_h = delta_time / mass * (delta_t * insulation_coefficient);  // time eqtn 6a
                            new_h = enthalpy + delta_h;

                            double new_x = quality;
                            double new_v = volume;

                            int new_region = region;

                            // either newly or all along
                            if (region == ThermoMath.region_twophase) {
                                new_x = ThermoMath.x_given_ph(new_p, new_h, region);
                                if (new_x <= 0.0f) {
                                    new_region = ThermoMath.region_liquid;

                                    double raw_v = ThermoMath.vliq_given_p(new_p, new_region);
                                    bool hit_stop;
                                    new_v = v_with_enforced_stops(raw_v, out hit_stop); // enforce volume stops
                                    double curr_v = volume;
                                    if (hit_stop) {
                                        // hit a stop (should actually be treating as constant v); roll back to previous volume
                                        volume = new_v;
                                        // accept new region? or stay in old?
                                        return true;
                                    }
                                }
                                else if (new_x > 1.0f) {
                                    // TODO: this sticks to the edge of two-phase and vapor instead of transitioning
                                    new_region = ThermoMath.region_vapor;
                                }
                                else {
                                    double raw_v = ThermoMath.v_given_px(new_p, new_x, region);
                                    bool hit_stop;
                                    new_v = v_with_enforced_stops(raw_v, out hit_stop); // enforce volume stops
                                    double curr_v = volume;
                                    if (hit_stop) {
                                        // hit a stop (should actually be treating as constant v); roll back to previous volume
                                        volume = new_v;
                                        return true;
                                    }
                                }
                            }

                            return false;
                        }
                    case ThermoMath.region_vapor: //superheated vapor
                    {
                            // TODO: here
                            //default guess
                            new_t = temperature;
                            new_u = internalenergy;
                            double new_v = volume;

                            double k = 1.27;
                            double raw_v = volume * Math.Pow(pressure / new_p, 1.0 / k);
                            bool hit_stop;
                            new_v = v_with_enforced_stops(raw_v, out hit_stop); // enforce volume stops
                            double curr_v = volume;
                            if (hit_stop) {
                                // hit a stop (should actually be treating as constant v); roll back to previous volume
                                volume = new_v;
                                return true;
                            }

                        }
                        break;
                }
            }
            catch (Exception e) { }
            return false;
        }

        /// <summary>
        /// Probes the simulation by projecting change in pressure. If enthalpy is out of bounds, the instigating function should not run.
        /// </summary>
        /// <param name="applied_heat"></param>
        /// <param name="insulation_coefficient"></param>
        /// <param name="delta_time"></param>
        /// <returns></returns>
        private bool enthalpy_bounded(double new_p, double curr_h) {
            try {
                ThermoMath.v_given_ph_projected(new_p, curr_h);
            }
            catch {
                // enthalpy out of range
                return true;
            }
            return false;
        }

        /// <summary>
        /// Probes the simulation by projecting change in temperature. If temperature is out of bounds, the instigating function should not run.
        /// </summary>
        /// <param name="applied_heat"></param>
        /// <param name="insulation_coefficient"></param>
        /// <param name="delta_time"></param>
        /// <returns></returns>
        private bool temperature_bounded(double applied_heat, double insulation_coefficient, double delta_time) {
            try {
                // project where v will be; if it would overshoot a volume stop, treat it as constant volume
                double delta_h = delta_time / mass * (applied_heat * insulation_coefficient);  // time eqtn 6a

                double new_h = enthalpy + delta_h;
                new_h = ClampEnthalpy(new_h, pressure);

                switch (region) {
                    case ThermoMath.region_liquid:
                    case ThermoMath.region_vapor:
                        // check that h is within bounds
                        // at this point, we have enough internal state to derive the rest
                        // clamp_state();
                        double raw_v = ThermoMath.v_given_ph_projected(pressure, new_h, region);
                        break;
                }
            }
            catch (Exception e) {
                return true;
            }

            return false;
        }

        /// <summary>
        /// Enforces volume stops. Returns the enforced value.
        /// </summary>
        /// <param name="applied_heat"></param>
        /// <param name="insulation_coefficient"></param>
        /// <param name="delta_time"></param>
        /// <returns></returns>
        private double v_with_enforced_stops(double projected_v, out bool hit_stop) {
            double curr_v = volume;
            bool volume_increasing = projected_v > curr_v;
            hit_stop = false;

            for (int i = 0; i < v_stops.Count; i++) {
                VolumeStop curr_stop = v_stops[i];
                double compare_v = curr_stop.Volume;

                if (volume_increasing) {
                    // when volume would increase, enforce stops above
                    if (curr_v <= compare_v && projected_v > compare_v) {
                        projected_v = Clampd(projected_v, curr_v, compare_v - STOP_BUFFER);
                        hit_stop = true;
                        Debug.Log("[Stops] enforced, preventing increase");
                    }
                }
                else {
                    // when volume would decrease, enforce stops below
                    if (curr_v >= compare_v && projected_v < compare_v) {
                        projected_v = Clampd(projected_v, curr_v, compare_v + STOP_BUFFER);
                        hit_stop = true;
                        Debug.Log("[Stops] enforced, preventing decrease");
                    }
                }
            }

            return projected_v;
        }

        #endregion // Projections

        #region Visualization

        private void update_vapor_vis(double delta_pressure, double insulation_coefficient) {
            if (Mathf.Abs((float)delta_pressure) > World.DELTA_PRESSURE_CUTOFF) {
                Tuple<double, double> vaporInfo = new Tuple<double, double>(delta_pressure, insulation_coefficient);
                GameMgr.Events.Dispatch(GameEvents.UpdateVaporFlow, vaporInfo);
            }
        }

        #endregion // Visualization

    }
}