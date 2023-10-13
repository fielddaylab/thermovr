using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR
{
    public static class MathUtility
    {
        //safe to call if not interactable, as it will just do nothing
        public static bool floatNumeric(float f) {
            if (double.IsNaN(f)) return false;
            if (double.IsInfinity(f)) return false;
            return true;
        }

        public static double Clampd(double v, double min, double max) {
            if (v < min) return min;
            if (v > max) return max;
            return v; 
        }
        public static float Clampf(float v, float min, float max) { 
            if (v < min) return min;
            if (v > max) return max;
            return v; 
        }

        public static Vector3 popVector() {
            return new Vector3(UnityEngine.Random.Range(-1f, 1f), 1f, UnityEngine.Random.Range(-1f, 1f));
        }
    }
}
