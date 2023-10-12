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
    }
}
