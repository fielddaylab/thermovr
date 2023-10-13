/*
DOCUMENTATION- phil, 12/16/19 [intended to be a description as of a point in time, NOT nec a prescription for how it should be evolved- feel free to uproot]

The various tools which can be variably engaged to the container, thrown, or stored.
*/

using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using TMPro;
using ThermoVR.Dials;
using System;

namespace ThermoVR.Tools
{
    [Serializable]
    public enum ToolType { 
        Burner,
        Coil,
        Weight,
        NegativeWeight,
        Insulator,
        SurroundingTemperature,
        SurroundingPressure,
        Stops
    }

    public struct VolumeStop
    {
        public double Volume;
        public Tool Source;

        public VolumeStop(double volume, Tool source) {
            Volume = volume;
            Source = source;
        }
    }

    public class Tool : MonoBehaviour
    {
        [System.NonSerialized]
        public bool engaged = false;

        private float val = 0.0f; // value of the tool

        [System.NonSerialized]
        public string display_unit = "";
        [System.NonSerialized]
        public float display_mul = 1.0f; //multiplied with map before displaying with display_unit

        #region Inspector

        public bool always_engaged = false;
        public ToolType tool_type;
        [SerializeField] private GameObject ModelContainer;

        #endregion // Inspector

        public void Init(string unit, float mul = 1) {
            this.display_unit = unit;
            this.display_mul = mul;

            engaged = always_engaged;
        }

        public void update_val(float new_val, Dial dial) {
            val = new_val;

            update_tool_text(dial);
        }

        public void update_text(Dial dial) {
            update_tool_text(dial);
        }

        public float get_val() {
            return val;
        }

        /// <summary>
        /// Moves the engaged tool position
        /// </summary>
        /// <param name="dv">the dial's value difference</param>
        public void move(float dv) {
            this.transform.position += new Vector3(0, dv, 0);
        }

        private void update_tool_text(Dial dial) {
            if (dial.textv_tmpro) dial.SetValText(display_unit, (float)(dial.map * display_mul));
        }
    }
}

