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

        public GameObject textl;
        [System.NonSerialized]
        public TextMeshPro textl_tmpro;
        [System.NonSerialized]
        public MeshRenderer textl_meshrenderer;

        public GameObject textd;
        [System.NonSerialized]
        public TextMeshPro textd_tmpro;
        [System.NonSerialized]
        public MeshRenderer textd_meshrenderer;

        public GameObject textn;

        [System.NonSerialized]
        public bool disabled;

        [System.NonSerialized]
        public string display_unit = "";
        [System.NonSerialized]
        public float display_mul = 1.0f; //multiplied with map before displaying with display_unit

        [System.NonSerialized]
        public AudioSource audioS;

        public bool always_engaged = false;

        public ToolType tool_type;

        private bool examined;

        public void Init(string unit, float mul = 1) {
            this.display_unit = unit;
            this.display_mul = mul;

            engaged = always_engaged;

            Setup();
        }

        private void Setup() {
            // TODO: separate textv from dial, move to tool
            // textv = dial_dial.textv.gameObject;
            //textv_tmpro = textv.GetComponent<TextMeshPro>();
            //textv_meshrenderer = textv.GetComponent<MeshRenderer>();
            textl_tmpro = textl.GetComponent<TextMeshPro>();
            textl_meshrenderer = textl.GetComponent<MeshRenderer>();
            if (textd) {
                textd_tmpro = textd.GetComponent<TextMeshPro>();
                textd_meshrenderer = textd.GetComponent<MeshRenderer>();
            }

            audioS = gameObject.GetComponent<AudioSource>();

            disabled = false;
            enabled = true;
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
            // find dial
            // if (!dial_dial) return;
            // do not apply to insulator
            // if (t == tool_insulator) return;

            string txt = string.Format("{0:3} " + display_unit, (float)(dial.map * display_mul));
            //t.textv_tmpro.SetText(txt);
            //if(t.textd_tmpro) t.textd_tmpro.SetText(txt);
            if (dial.textv_tmpro) dial.SetValText(display_unit, (float)(dial.map * display_mul));
            if (textd_tmpro) textd_tmpro.SetText("{0:3} " + display_unit, (float)(dial.map * display_mul));
        }
    }
}

