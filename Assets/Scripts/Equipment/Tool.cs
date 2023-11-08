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
using BeauRoutine;
using UnityEngine.Events;

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

    public abstract class Tool : MonoBehaviour, ITool
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

        [HideInInspector] public UnityEvent ValUpdated;

        #endregion // Inspector

        #region ITool

        protected Routine m_TransitionRoutine;

        protected abstract IEnumerator ActivationRoutine();
        protected abstract IEnumerator DeactivationRoutine();

        public void TriggerActivation() {
            m_TransitionRoutine.Replace(this, ActivationRoutine()).ExecuteWhileDisabled();
        }

        public void TriggerDeactivation() {
            m_TransitionRoutine.Replace(this, DeactivationRoutine()).ExecuteWhileDisabled();
        }

        #endregion // ITool

        public void Init(string unit, float mul = 1) {
            this.display_unit = unit;
            this.display_mul = mul;

            engaged = always_engaged;
        }

        public void UpdateVal(float new_val, Dial dial) {
            val = new_val;

            ValUpdated?.Invoke();

            UpdateToolText(dial);
        }

        public void UpdateText(Dial dial) {
            UpdateToolText(dial);
        }

        public float GetVal() {
            return val;
        }

        /// <summary>
        /// Moves the engaged tool position
        /// </summary>
        /// <param name="dv">the dial's value difference</param>
        public void Move(float dv) {
            this.transform.position += new Vector3(0, dv, 0);
        }


        private void UpdateToolText(Dial dial) {
            if (dial.textv_tmpro) dial.SetValText((float)(dial.map * display_mul));
        }
    }
}

