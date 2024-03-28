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

    public enum StopDir
    {
        None,
        Increase,
        Decrease
    }

    public struct VolumeStop
    {
        public double Volume;
        public Tool Source;
        public StopDir StoppedDir;

        public VolumeStop(double volume, Tool source) {
            Volume = volume;
            Source = source;
            StoppedDir = StopDir.None;
        }
    }

    public abstract class Tool : MonoBehaviour, ITool
    {
        [System.NonSerialized]
        public bool engaged = false;

        [HideInInspector] public bool allowed; // whether the player can use this tool

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

        public int unique_id; // used to distinguish tools of the same type

        #endregion // Inspector

        [SerializeField] protected float m_RoutineSpeed = 1; // one control for overall routine speeds
        protected Vector3 m_ActivatedBasePos; // base position tool moves to upon activation
        protected Vector3 m_DeactivatedBasePos; // base position tool moves to upon deactivation

        protected List<GameObject> m_Elements; // individual pieces that move on dial move

        #region ITool

        protected Routine m_ActivationRoutineControl;
        protected Routine m_EngagementRoutineControl;
        protected Routine m_AdjustRoutineControl;

        protected abstract void InitializeRoutines_Impl();
        protected abstract IEnumerator ActivationRoutine(); // when tool is activated
        protected abstract IEnumerator DeactivationRoutine(); // when tool is deactivated

        protected abstract IEnumerator EngageRoutine(); // when tool dial is set to a non-zero number
        protected abstract IEnumerator DisengageRoutine(); // when tool dial is set to 0

        protected abstract IEnumerator BeginAdjustRoutine(); // when player has grabbed the dial
        protected abstract IEnumerator EndAdjustRoutine(); // when player has released the dial

        public void InitializeRoutines() {
            InitializeRoutines_Impl();
        }

        public void TriggerActivation() {
            m_ActivationRoutineControl.Replace(this, ActivationRoutine()).ExecuteWhileDisabled();
        }

        public void TriggerDeactivation() {
            m_ActivationRoutineControl.Replace(this, DeactivationRoutine()).ExecuteWhileDisabled();
        }

        public void TriggerEngage() {
            m_EngagementRoutineControl.Replace(this, EngageRoutine()).ExecuteWhileDisabled();
        }

        public void TriggerDisengage() {
            m_EngagementRoutineControl.Replace(this, DisengageRoutine()).ExecuteWhileDisabled();
        }

        public void TriggerBeginAdjust() {
            m_AdjustRoutineControl.Replace(this, BeginAdjustRoutine()).ExecuteWhileDisabled();
        }

        public void TriggerEndAdjust() {
            m_AdjustRoutineControl.Replace(this, EndAdjustRoutine()).ExecuteWhileDisabled();
        }

        #endregion // ITool

        public void Init(string unit, float mul = 1) {
            this.display_unit = unit;
            this.display_mul = mul;

            engaged = always_engaged;

            GameMgr.Events?.Register<Collider>(GameEvents.ColliderReleased, HandleColliderReleased);
            GameMgr.Events?.Register<Collider>(GameEvents.ColliderGrabbed, HandleColliderGrabbed);

            m_Elements = new List<GameObject>();

            allowed = true;

            InitializeRoutines();
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

        private void UpdateToolText(Dial dial) {
            if (dial.textv_tmpro) dial.SetValText((float)(dial.map * display_mul));
        }

        /// <summary>
        /// Moves the engaged tool position
        /// </summary>
        /// <param name="dv">the dial's value difference</param>
        public void Move(float dv) {
            foreach(var obj in m_Elements) {
                obj.transform.position += new Vector3(0, dv, 0);
            }
        }

        #region Handlers

        private void HandleColliderGrabbed(Collider col) {

        }

        private void HandleColliderReleased(Collider col) {

        }

        #endregion // Handlers
    }
}

