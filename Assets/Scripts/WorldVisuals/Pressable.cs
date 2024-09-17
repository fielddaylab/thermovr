using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using ThermoVR.Controls;
using BeauUtil.Extensions;

namespace ThermoVR
{
    [RequireComponent(typeof(Touchable))]
    [RequireComponent(typeof(FingerToggleable))]
    [RequireComponent(typeof(AudioSource))]
    public class Pressable : MonoBehaviour
    {
        private Touchable m_touchable;
        private FingerToggleable m_fingerToggleable;
        public AudioSource ClickAudio;

        [SerializeField] private double m_touchTime = 0.1; // min time before new press is registered
        private double m_touchTimer; // time remaining before new press may be registered

        public event EventHandler OnPress;
        public event EventHandler PressCompleted;
        // public event EventHandler OnRelease;

        #region Unity Callbacks

        private void Awake() {
            m_touchable = GetComponent<Touchable>();
            m_fingerToggleable = GetComponent<FingerToggleable>();
            ClickAudio = GetComponent<AudioSource>();

            m_touchTimer = 0;
        }

        private void Start() {
            EventMgr.Events?.Register(GameEvents.GatherPressables, HandleGatherPressables);
            EventMgr.Events?.Register<bool>(GameEvents.CheckForPress, HandleCheckForPress);
        }

        private void FixedUpdate() {
            if (m_touchTimer > 0) {
                m_touchTimer -= Time.fixedDeltaTime;
            }
        }

        #endregion // Unity Callbacks

        public void OnEnable() {
            OnPress += HandlePress;
            PressCompleted += HandlePressCompleted;
        }

        public void OnDisable() {
            OnPress += HandlePress;
            PressCompleted += HandlePressCompleted;
        }

        /// <summary>
        /// Triggers the button press
        /// </summary>
        /// <param name="cooldown">Cooldown if in VR, none if in desktop</param>
        public void Press(bool cooldown, Hand inputType) {
            if (m_touchTimer <= 0) {
                EventMgr.Events.Dispatch(GameEvents.HandStartPress, EvtArgs.Create(inputType));
                OnPress?.Invoke(this, EventArgs.Empty);
                PressCompleted?.Invoke(this, EventArgs.Empty);
                if (cooldown) {
                    m_touchTimer = m_touchTime;
                }
            }
        }

        /*
        public void PressByProxy() {
            Press();
        }
        */

        #region Queries

        /// <summary>
        /// 
        /// </summary>
        /// <param name="queryingLeft">true if querying the left finger, false if the right finger</param>
        /// <returns></returns>
        public void SetFingerTouches(ref bool ltouch, ref bool rtouch) {
            m_touchable.SetFingerTouches(ref ltouch, ref rtouch);
        }

        #endregion // Queries

        #region Handlers

        private void HandlePress(object sender, EventArgs args) {
            // m_audioSrc.Play();
            m_fingerToggleable.on = false;
        }

        private void HandlePressCompleted(object sender, EventArgs args) {
            // nothing yet
        }

        private void HandleGatherPressables() {
            EventMgr.Events.Dispatch(GameEvents.RegisterPressable, EvtArgs.Ref(this));
        }

        private void HandleCheckForPress(bool left_hand) {
            if (m_fingerToggleable.finger) {  //finger hitting pressable object
                // check if the correct hand
                if ((left_hand && m_fingerToggleable.lfinger)
                    || (!left_hand && m_fingerToggleable.rfinger)) {
                    // trigger button effect
                    bool isLeft = left_hand && m_fingerToggleable.lfinger;
                    if (m_fingerToggleable.on) {
                        Hand pressType = isLeft ? Hand.LEFT : Hand.RIGHT;
                        Press(true, pressType);
                        if (isLeft) { m_fingerToggleable.lfinger = false; }
                        else { m_fingerToggleable.rfinger = false; }
                    }
                }
            }
            else {
                m_fingerToggleable.on = false;
            }
        }

        #endregion // Handlers
    }
}
