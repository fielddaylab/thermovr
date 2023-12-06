using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

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
            GameMgr.Events?.Register(GameEvents.GatherPressables, HandleGatherPressables);
            GameMgr.Events?.Register<bool>(GameEvents.CheckForPress, HandleCheckForPress);

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
        public void Press(bool cooldown) {
            if (m_touchTimer <= 0) {
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
            GameMgr.Events.Dispatch(GameEvents.RegisterPressable, this);
        }

        private void HandleCheckForPress(bool left_hand) {
            if (m_fingerToggleable.finger) {  //finger hitting pressable object
                // check if the correct hand
                if ((left_hand && m_fingerToggleable.lfinger)
                    || (!left_hand && m_fingerToggleable.rfinger)) {
                    // trigger button effect
                    if (m_fingerToggleable.on) {
                        Press(true);
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
