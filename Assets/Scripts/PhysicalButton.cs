using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Controls
{
    [RequireComponent(typeof(Touchable))]
    [RequireComponent(typeof(FingerToggleable))]
    [RequireComponent(typeof(AudioSource))]
    public class PhysicalButton : MonoBehaviour
    {
        private Touchable m_touchable;
        private FingerToggleable m_fingerToggleable;
        private AudioSource m_audioSrc;


        public event EventHandler OnPress;
        // public event EventHandler OnRelease;

        public void Init() {
            OnPress += HandlePress;

            m_touchable = GetComponent<Touchable>();
            m_fingerToggleable = GetComponent<FingerToggleable>();
            m_audioSrc = GetComponent<AudioSource>();
        }

        public void CheckForPress(bool left_hand) {
            if (m_fingerToggleable.finger) {  //finger hitting button object
                // check if the correct hand
                if ((left_hand && m_fingerToggleable.lfinger)
                    || (!left_hand && m_fingerToggleable.rfinger)) {
                    // trigger button effect
                    if (m_fingerToggleable.on) {
                        Press();
                    }
                }
            }
            else {
                m_fingerToggleable.on = false;
            }
        }

        private void Press() {
            OnPress?.Invoke(this, EventArgs.Empty);
        }

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
            m_audioSrc.Play();
            m_fingerToggleable.on = false;
        }

        #endregion // Handlers
    }
}
