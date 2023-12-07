using Oculus.Interaction.Input;
using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR.Controls;
using ThermoVR.State;
using ThermoVR.UI;
using TMPro;
using UnityEngine;
using static OVRInput;
using static UnityEngine.InputSystem.HID.HID;

namespace ThermoVR
{

    public class Tablet : MonoBehaviour
    {
        public Touchable touchable;

        [SerializeField] private UIHub m_hub;

        [Space(5)]
        [Header("Tabs")]
        [SerializeField] private Pressable m_sandboxTabButton;
        [SerializeField] private Pressable m_labModeButton;
        [SerializeField] private Pressable m_graphTabButton;

        [Space(5)]
        [Header("Functions")]
        [SerializeField] private Pressable m_resetButton;

        private List<Pressable> m_tabButtons;

        public void Init() {
            // Add buttons to list
            m_tabButtons = new List<Pressable> {
                m_sandboxTabButton,
                m_labModeButton,
                m_graphTabButton,
                m_resetButton
            };

            // Register button press responses
            m_sandboxTabButton.OnPress += HandleSandboxTabPress;
            m_labModeButton.OnPress += HandleQuizTabPress;
            m_graphTabButton.OnPress += HandleGraphTabPress;
            m_resetButton.OnPress += HandleResetPress;

        }

        #region World Interactions

        public void SetFingerTouches(ref bool ltouch, ref bool rtouch) {
            for (int i = 0; i < m_tabButtons.Count; i++) {
                m_tabButtons[i].SetFingerTouches(ref ltouch, ref rtouch);
            }
        }

        #endregion // World Interactions

        #region Handlers

        private void HandleSandboxTabPress(object sender, EventArgs args) {
            PlayClick(m_sandboxTabButton);

            // Open Sandbox UI
            m_hub.OpenUI(UIID.Sandbox);
        }

        private void HandleQuizTabPress(object sender, EventArgs args) {
            PlayClick(m_labModeButton);

            // Open Quiz UI
            m_hub.OpenUI(UIID.Lab);
        }

        private void HandleGraphTabPress(object sender, EventArgs args) {
            PlayClick(m_graphTabButton);

            // Open Graph UI
            m_hub.OpenUI(UIID.Graph);
        }

        private void HandleResetPress(object sender, EventArgs args)
        {
            PlayClick(m_graphTabButton);
        }

        private void PlayClick(Pressable pressable) {
            pressable.ClickAudio.Play();
        }

        #endregion // Handlers
    }
}
