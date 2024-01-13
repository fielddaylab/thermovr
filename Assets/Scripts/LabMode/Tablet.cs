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

        [Space(5)]
        [Header("Pullout Readout")]
        [SerializeField] private GameObject m_pulloutReadoutModel;
        [SerializeField] private GameObject m_pulloutReadoutScreen;

        private List<Pressable> m_tabButtons;

        private UIID m_currID;

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

            m_currID = UIID.Sandbox;

            HidePullout();
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

            if (m_currID == UIID.Sandbox) { return; }
            m_currID = UIID.Sandbox;

            // Open Sandbox UI
            m_hub.OpenUI(UIID.Sandbox);

            HidePullout();
        }

        private void HandleQuizTabPress(object sender, EventArgs args) {
            PlayClick(m_labModeButton);

            if (m_currID == UIID.Lab) { return; }
            m_currID = UIID.Lab;

            // Open Quiz UI
            m_hub.OpenUI(UIID.Lab);

            ShowPullout();
        }

        private void HandleGraphTabPress(object sender, EventArgs args) {
            PlayClick(m_graphTabButton);

            if (m_currID == UIID.Graph) { return; }
            m_currID = UIID.Graph;

            // Open Graph UI
            m_hub.OpenUI(UIID.Graph);

            ShowPullout();
        }

        private void HandleResetPress(object sender, EventArgs args)
        {
            PlayClick(m_graphTabButton);
        }

        private void PlayClick(Pressable pressable) {
            pressable.ClickAudio.Play();
        }

        #endregion // Handlers

        private void ShowPullout()
        {
            m_pulloutReadoutModel.SetActive(true);
            m_pulloutReadoutScreen.SetActive(true);
        }

        private void HidePullout()
        {
            m_pulloutReadoutModel.SetActive(false);
            m_pulloutReadoutScreen.SetActive(false);
        }
    }
}
