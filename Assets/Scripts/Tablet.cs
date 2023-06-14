using Oculus.Interaction.Input;
using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR.Controls;
using ThermoVR.State;
using ThermoVR.UI;
using UnityEngine;
using static OVRInput;

namespace ThermoVR
{

    public class Tablet : MonoBehaviour
    {
        public Touchable touchable;

        [SerializeField] private UIHub m_hub;

        [Space(5)]
        [Header("Tabs")]
        [SerializeField] private PhysicalButton m_readoutTabButton;
        [SerializeField] private PhysicalButton m_sandboxTabButton;
        [SerializeField] private PhysicalButton m_quizTabButton;
        [SerializeField] private PhysicalButton m_graphTabButton;

        private List<PhysicalButton> m_tabButtons;

        public void Init() {
            // Add buttons to list
            m_tabButtons = new List<PhysicalButton> {
                m_readoutTabButton,
                m_sandboxTabButton,
                m_quizTabButton,
                m_graphTabButton
            };

            // Initialize buttons
            for (int i = 0; i < m_tabButtons.Count; i++) {
                m_tabButtons[i].Init();
            }

            // Register button press responses
            m_readoutTabButton.OnPress += HandleReadoutTabPress;
            m_sandboxTabButton.OnPress += HandleSandboxTabPress;
            m_quizTabButton.OnPress += HandleQuizTabPress;
            m_graphTabButton.OnPress += HandleGraphTabPress;

            // Initialize UI hub, display starting screen
            m_hub.Init();
        }

        #region World Interactions

        public void SetFingerTouches(ref bool ltouch, ref bool rtouch) {
            for (int i = 0; i < m_tabButtons.Count; i++) {
                m_tabButtons[i].SetFingerTouches(ref ltouch, ref rtouch);
            }
        }

        public void CheckButtonsForPress(bool left_hand) {
            for (int i = 0; i < m_tabButtons.Count; i++) {
                m_tabButtons[i].CheckForPress(left_hand);
            }
        }

        #endregion // World Interactions

        #region Handlers

        private void HandleReadoutTabPress(object sender, EventArgs args) {
            // Open Readout UI
            m_hub.OpenUI(UIID.Readout);
        }

        private void HandleSandboxTabPress(object sender, EventArgs args) {
            // Open Sandbox UI
            m_hub.OpenUI(UIID.Sandbox);
        }

        private void HandleQuizTabPress(object sender, EventArgs args) {
            // Open Quiz UI
            m_hub.OpenUI(UIID.Quiz);
        }

        private void HandleGraphTabPress(object sender, EventArgs args) {
            // Open Graph UI
            m_hub.OpenUI(UIID.Graph);
        }

        #endregion // Handlers
    }
}
