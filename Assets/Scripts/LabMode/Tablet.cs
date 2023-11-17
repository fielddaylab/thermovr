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
        [SerializeField] private Pressable m_readoutTabButton;
        [SerializeField] private Pressable m_sandboxTabButton;
        [SerializeField] private Pressable m_quizTabButton;
        [SerializeField] private Pressable m_graphTabButton;

        [Header("Pullout")]
        [SerializeField] private Pressable m_expandToggleButton;
        [SerializeField] private GameObject m_expandedModel;
        [SerializeField] private GameObject m_collapsedModel;
        [SerializeField] private TextMeshPro m_toggleText;


        private bool m_isExpanded;

        private List<Pressable> m_tabButtons;

        public void Init() {
            // Add buttons to list
            m_tabButtons = new List<Pressable> {
                m_readoutTabButton,
                m_sandboxTabButton,
                m_quizTabButton,
                m_graphTabButton
            };

            // Initialize buttons
            for (int i = 0; i < m_tabButtons.Count; i++) {
                // m_tabButtons[i].Init();
            }

            // Register button press responses
            m_readoutTabButton.OnPress += HandleReadoutTabPress;
            m_sandboxTabButton.OnPress += HandleSandboxTabPress;
            m_quizTabButton.OnPress += HandleQuizTabPress;
            m_graphTabButton.OnPress += HandleGraphTabPress;

            m_expandToggleButton.OnPress += HandleExpandToggleButtonPress;

            m_isExpanded = false;

            ToggleExpandedModel();
        }

        #region World Interactions

        public void SetFingerTouches(ref bool ltouch, ref bool rtouch) {
            for (int i = 0; i < m_tabButtons.Count; i++) {
                m_tabButtons[i].SetFingerTouches(ref ltouch, ref rtouch);
            }
        }

        /*
        public void CheckButtonsForPress(bool left_hand) {
            for (int i = 0; i < m_tabButtons.Count; i++) {
                m_tabButtons[i].CheckForPress(left_hand);
            }
        }
        */

        #endregion // World Interactions

        #region Handlers

        private void HandleReadoutTabPress(object sender, EventArgs args) {
            PlayClick(m_readoutTabButton);

            // Open Readout UI
            m_hub.OpenUI(UIID.Readout);
        }

        private void HandleSandboxTabPress(object sender, EventArgs args) {
            PlayClick(m_sandboxTabButton);

            // Open Sandbox UI
            m_hub.OpenUI(UIID.Sandbox);
        }

        private void HandleQuizTabPress(object sender, EventArgs args) {
            PlayClick(m_quizTabButton);

            // Open Quiz UI
            m_hub.OpenUI(UIID.Quiz);
        }

        private void HandleGraphTabPress(object sender, EventArgs args) {
            PlayClick(m_graphTabButton);

            // Open Graph UI
            m_hub.OpenUI(UIID.Graph);
        }

        private void HandleExpandToggleButtonPress(object sender, EventArgs args) {
            PlayClick(m_expandToggleButton);

            ToggleExpandedModel();
        }

        private void PlayClick(Pressable pressable) {
            pressable.ClickAudio.Play();
        }

        #endregion // Handlers

        private void ToggleExpandedModel() {
            // expand / collapse readouts
            m_isExpanded = !m_isExpanded;

            m_expandedModel.SetActive(m_isExpanded);
            m_collapsedModel.SetActive(!m_isExpanded);

            if (m_isExpanded) {
                m_toggleText.SetText(">");
            }
            else
            {
                m_toggleText.SetText("<");
            }
        }
    }
}
