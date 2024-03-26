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
    [RequireComponent(typeof(AudioSource))]
    public class Tablet : MonoBehaviour
    {
        public static Tablet Instance;

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

        [Space(5)]
        [Header("Analytics")]
        [SerializeField] private TMP_Text[] m_PlayerCodeDisplays;

        private List<Pressable> m_tabButtons;

        private UIID m_currID;

        private AudioSource m_audioSource;

        public void Init() {
            if (Instance == null)
            {
                Instance = this;
            }
            else if (this != Instance)
            {
                return;
            }

            touchable.OnGrab += HandleGrabbed;
            touchable.OnRelease += HandleReleased;

            m_audioSource = this.GetComponent<AudioSource>();

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

            GameMgr.Events.Register(GameEvents.UISwitched, HandleUISwitched)
                .Register<string>(GameEvents.NewNameGenerated, SetUserCode, this);


            if (!GameMgr.I.IsDesktop && OVRManager.display != null)
            {
                OVRManager.display.RecenteredPose += DisconnectGrab;
            }

            m_currID = UIID.Sandbox;

            HidePullout();
        }

        #region World Interactions

        // Constant source for UI elements that often show/hide on click
        public void PlayUIAudio(AudioClip clip)
        {
            m_audioSource.PlayOneShot(clip);
        }

        #endregion // World Interactions

        #region Handlers

        private void HandleSandboxTabPress(object sender, EventArgs args) {
            PlayClick(m_sandboxTabButton);

            if (m_currID == UIID.Sandbox) { return; }
            m_currID = UIID.Sandbox;

            // Open Sandbox UI
            m_hub.OpenUI(UIID.Sandbox);

            GameMgr.Events?.Dispatch(GameEvents.SandboxModeClicked);

            HidePullout();
        }

        private void HandleQuizTabPress(object sender, EventArgs args) {
            PlayClick(m_labModeButton);

            if (m_currID == UIID.Lab) { return; }
            m_currID = UIID.Lab;

            // Open Quiz UI
            m_hub.OpenUI(UIID.Lab);

            GameMgr.Events?.Dispatch(GameEvents.LabModeClicked);

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

        private void HandleUISwitched()
        {
            DisconnectGrab();
        }

        private void SetUserCode(string userCode)
        {
            for (int i = 0; i < m_PlayerCodeDisplays.Length; i++)
            {
                if (m_PlayerCodeDisplays[i])
                {
                    m_PlayerCodeDisplays[i].SetText(userCode);
                }
            }
        }

        private void DisconnectGrab()
        {
            bool wasAny = touchable.ltouch || touchable.rtouch;
            bool wasLeft = touchable.ltouch;
            touchable.rtouch = false;
            touchable.ltouch = false;
            touchable.touch = false;
            if (wasAny) { touchable.SetGrabbed(false, wasLeft); }
        }

        private void PlayClick(Pressable pressable) {
            if (GameMgr.I.AudioEnabled) { pressable.ClickAudio.Play(); }

        }

        private void HandleGrabbed(object sender, bool arg)
        {
            GameMgr.Events?.Dispatch(GameEvents.TabletGrabbed, this.transform);
        }

        private void HandleReleased(object sender, bool arg)
        {
            GameMgr.Events?.Dispatch(GameEvents.TabletReleased, this.transform);
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
