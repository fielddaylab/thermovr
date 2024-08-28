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
using Hand = ThermoVR.Controls.Hand;

namespace ThermoVR
{
    [RequireComponent(typeof(AudioSource))]
    public class Tablet : MonoBehaviour
    {
        public static Tablet Instance;

        public Touchable touchable;

        [SerializeField] private UIHub m_hub;
        [SerializeField] private BoxCollider m_failsafeBounds; // prevent bug where player hand gets stuck grabbing tablet

        [Space(5)]
        [Header("Tabs")]
        [SerializeField] private Pressable m_sandboxTabButton;
        [SerializeField] private Pressable m_labModeButton;
        [SerializeField] private Pressable m_graphTabButton;
        [SerializeField] private Pressable m_gameTabButton;

        [SerializeField] private ButtonPressMovement m_sandboxTabMovement;
        [SerializeField] private ButtonPressMovement m_labModeMovement;
        [SerializeField] private ButtonPressMovement m_graphTabMovement;
        [SerializeField] private ButtonPressMovement m_gameTabMovement;

        [Space(5)]
        [Header("Functions")]
        [SerializeField] private Pressable m_resetButton;
        [SerializeField] private ButtonPressMovement m_resetButtonMovement;

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
                m_gameTabButton,
                m_resetButton
            };

            // Register button press responses
            m_sandboxTabButton.OnPress += HandleSandboxTabPress;
            m_labModeButton.OnPress += HandleQuizTabPress;
            m_graphTabButton.OnPress += HandleGraphTabPress;
            m_gameTabButton.OnPress += HandleGameTabPress;
            m_resetButton.OnPress += HandleResetPress;

            EventMgr.Events.Register(GameEvents.UISwitched, HandleUISwitched)
                .Register<string>(GameEvents.NewNameGenerated, SetUserCode, this);


            if (!ModeMgr.Instance.IsDesktop && OVRManager.display != null)
            {
                OVRManager.display.RecenteredPose += DisconnectGrab;
            }

            m_currID = UIID.Sandbox;

            // HidePullout();
        }

        public bool IsObjWithinBounds(GameObject obj)
        {
            if (obj == null)
            {
                return true;
            }

            if (m_failsafeBounds.bounds.Contains(obj.transform.position))
            {
                return true;
            }

            var point = m_failsafeBounds.ClosestPoint(obj.transform.position);
            float dist = Vector3.Distance(obj.transform.position, point);
            return dist <= 0.03f;
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

            EventMgr.Events?.Dispatch(GameEvents.SandboxModeClicked);

            // HidePullout();

            m_graphTabMovement.ResetPosition();
            m_labModeMovement.ResetPosition();
            m_sandboxTabMovement.Indent();
            m_gameTabMovement.ResetPosition();
        }

        private void HandleQuizTabPress(object sender, EventArgs args) {
            PlayClick(m_labModeButton);

            if (m_currID == UIID.Lab) { return; }
            m_currID = UIID.Lab;

            // Open Quiz UI
            m_hub.OpenUI(UIID.Lab);

            EventMgr.Events?.Dispatch(GameEvents.LabModeClicked);

            ShowPullout();

            m_graphTabMovement.ResetPosition();
            m_labModeMovement.Indent();
            m_sandboxTabMovement.ResetPosition();
            m_gameTabMovement.ResetPosition();
        }

        private void HandleGraphTabPress(object sender, EventArgs args) {
            PlayClick(m_graphTabButton);

            if (m_currID == UIID.Graph) { return; }
            m_currID = UIID.Graph;

            // Open Graph UI
            m_hub.OpenUI(UIID.Graph);

            EventMgr.Events?.Dispatch(GameEvents.SettingsViewClicked);

            ShowPullout();

            m_graphTabMovement.Indent();
            m_labModeMovement.ResetPosition();
            m_sandboxTabMovement.ResetPosition();
            m_gameTabMovement.ResetPosition();
        }

        private void HandleGameTabPress(object sender, EventArgs args)
        {
            PlayClick(m_gameTabButton);

            if (m_currID == UIID.Game) { return; }
            m_currID = UIID.Game;

            // Open Graph UI
            m_hub.OpenUI(UIID.Game);

            EventMgr.Events?.Dispatch(GameEvents.GameModeClicked);

            m_graphTabMovement.ResetPosition();
            m_labModeMovement.ResetPosition();
            m_sandboxTabMovement.ResetPosition();
            m_gameTabMovement.Indent();
        }

        private void HandleResetPress(object sender, EventArgs args)
        {
            PlayClick(m_graphTabButton);

            EventMgr.Events.Dispatch(GameEvents.ResetPressed);
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
            Hand handType = wasLeft ? Hand.LEFT : Hand.RIGHT;
            if (wasAny) { touchable.SetGrabbed(false, handType); }
        }

        private void PlayClick(Pressable pressable) {
            if (GameMgr.I.AudioEnabled) { pressable.ClickAudio.Play(); }

        }

        private void HandleGrabbed(object sender, Hand arg)
        {
            PositionDataFrame currPos = new PositionDataFrame();
            currPos.pos = new float[] { this.transform.position.x, this.transform.position.y, this.transform.position.z };
            currPos.rot = new float[] { this.transform.rotation.x, this.transform.rotation.y, this.transform.rotation.z, this.transform.rotation.w };

            EventMgr.Events?.Dispatch(GameEvents.TabletGrabbed, new Tuple<PositionDataFrame, Hand>(currPos, arg));
        }

        private void HandleReleased(object sender, Hand arg)
        {
            PositionDataFrame currPos = new PositionDataFrame();
            currPos.pos = new float[] { this.transform.position.x, this.transform.position.y, this.transform.position.z };
            currPos.rot = new float[] { this.transform.rotation.x, this.transform.rotation.y, this.transform.rotation.z, this.transform.rotation.w };

            EventMgr.Events?.Dispatch(GameEvents.TabletReleased, new Tuple<PositionDataFrame, Hand>(currPos, arg));
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
