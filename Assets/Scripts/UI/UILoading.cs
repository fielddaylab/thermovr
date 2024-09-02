using System.Collections;
using System.Collections.Generic;
using ThermoVR.Lab;
using ThermoVR.UI.GraphElements;
using UnityEngine;
using UnityEngine.UI;
using static ThermoVR.Analytics.AnalyticsService;

namespace ThermoVR.UI
{
    [DefaultExecutionOrder(5)]
    public class UILoading : MonoBehaviour
    {
        [SerializeField] private Canvas m_canvas;
        [SerializeField] private Image m_loadIcon;
        [SerializeField] private Vector3 m_rotation;

        [SerializeField] private float m_minTimer = 3;

        [SerializeField] private CameraClampToVirtualViewport m_clampCam;

        [Header("Elements")]
        [SerializeField] private RectTransform m_wrapperGroup;
        [SerializeField] private GameObject m_desktopGroup;
        [SerializeField] private GameObject m_vrGroup;
        [SerializeField] private RectTransform m_loadGroup;
        [SerializeField] private Button m_startButton;
        [SerializeField] private UICredits m_credits;
        [SerializeField] private Button m_creditsButton;

        [Header("Settings")]
        [SerializeField] private Button m_musicOption;
        [SerializeField] private Image m_musicFill;
        private bool m_musicSelected;
        [SerializeField] private Button m_fullscreenOption;
        [SerializeField] private Image m_fullscreenFill;
        private bool m_fullscreenSelected;

        [Header("VR layout")]
        [SerializeField] private Vector3 m_vrLoadPos;
        [SerializeField] private Vector3 m_vrLoadScale;
        [SerializeField] private Vector3 m_vrWrapperPos;
        [SerializeField] private Vector3 m_vrWrapperScale;

        private void Awake()
        {
            if (ModeMgr.Instance.IsDesktop)
            {
                m_canvas.renderMode = RenderMode.ScreenSpaceOverlay;

                m_desktopGroup.SetActive(true);
                m_vrGroup.SetActive(false);

                EventMgr.Events.Register(GameEvents.AsyncLoadComplete, HandleAsyncLoadComplete);
                m_startButton.onClick.AddListener(HandleStartGameClicked);
                m_startButton.interactable = false;

                m_musicSelected = true;
                m_musicOption.onClick.AddListener(HandleMusicOptionSelected);
                m_musicFill.enabled = m_musicSelected;

                m_fullscreenOption.onClick.AddListener(HandleFullscreenOptionSelected);
                m_fullscreenFill.enabled = m_fullscreenSelected;

                m_creditsButton.onClick.AddListener(HandleTitleCreditsClicked);
            }
            else
            {
                m_canvas.renderMode = RenderMode.ScreenSpaceCamera;
                m_clampCam.enabled = false;

                m_desktopGroup.SetActive(false);
                m_vrGroup.SetActive(true);

                m_loadGroup.anchoredPosition = m_vrLoadPos;
                m_loadGroup.localScale = m_vrLoadScale;
                m_wrapperGroup.anchoredPosition = m_vrWrapperPos;
                m_wrapperGroup.localScale = m_vrWrapperScale;
            }

            m_credits.Init(true);
            m_credits.gameObject.SetActive(false);

            EventMgr.Events.Register(GameEvents.TitleCreditsClosed, HandleTitleCreditsClosed);

            EventMgr.Events.Dispatch(GameEvents.TitleScreenDisplayed);
        }

        private void Update()
        {
            if (m_minTimer > 0)
            {
                m_minTimer -= Time.deltaTime;
            }

            if (ModeMgr.Instance.IsDesktop && m_startButton.interactable)
            {
                m_loadGroup.gameObject.SetActive(false);
            }
            else
            {
                m_loadIcon.transform.Rotate(-m_rotation * Time.deltaTime, Space.Self);
            }
        }

        public bool MinLoadTimeCompleted()
        {
            return m_minTimer <= 0;
        }
        
        public bool GetFullscreen()
        {
            return m_fullscreenSelected;
        }

        public bool GetMusic()
        {
            return m_musicSelected;
        }

        #region Handlers

        private void HandleAsyncLoadComplete()
        {
            m_startButton.interactable = true;
        }

        private void HandleStartGameClicked()
        {
            EventMgr.Events.Dispatch(GameEvents.StartGameClicked);

            EventMgr.Events.Dispatch(GameEvents.ClickCloseTitleScreen);

            EventMgr.Events.Dispatch(GameEvents.TitleScreenClosed);
        }

        private void HandleMusicOptionSelected()
        {
            m_musicSelected = !m_musicSelected;
            m_musicFill.enabled = m_musicSelected;

            EventMgr.Events.Dispatch(GameEvents.ClickToggleTitleSetting, new GraphSettingUpdate(GraphElementID.Music, m_fullscreenSelected));
        }

        private void HandleFullscreenOptionSelected()
        {
            m_fullscreenSelected = !m_fullscreenSelected;
            m_fullscreenFill.enabled = m_fullscreenSelected;

            EventMgr.Events.Dispatch(GameEvents.ClickToggleTitleSetting, new GraphSettingUpdate(GraphElementID.Fullscreen, m_fullscreenSelected));
        }

        private void HandleTitleCreditsClosed()
        {
            m_wrapperGroup.gameObject.SetActive(true);

            EventMgr.Events.Dispatch(GameEvents.ClickCloseCredits, true);
        }

        private void HandleTitleCreditsClicked()
        {
            m_wrapperGroup.gameObject.SetActive(false);
            EventMgr.Events.Dispatch(GameEvents.TitleCreditsOpened);

            EventMgr.Events.Dispatch(GameEvents.ClickDisplayCredits, true);

        }

        #endregion // Handlers
    }
}