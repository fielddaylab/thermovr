using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR;
using ThermoVR.Controls;
using ThermoVR.UI;
using ThermoVR.UI.GraphElements;
using ThermoVR.UI.Interfaces;
using UnityEngine;
using UnityEngine.UI;

namespace ThermoVR
{
    public enum TabletMode
    {
        NONE,
        SANDBOX,
        LAB,
        GAME,
        SETTINGS,
    }

    public class GraphModule : UIModule
    {
        #region Inspector

        [Header("Groups")]
        [SerializeField] private CanvasGroup m_settingsGroup;
        [SerializeField] private CanvasGroup m_configureGroup;
        [SerializeField] private CanvasGroup m_controlsGroup;
        [SerializeField] private CanvasGroup m_creditsGroup;

        [Header("Settings")]
        [SerializeField] private ThermoButton m_creditsButton;
        [SerializeField] private ThermoButton m_configureButton;
        [SerializeField] private ThermoButton m_controlsButton;
        [SerializeField] private Image m_configureGraphic;
        [SerializeField] private Image m_controlsGraphic;
        [SerializeField] private RectTransform m_configureTextRect;
        [SerializeField] private RectTransform m_controlsTextRect;

        [SerializeField] private float m_selectedTabY;
        [SerializeField] private float m_unselectedTabY;
        [SerializeField] private Sprite m_selectedTabSprite;
        [SerializeField] private Sprite m_unselectedTabSprite;
        [SerializeField] private float m_selectedTextRectOffset = 3.11f;

        [Header("Configure")]
        [SerializeField] private ThermoToggle m_axisNumbersToggle;
        [SerializeField] private ThermoToggle m_gridLinesToggle;
        [SerializeField] private ThermoToggle m_regionLabelsToggle;
        [SerializeField] private ThermoToggle m_constantLinesToggle;
        [SerializeField] private ThermoToggle m_axisTrackersToggle;
        [SerializeField] private ThermoToggle m_musicToggle;

        [Header("Credits")]
        [SerializeField] private ThermoButton m_closeCreditsButton;
        [SerializeField] private UICredits m_credits;

        private ThermoToggle[] m_toggles;

        public override void Init()
        {
            base.Init();

            m_toggles = new ThermoToggle[6] {
            m_axisNumbersToggle,
            m_gridLinesToggle,
            m_regionLabelsToggle,
            m_constantLinesToggle,
            m_axisTrackersToggle,
            m_musicToggle
        };

            if (PersistentState.Instance.Bools.ContainsKey(PersistentVars.MusicOnStart))
            {
                m_musicToggle.SetToggleState(PersistentState.Instance.Bools[PersistentVars.MusicOnStart]);
            }

            for (int i = 0; i < m_toggles.Length; i++)
            {
                m_toggles[i].Init();
            }

            this.gameObject.SetActive(true);
            m_credits.Init(false);
            this.gameObject.SetActive(false);

            if (ModeMgr.Instance.IsDesktop)
            {
                // Hide Controls Button
                m_controlsButton.gameObject.SetActive(false);

                // center the only remaining button
                m_configureButton.Rect.anchoredPosition = new Vector2(
                    0,
                    m_configureButton.Rect.anchoredPosition.y);
            }

            Array settings = Enum.GetValues(typeof(GraphElementID));
            for (int i = 0; i < settings.Length; i++)
            {
                DispatchSettingUpdate((GraphElementID)settings.GetValue(i));
            }
        }
        #endregion // Inspector

        #region IUIModule

        public override void Open()
        {
            this.gameObject.SetActive(true);

            SetMainPanelVisible(true);
            SetCreditsVisible(false);
            SelectConfigureTab();

            AddListeners();
        }

        public override void Close()
        {
            this.gameObject.SetActive(false);

            RemoveListeners();
        }

        #endregion // IUIModule

        #region Helpers

        private void AddListeners()
        {
            m_axisNumbersToggle.Pressable.PressCompleted += HandleAxisNumbersToggle;
            m_gridLinesToggle.Pressable.PressCompleted += HandleGridLinesToggle;
            m_regionLabelsToggle.Pressable.PressCompleted += HandleRegionLabelsToggle;
            m_constantLinesToggle.Pressable.PressCompleted += HandleConstantLinesToggle;
            m_axisTrackersToggle.Pressable.PressCompleted += HandleAxisTrackersToggle;
            m_musicToggle.Pressable.PressCompleted += HandleMusicToggle;

            m_creditsButton.OnButtonPressed += HandleCreditsButtonPressed;
            m_closeCreditsButton.OnButtonPressed += HandleCloseCreditsButtonPressed;
            m_configureButton.OnButtonPressed += HandleConfigureButtonPressed;
            m_controlsButton.OnButtonPressed += HandleControlsButtonPressed;
        }

        private void RemoveListeners()
        {
            m_axisNumbersToggle.Pressable.PressCompleted -= HandleAxisNumbersToggle;
            m_gridLinesToggle.Pressable.PressCompleted -= HandleGridLinesToggle;
            m_regionLabelsToggle.Pressable.PressCompleted -= HandleRegionLabelsToggle;
            m_constantLinesToggle.Pressable.PressCompleted -= HandleConstantLinesToggle;
            m_axisTrackersToggle.Pressable.PressCompleted -= HandleAxisTrackersToggle;
            m_musicToggle.Pressable.PressCompleted -= HandleMusicToggle;

            m_creditsButton.OnButtonPressed -= HandleCreditsButtonPressed;
            m_closeCreditsButton.OnButtonPressed -= HandleCloseCreditsButtonPressed;
            m_configureButton.OnButtonPressed -= HandleConfigureButtonPressed;
            m_controlsButton.OnButtonPressed -= HandleControlsButtonPressed;
        }

        private void DispatchSettingUpdate(GraphElementID id)
        {
            ThermoToggle toggle;

            switch (id)
            {
                case GraphElementID.AxisNumbers:
                    toggle = m_axisNumbersToggle;
                    break;
                case GraphElementID.RegionLabels:
                    toggle = m_regionLabelsToggle;
                    break;
                case GraphElementID.GridLines:
                    toggle = m_gridLinesToggle;
                    break;
                case GraphElementID.ConstantLines:
                    toggle = m_constantLinesToggle;
                    break;
                case GraphElementID.AxisTrackers:
                    toggle = m_axisTrackersToggle;
                    break;
                case GraphElementID.Music:
                    toggle = m_musicToggle;
                    break;
                default:
                    return;
            }

            bool val = toggle.IsOn();
            EventMgr.Events.Dispatch(GameEvents.UpdateGraphSetting, new GraphSettingUpdate(id, val));
        }

        private void SetCreditsVisible(bool isVisible)
        {
            m_creditsGroup.alpha = isVisible ? 1 : 0;
            m_creditsGroup.interactable = isVisible;
            m_creditsGroup.blocksRaycasts = isVisible;
        }

        private void SetMainPanelVisible(bool isVisible)
        {
            m_settingsGroup.alpha = isVisible ? 1 : 0;
            m_settingsGroup.interactable = isVisible;
            m_settingsGroup.blocksRaycasts = isVisible;
        }

        private void SelectConfigureTab()
        {
            EventMgr.Events.Dispatch(GameEvents.ClickConfigTab);

            m_configureGroup.alpha = 1;
            m_configureGroup.interactable = true;
            m_configureGroup.blocksRaycasts = true;

            m_controlsGroup.alpha = 0;
            m_controlsGroup.interactable = false;
            m_controlsGroup.blocksRaycasts = false;

            m_configureGraphic.sprite = m_selectedTabSprite;
            m_configureGraphic.SetNativeSize();
            m_controlsGraphic.sprite = m_unselectedTabSprite;
            m_controlsGraphic.SetNativeSize();

            m_configureTextRect.offsetMin = new Vector2(m_configureTextRect.offsetMin.x, m_selectedTextRectOffset);
            m_controlsTextRect.offsetMin = new Vector2(m_controlsTextRect.offsetMin.x, 0);

            m_configureButton.Rect.anchoredPosition = new Vector2(
                m_configureButton.Rect.anchoredPosition.x,
                m_selectedTabY
                );

            m_controlsButton.Rect.anchoredPosition = new Vector2(
                 m_controlsButton.Rect.anchoredPosition.x,
                 m_unselectedTabY
                 );


        }

        private void SelectControlsTab()
        {
            EventMgr.Events.Dispatch(GameEvents.ClickControlsTab);

            m_configureGroup.alpha = 0;
            m_configureGroup.interactable = false;
            m_configureGroup.blocksRaycasts = false;

            m_controlsGroup.alpha = 1;
            m_controlsGroup.interactable = true;
            m_controlsGroup.blocksRaycasts = true;

            m_configureGraphic.sprite = m_unselectedTabSprite;
            m_configureGraphic.SetNativeSize();
            m_controlsGraphic.sprite = m_selectedTabSprite;
            m_controlsGraphic.SetNativeSize();

            m_configureTextRect.offsetMin = new Vector2(m_configureTextRect.offsetMin.x, 0);
            m_controlsTextRect.offsetMin = new Vector2(m_controlsTextRect.offsetMin.x, m_selectedTextRectOffset);

            m_configureButton.Rect.anchoredPosition = new Vector2(
                m_configureButton.Rect.anchoredPosition.x,
                m_unselectedTabY
                );

            m_controlsButton.Rect.anchoredPosition = new Vector2(
                m_controlsButton.Rect.anchoredPosition.x,
                m_selectedTabY
                );
        }

        #endregion // Helpers

        #region Handlers

        // Also see if you can disconnect toggles from general touch checking; otherwise incorporate into touch pipeline

        private void HandleAxisNumbersToggle(object sender, EventArgs args)
        {
            DispatchSettingUpdate(GraphElementID.AxisNumbers);
        }

        private void HandleGridLinesToggle(object sender, EventArgs args)
        {
            DispatchSettingUpdate(GraphElementID.GridLines);
        }

        private void HandleRegionLabelsToggle(object sender, EventArgs args)
        {
            DispatchSettingUpdate(GraphElementID.RegionLabels);
        }

        private void HandleConstantLinesToggle(object sender, EventArgs args)
        {
            DispatchSettingUpdate(GraphElementID.ConstantLines);
        }

        private void HandleAxisTrackersToggle(object sender, EventArgs args)
        {
            DispatchSettingUpdate(GraphElementID.AxisTrackers);
        }

        private void HandleMusicToggle(object sender, EventArgs args)
        {
            DispatchSettingUpdate(GraphElementID.Music);
        }

        private void HandleCreditsButtonPressed(object sender, EventArgs args)
        {
            EventMgr.Events.Dispatch(GameEvents.ClickDisplayCredits, false);

            SetMainPanelVisible(false);
            EventMgr.Events.Dispatch(GameEvents.TitleCreditsOpened);
            SetCreditsVisible(true);
        }

        private void HandleCloseCreditsButtonPressed(object sender, EventArgs args)
        {
            EventMgr.Events.Dispatch(GameEvents.ClickCloseCredits, false);

            m_credits.ManualReturnClick();
            SetMainPanelVisible(true);
            SetCreditsVisible(false);
        }

        private void HandleConfigureButtonPressed(object sender, EventArgs args)
        {
            SelectConfigureTab();
        }

        private void HandleControlsButtonPressed(object sender, EventArgs args)
        {
            SelectControlsTab();
        }

        #endregion // Handlers
    }

}