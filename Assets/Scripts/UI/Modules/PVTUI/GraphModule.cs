using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR;
using ThermoVR.Controls;
using ThermoVR.UI;
using ThermoVR.UI.GraphElements;
using ThermoVR.UI.Interfaces;
using UnityEngine;

public class GraphModule : UIModule
{
    #region Inspector

    // TODO: make these toggles and assign in the editor

    [SerializeField] private ThermoToggle m_axisNumbersToggle;
    [SerializeField] private ThermoToggle m_gridLinesToggle;
    [SerializeField] private ThermoToggle m_regionLabelsToggle;
    [SerializeField] private ThermoToggle m_constantLinesToggle;
    [SerializeField] private ThermoToggle m_axisTrackersToggle;

    private ThermoToggle[] m_toggles;

    public override void Init() {
        base.Init();

        m_toggles = new ThermoToggle[5] {
            m_axisNumbersToggle,
            m_gridLinesToggle,
            m_regionLabelsToggle,
            m_constantLinesToggle,
            m_axisTrackersToggle
        };

        for (int i = 0; i < m_toggles.Length; i++) {
            m_toggles[i].Init();
        }

        Array settings = Enum.GetValues(typeof(GraphElementID));
        for (int i = 0; i < settings.Length; i++) {
            DispatchSettingUpdate((GraphElementID)settings.GetValue(i));
        }
    }
    #endregion // Inspector

    #region IUIModule

    public override void Open() {
        this.gameObject.SetActive(true);

        AddListeners();
    }

    public override void Close() {
        this.gameObject.SetActive(false);

        RemoveListeners();
    }

    #endregion // IUIModule

    #region Helpers

    private void AddListeners() {
        m_axisNumbersToggle.Pressable.PressCompleted += HandleAxisNumbersToggle;
        m_gridLinesToggle.Pressable.PressCompleted += HandleGridLinesToggle;
        m_regionLabelsToggle.Pressable.PressCompleted += HandleRegionLabelsToggle;
        m_constantLinesToggle.Pressable.PressCompleted += HandleConstantLinesToggle;
        m_axisTrackersToggle.Pressable.PressCompleted += HandleAxisTrackersToggle;
    }

    private void RemoveListeners() {
        m_axisNumbersToggle.Pressable.PressCompleted -= HandleAxisNumbersToggle;
        m_gridLinesToggle.Pressable.PressCompleted -= HandleGridLinesToggle;
        m_regionLabelsToggle.Pressable.PressCompleted -= HandleRegionLabelsToggle;
        m_constantLinesToggle.Pressable.PressCompleted -= HandleConstantLinesToggle;
        m_axisTrackersToggle.Pressable.PressCompleted -= HandleAxisTrackersToggle;
    }

    private void DispatchSettingUpdate(GraphElementID id) {
        ThermoToggle toggle;

        switch (id) {
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
            default:
                return;
        }

        bool val = toggle.IsOn();
        GameMgr.Events.Dispatch(GameEvents.UpdateGraphSetting, new GraphSettingUpdate(id, val));
    }

    #endregion // Helpers

    #region Handlers

    // Also see if you can disconnect toggles from general touch checking; otherwise incorporate into touch pipeline

    private void HandleAxisNumbersToggle(object sender, EventArgs args) {
        DispatchSettingUpdate(GraphElementID.AxisNumbers);
    }

    private void HandleGridLinesToggle(object sender, EventArgs args) {
        DispatchSettingUpdate(GraphElementID.GridLines);
    }

    private void HandleRegionLabelsToggle(object sender, EventArgs args) {
        DispatchSettingUpdate(GraphElementID.RegionLabels);
    }

    private void HandleConstantLinesToggle(object sender, EventArgs args) {
        DispatchSettingUpdate(GraphElementID.ConstantLines);
    }

    private void HandleAxisTrackersToggle(object sender, EventArgs args) {
        DispatchSettingUpdate(GraphElementID.AxisTrackers);
    }

    #endregion // Handlers
}
