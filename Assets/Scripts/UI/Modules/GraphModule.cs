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

    private ThermoToggle[] m_toggles;

    public override void Init() {
        base.Init();

        m_toggles = new ThermoToggle[4] {
            m_axisNumbersToggle,
            m_gridLinesToggle,
            m_regionLabelsToggle,
            m_constantLinesToggle
        };

        for (int i = 0; i < m_toggles.Length; i++) {
            m_toggles[i].Init();
        }
    }
    #endregion // Inspector

    #region IUIModule

    public override void Open() {
        base.Open();

        AddListeners();
    }

    public override void Close() {
        base.Close();

        RemoveListeners();
    }

    #endregion // IUIModule

    #region Helpers

    private void AddListeners() {
        m_axisNumbersToggle.Pressable.PressCompleted += HandleAxisNumbersToggle;
        m_gridLinesToggle.Pressable.PressCompleted += HandleGridLinesToggle;
        m_regionLabelsToggle.Pressable.PressCompleted += HandleRegionLabelsToggle;
        m_constantLinesToggle.Pressable.PressCompleted += HandleConstantLinesToggle;
    }

    private void RemoveListeners() {
        m_axisNumbersToggle.Pressable.PressCompleted -= HandleAxisNumbersToggle;
        m_gridLinesToggle.Pressable.PressCompleted -= HandleGridLinesToggle;
        m_regionLabelsToggle.Pressable.PressCompleted -= HandleRegionLabelsToggle;
        m_constantLinesToggle.Pressable.PressCompleted -= HandleConstantLinesToggle;
    }

    #endregion // Helpers

    #region Handlers

    // Also see if you can disconnect toggles from general touch checking; otherwise incorporate into touch pipeline

    private void HandleAxisNumbersToggle(object sender, EventArgs args) {
        bool val = m_axisNumbersToggle.IsOn();
        GameMgr.Events.Dispatch(GameEvents.UpdateGraphSetting, new GraphSettingUpdate(GraphElementID.AxisNumbers, val));
    }

    private void HandleGridLinesToggle(object sender, EventArgs args) {
        bool val = m_gridLinesToggle.IsOn();
        GameMgr.Events.Dispatch(GameEvents.UpdateGraphSetting, new GraphSettingUpdate(GraphElementID.GridLines, val));
    }

    private void HandleRegionLabelsToggle(object sender, EventArgs args) {
        bool val = m_regionLabelsToggle.IsOn();
        GameMgr.Events.Dispatch(GameEvents.UpdateGraphSetting, new GraphSettingUpdate(GraphElementID.RegionLabels, val));
    }

    private void HandleConstantLinesToggle(object sender, EventArgs args) {
        bool val = m_constantLinesToggle.IsOn();
        GameMgr.Events.Dispatch(GameEvents.UpdateGraphSetting, new GraphSettingUpdate(GraphElementID.ConstantLines, val));
    }

    #endregion // Handlers
}
