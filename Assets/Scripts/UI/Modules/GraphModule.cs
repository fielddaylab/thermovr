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

    [SerializeField] private PhysicalButton m_axisNumbersToggle;
    [SerializeField] private PhysicalButton m_gridLinesToggle;
    [SerializeField] private PhysicalButton m_regionLabelsToggle;
    [SerializeField] private PhysicalButton m_constantLinesToggle;

    #endregion // Inspector

    #region IUIModule

    public override void Open() {
        base.Open();
    }

    public override void Close() {
        base.Close();
    }

    #endregion // IUIModule

    #region Handlers

    // TODO: store this data in the toggle.
    // Also see if you can disconnect toggles from general touch checking; otherwise incorporate into touch pipeline

    private void HandleAxisNumbersToggle() {
        // TODO: dispatch relevant toggle value, not always 'true'
        GameMgr.Events.Dispatch(GameEvents.UpdateGraphSetting, new GraphSettingUpdate(GraphElementID.AxisNumbers, true));
    }

    private void HandleGridLinesToggle() {
        // TODO: dispatch relevant toggle value, not always 'true'
        GameMgr.Events.Dispatch(GameEvents.UpdateGraphSetting, new GraphSettingUpdate(GraphElementID.GridLines, true));
    }

    private void HandleRegionLabelsToggle() {
        // TODO: dispatch relevant toggle value, not always 'true'
        GameMgr.Events.Dispatch(GameEvents.UpdateGraphSetting, new GraphSettingUpdate(GraphElementID.RegionLabels, true));
    }

    private void HandleConstantLinesToggle() {
        // TODO: dispatch relevant toggle value, not always 'true'
        GameMgr.Events.Dispatch(GameEvents.UpdateGraphSetting, new GraphSettingUpdate(GraphElementID.ConstantLines, true));
    }

    #endregion // Handlers
}
