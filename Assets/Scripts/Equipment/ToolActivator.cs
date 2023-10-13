using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR;
using ThermoVR.Tools;
using UnityEngine;

public class ToolActivator : MonoBehaviour
{
    [SerializeField] private Pressable m_button;
    private List<Tool> m_tools;

    private void Start() {
        m_button.OnPress += HandleActivatorPressed;
    }

    public void SetTools(List<Tool> tools) {
        m_tools = tools;
    }

    #region Handlers

    private void HandleActivatorPressed(object sender, EventArgs args) {
        if (m_tools == null || m_tools.Count == 0) {
            return;
        }
        foreach(var tool in m_tools) {
            GameMgr.Events.Dispatch(GameEvents.PressedToolToggle, tool);
        }
    }

    #endregion // Handlers
}
