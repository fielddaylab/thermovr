using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR;
using ThermoVR.Tools;
using UnityEngine;

public class ToolActivator : MonoBehaviour
{
    [SerializeField] private Pressable m_button;
    [SerializeField] private Tool m_tool;

    private void Start() {
        m_button.OnPress += HandleActivatorPressed;
    }

    #region Handlers

    private void HandleActivatorPressed(object sender, EventArgs args) {
        GameMgr.Events.Dispatch(GameEvents.PressedToolToggle, m_tool);
    }

    #endregion // Handlers
}
