using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR;
using ThermoVR.Tools;
using UnityEngine;

public class ToolActivator : MonoBehaviour
{
    private const int BASE_MAT_INDEX = 1;

    [SerializeField] private Pressable m_button;
    [SerializeField] private MeshRenderer m_mesh;

    private List<Tool> m_tools;

    private void Start() {
        m_button.OnPress += HandleActivatorPressed;
    }

    public void SetTools(List<Tool> tools) {
        m_tools = tools;
    }

    public void UpdateActiveMaterial() {
        Material toSet = GameDB.Instance.InactiveButtonMaterial; // default

        foreach (var tool in m_tools) {
            if (tool.engaged) {
                toSet = GameDB.Instance.ActiveButtonMaterial; // default
                break;
            }
        }

        Material[] materials = m_mesh.materials;
        materials[BASE_MAT_INDEX] = toSet;
        m_mesh.materials = materials;
    }

    #region Handlers

    private void HandleActivatorPressed(object sender, EventArgs args) {
        if (m_tools == null || m_tools.Count == 0) {
            return;
        }

        foreach(var tool in m_tools) {
            UpdateActiveMaterial();

            GameMgr.Events.Dispatch(GameEvents.PressedToolToggle, tool);
        }
    }

    #endregion // Handlers
}
