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
    [SerializeField] private BoxCollider m_collider;
    [SerializeField] private ButtonPressMovement m_buttonMovement;

    private List<Tool> m_tools;

    private void Start() {
        m_button.OnPress += HandleActivatorPressed;
    }

    public void SetTools(List<Tool> tools) {
        m_tools = tools;
    }

    public void UpdateActiveMaterial() {
        bool active = false;
        Material toSet = GameDB.Instance.InactiveButtonMaterial; // default

        foreach (var tool in m_tools) {
            if (tool.engaged) {
                active = true;
                toSet = GameDB.Instance.ActiveButtonMaterial; // default
                break;
            }
        }

        Material[] materials = m_mesh.materials;
        materials[BASE_MAT_INDEX] = toSet;
        m_mesh.materials = materials;

        // update button indentation
        if (m_buttonMovement)
        {
            if (active)
            {
                m_buttonMovement.Indent();
            }
            else
            {
                m_buttonMovement.ResetPosition();
            }
        }
    }

    #region Handlers

    private void HandleActivatorPressed(object sender, EventArgs args) {
        if (m_tools == null || m_tools.Count == 0) {
            return;
        }

        foreach(var tool in m_tools) {
            if (!tool.allowed) {
                continue;
            }

            UpdateActiveMaterial();

            if (GameMgr.I.AudioEnabled) { m_button.ClickAudio.Play(); }
            
            EventMgr.Events.Dispatch(GameEvents.PressedToolToggle, tool);
        }
    }

    #endregion // Handlers

#if UNITY_EDITOR
    [ContextMenu("Apply Desktop Collider Size")]
    private void ApplyDesktopColliderSize()
    {
        m_collider.center = Vector3.zero;
        m_collider.size = new Vector3(0.05139474f, 0.016f, 0.07354519f);
    }
#endif
}
