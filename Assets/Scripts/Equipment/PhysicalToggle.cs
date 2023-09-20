using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR;
using UnityEngine;

public class PhysicalToggle : MonoBehaviour
{
    [SerializeField] private Pressable m_button;
    [SerializeField] private Material m_defaultMat, m_activeMat;
    [SerializeField] private MeshRenderer m_renderer;

    private bool m_isOn;

    private void Awake() {
        m_isOn = false;

        m_button.OnPress += HandleTogglePressed;
    }

    public bool IsOn() {
        return m_isOn;
    }

    #region Handlers

    private void HandleTogglePressed(object sender, EventArgs args) {
        m_isOn = !m_isOn;
        if (m_isOn) {
            if (m_renderer != null && m_renderer.materials.Length > 1) {
                var mats = m_renderer.materials;
                mats[1] = m_activeMat;
                m_renderer.materials = mats;
            }
        }
        else {
            if (m_renderer != null && m_renderer.materials.Length > 1) {
                var mats = m_renderer.materials;
                mats[1] = m_defaultMat;
                m_renderer.materials = mats;
            }
        }
    }

    #endregion // Handlers
}
