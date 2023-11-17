using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR;
using UnityEngine;

public class PhysicalToggle : MonoBehaviour
{
    [SerializeField] private Pressable m_button;
    [SerializeField] private MeshRenderer m_renderer;

    private bool m_isOn;

    private void Awake() {
        m_isOn = false;

        m_button.OnPress += HandleTogglePressed;
    }

    public bool IsOn() {
        return m_isOn;
    }

    public void ResetToggle() {
        m_isOn = false;
        UpdateActiveMaterial();
    }

    private void UpdateActiveMaterial() {
        if (m_isOn) {
            if (m_renderer != null && m_renderer.materials.Length > 1) {
                var mats = m_renderer.materials;
                mats[1] = GameDB.Instance.ActiveButtonMaterial;
                m_renderer.materials = mats;
            }
        }
        else {
            if (m_renderer != null && m_renderer.materials.Length > 1) {
                var mats = m_renderer.materials;
                mats[1] = GameDB.Instance.InactiveButtonMaterial;
                m_renderer.materials = mats;
            }
        }
    }

    #region Handlers

    private void HandleTogglePressed(object sender, EventArgs args) {
        m_isOn = !m_isOn;
        UpdateActiveMaterial();
        m_button.ClickAudio.Play();
    }

    #endregion // Handlers
}
