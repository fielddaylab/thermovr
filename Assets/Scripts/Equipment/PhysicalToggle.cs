using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR;
using UnityEngine;

public class PhysicalToggle : MonoBehaviour
{
    [SerializeField] private Pressable m_button;

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
    }

    #endregion // Handlers
}
