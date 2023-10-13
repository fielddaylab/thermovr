using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

namespace ThermoVR.UI
{
    [RequireComponent(typeof(Pressable))]
    public class ThermoToggle : MonoBehaviour
    {
        [HideInInspector] public Pressable Pressable;
        [SerializeField] private Toggle m_toggle;

        public void Init() {
            Pressable = GetComponent<Pressable>();

            Pressable.OnPress += HandlePress;
        }

        public bool IsOn() {
            return m_toggle.isOn;
        }

        #region Handlers

        private void HandlePress(object sender, EventArgs args) {
            m_toggle.isOn = !m_toggle.isOn;
            m_toggle.onValueChanged?.Invoke(m_toggle.isOn);
        }

        #endregion // Handlers
    }
}