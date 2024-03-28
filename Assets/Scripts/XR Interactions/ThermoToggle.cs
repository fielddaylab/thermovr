using BeauUtil;
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
        [Required] [SerializeField] private ClickToggle m_toggle;

        public void Init() {
            Pressable = GetComponent<Pressable>();

            Pressable.OnPress += HandlePress;
            m_toggle.PointerClicked.AddListener(HandlePointerClick);
        }

        public bool IsOn() {
            return m_toggle.isOn;
        }

        #region Handlers

        private void HandlePress(object sender, EventArgs args) {
            m_toggle.isOn = !m_toggle.isOn;
            m_toggle.onValueChanged?.Invoke(m_toggle.isOn);
        }

        private void HandlePointerClick()
        {
            // TODO: handle pointer clicks
            Pressable.Press(false, true, false);
        }

        #endregion // Handlers
    }
}