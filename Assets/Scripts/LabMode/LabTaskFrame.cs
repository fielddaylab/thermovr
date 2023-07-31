using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR.UI;
using UnityEngine;

namespace ThermoVR.Lab
{
    public class LabTaskFrame : MonoBehaviour
    {
        [SerializeField] private ThermoButton m_taskReset;
        [SerializeField] private Evaluable[] m_evaluables;

        private void OnEnable() {
            m_taskReset.OnButtonPressed += HandleResetPressed;
        }

        private void OnDisable() {
            m_taskReset.OnButtonPressed -= HandleResetPressed;
        }

        private void HandleResetPressed(object sender, EventArgs args) {
            for (int i = 0; i < m_evaluables.Length; i++) {
                m_evaluables[i].ResetState();
            }
        }
    }
}

