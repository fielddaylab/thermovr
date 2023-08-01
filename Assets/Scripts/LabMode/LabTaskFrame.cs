using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR.UI;
using UnityEngine;

namespace ThermoVR.Lab
{
    public class LabTaskFrame : MonoBehaviour
    {
        public AnswerEvaluator AnswerEvaluator;
        public ThermoButton TaskResetButton;

        [SerializeField] private Evaluable[] m_evaluables;

        private void OnEnable() {
            TaskResetButton.OnButtonPressed += HandleResetPressed;
        }

        private void OnDisable() {
            TaskResetButton.OnButtonPressed -= HandleResetPressed;
        }

        public Evaluable[] GetEvaluables() {
            return m_evaluables;
        }

        private void HandleResetPressed(object sender, EventArgs args) {
            for (int i = 0; i < m_evaluables.Length; i++) {
                m_evaluables[i].ResetState();
            }
        }
    }
}

