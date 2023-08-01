using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR.UI;
using UnityEngine;

namespace ThermoVR.Lab
{
    public class AnswerEvaluator : MonoBehaviour
    {
        [SerializeField] private Evaluable[] m_toEvaluate;
        [SerializeField] private ThermoButton m_submitButton;

        private void OnEnable() {
            m_submitButton.OnButtonPressed += HandleSubmitPressed;
        }

        private void OnDisable() {
            m_submitButton.OnButtonPressed -= HandleSubmitPressed;
        }

        private void HandleSubmitPressed(object sender, EventArgs args) {
            Debug.Log("[Q] submit pressed");
            // Only submit when answers have been selected
            for (int i = 0; i < m_toEvaluate.Length; i++) {
                if (!m_toEvaluate[i].AnswerSelected()) {
                    return;
                }
            }

            bool allCorrect = true;

            for (int i = 0; i < m_toEvaluate.Length; i++) {
                if (!m_toEvaluate[i].IsCorrect()) {
                    allCorrect = false;
                }

                m_toEvaluate[i].HandleEvaluation(m_toEvaluate[i].IsCorrect());
            }

            if (allCorrect) {
                Debug.Log("[Q] Correct!");
            }
            else {
                Debug.Log("[Q] Not crrect!");
            }
        }
    }
}