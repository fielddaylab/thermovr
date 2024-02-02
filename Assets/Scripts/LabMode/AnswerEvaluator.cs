using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR.UI;
using UnityEngine;

namespace ThermoVR.Lab
{
    public class BoolEventArgs : EventArgs {
        public bool Value;

        public BoolEventArgs(bool value) {
            Value = value; 
        }
    }

    public class AnswerEvaluator : MonoBehaviour
    {
        [SerializeField] private Evaluable[] m_toEvaluate;
        [SerializeField] private bool m_constantCheck = false;
        [SerializeField] private ThermoButton m_submitButton;

        public EventHandler<BoolEventArgs> OnEvaluationUpdated;

        private void OnEnable() {
            if (!m_constantCheck) {
                m_submitButton.OnButtonPressed += HandleSubmitPressed;
                m_submitButton.SetInteractable(false);
            }
        }

        private void OnDisable() {
            if (!m_constantCheck) {
                m_submitButton.OnButtonPressed -= HandleSubmitPressed;
            }
        }

        private void Update() {
            if (m_constantCheck) {
                bool allCorrect = true;

                for (int i = 0; i < m_toEvaluate.Length; i++) {
                    if (m_toEvaluate[i].HasBeenEvaluated()) {
                        // no need for evaluation again until reset
                        continue;
                    }
                    if (!m_toEvaluate[i].IsCorrect()) {
                        allCorrect = false;
                    }

                    // Only update result when requirements have been met
                    // if (m_toEvaluate[i].IsCorrect()) {
                        m_toEvaluate[i].HandleEvaluation(m_toEvaluate[i].IsCorrect());
                        OnEvaluationUpdated?.Invoke(this, new BoolEventArgs(m_toEvaluate[i].IsCorrect()));
                    // }
                }
            }

            bool allSelected = true;
            for (int i = 0; i < m_toEvaluate.Length; i++)
            {
                if (!m_toEvaluate[i].AnswerSelected())
                {
                    allSelected = false;
                }
            }

            m_submitButton.SetInteractable(allSelected);
        }

        private void HandleSubmitPressed(object sender, EventArgs args) {
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
                OnEvaluationUpdated?.Invoke(this, new BoolEventArgs(m_toEvaluate[i].IsCorrect()));
            }

            if (allCorrect) {
                Debug.Log("[Q] Correct!");
            }
            else {
                Debug.Log("[Q] Not crrect!");
            }
        }

        public bool HasBeenEvaluated()
        {
            bool eval = true;

            for (int i = 0; i < m_toEvaluate.Length; i++)
            {
                if (!m_toEvaluate[i].HasBeenEvaluated())
                {
                    eval = false;
                }
            }

            return eval;
        }

        public void LoadCompleted(bool completed)
        {
            for (int i = 0; i < m_toEvaluate.Length; i++)
            {
                m_toEvaluate[i].LoadCompleted(completed);
            }
        }
    }
}