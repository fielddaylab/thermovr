using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR.Lab;
using ThermoVR.UI;
using TMPro;
using UnityEngine;
using UnityEngine.UI;

namespace ThermoVR
{
    public class MultipleChoiceOption : MonoBehaviour
    {
        public event EventHandler<IDEventArgs> OnChoiceSelected; // used by MultipleChoiceHub

        [SerializeField] private ThermoButton m_button;
        [SerializeField] private TMP_Text m_optionText;
        [SerializeField] private Image m_fill; // fill in bubble

        private uint m_choiceID; // unique identifier used for determining correct answer

        #region Unity Callbacks

        private void OnEnable() {
            m_button.OnButtonPressed += HandleButtonPressed;
        }

        private void OnDisable() {
            m_button.OnButtonPressed -= HandleButtonPressed;
        }

        #endregion // Unity Callbacks

        public void SetChoiceID(uint id) {
            m_choiceID = id;
        }

        public uint GetChoiceID() {
            return m_choiceID;
        }

        public void SetOptionText(string text) {
            m_optionText.text = text;
        }

        public void SetSelected(bool selected) {
            m_fill.enabled = selected;
        }

        public void SetEvaluatedState(EvalState state) { // black, green, or red img depending on evaluation state
            switch (state) {
                case EvalState.Pending:
                    m_fill.sprite = GameDB.Instance?.MCFill;
                    break;
                case EvalState.Correct:
                    m_fill.sprite = GameDB.Instance?.Correct;
                    break;
                case EvalState.Incorrect:
                    m_fill.sprite = GameDB.Instance?.Incorrect;
                    break;
                case EvalState.Missed:
                    m_fill.sprite = GameDB.Instance?.Missed;
                    break;
                default:
                    break;
            }
        }

        #region Handlers

        private void HandleButtonPressed(object sender, EventArgs args) {
            OnChoiceSelected?.Invoke(this, new IDEventArgs(m_choiceID));
        }

        #endregion // Handlers
    }
}