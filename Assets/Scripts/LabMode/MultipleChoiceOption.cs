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

        [SerializeField] private BoxCollider m_collider;
        [SerializeField] private Vector3 m_ExemplaryMCOptionCenter;
        [SerializeField] private Vector3 m_ExemplaryMCOptionSize;

        [SerializeField] private ThermoButton m_button;
        [SerializeField] private TMP_Text m_optionText;
        [SerializeField] private Image m_fill; // fill in bubble
        [SerializeField] private Graphic m_bg; // background rect
        [SerializeField] private Graphic m_thermoBtnImg;

        private uint m_choiceID; // unique identifier used for determining correct answer

        #region Unity Callbacks

        private void OnEnable() {
            m_button.OnButtonPressed += HandleButtonPressed;

            FixCollider();
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
            m_bg.color = selected ? GameDB.Instance.MCSelectedBG : GameDB.Instance.MCUnselectedBG;
        }

        public void SetEvaluatedState(EvalState state) { // black, green, or red img depending on evaluation state
            switch (state) {
                case EvalState.Pending:
                    m_fill.sprite = GameDB.Instance?.MCFill;
                    m_bg.color = GameDB.Instance.MCSelectedBG;
                    m_thermoBtnImg.enabled = true;
                    break;
                case EvalState.Correct:
                    m_fill.sprite = GameDB.Instance?.Correct;
                    m_bg.color = GameDB.Instance.MCSelectedBG;
                    m_thermoBtnImg.enabled = false;
                    break;
                case EvalState.Incorrect:
                    m_fill.sprite = GameDB.Instance?.Incorrect;
                    m_bg.color = GameDB.Instance.MCIncorrectBG;
                    m_thermoBtnImg.enabled = false;
                    break;
                case EvalState.Missed:
                    m_fill.sprite = GameDB.Instance?.Missed;
                    m_bg.color = GameDB.Instance.MCIncorrectBG;
                    m_thermoBtnImg.enabled = false;
                    break;
                default:
                    break;
            }
        }

        public void FixCollider()
        {
            m_collider.center = m_ExemplaryMCOptionCenter;
            m_collider.size = m_ExemplaryMCOptionSize;
        }

        #region Handlers

        private void HandleButtonPressed(object sender, EventArgs args) {
            OnChoiceSelected?.Invoke(this, new IDEventArgs(m_choiceID));
        }

        #endregion // Handlers
    }
}