using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR.Lab;
using ThermoVR.UI;
using TMPro;
using UnityEngine;

namespace ThermoVR.Lab
{
    public class WordBankOption : MonoBehaviour
    {
        public event EventHandler<IDEventArgs> OnChoiceSelected; // used by WordBankHub

        [SerializeField] private ThermoButton m_button;


        private uint m_choiceID; // unique identifier used for determining correct answer

        public void SetChoiceID(uint id) {
            m_choiceID = id;
        }

        public uint GetChoiceID() {
            return m_choiceID;
        }

        public void SetOptionText(string text) {
            m_button.SetText(text);
        }

        public string GetOptionText()
        {
            return m_button.GetText();
        }

        #region Unity Callbacks

        private void OnEnable() {
            m_button.OnButtonPressed += HandleButtonPressed;
        }

        private void OnDisable() {
            m_button.OnButtonPressed -= HandleButtonPressed;
        }

        #endregion // Unity Callbacks

        #region Handlers

        private void HandleButtonPressed(object sender, EventArgs args) {
            OnChoiceSelected?.Invoke(this, new IDEventArgs(m_choiceID));
        }

        #endregion // Handlers
    }
}
