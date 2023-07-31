using BeauUtil;
using JetBrains.Annotations;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using ThermoVR.UI;
using TMPro;
using UnityEngine;

namespace ThermoVR.Lab
{
    public class IDEventArgs : EventArgs
    {
        public uint ID;

        public IDEventArgs(uint id) {
            ID = id;
        }
    }

    public enum EvalState {
        Pending,
        Correct,
        Incorrect
    }

    public struct MultipleChoiceDefinition
    {
        public string QuestionText;
        public string[] OptionTexts;
        public uint CorrectID; // between 0 and number of options - 1

        public MultipleChoiceDefinition(string qText, string[] optionTexts, uint id) {
            QuestionText = qText;
            OptionTexts = optionTexts;
            CorrectID = id;
        }
    }

    public class MultipleChoiceHub : Evaluable
    {

        [SerializeField] private TMP_Text m_questionText;
        [SerializeField] private MultipleChoiceOption[] m_options; // option "slots"; not all questions will use all slots
        [SerializeField] private bool m_randomOrder = false;

        private MultipleChoiceDefinition m_definition;

        private uint m_selectedID;

        private uint[] m_order;

        #region Unity Callbacks

        private void OnEnable() {
            if (m_evaluated) {
                // no need to re-order or re-load
            }
            else {
                // ResetState();
            }
        }

        #endregion // Unity Callbacks

        private void UpdateOptions(uint[] order) {
            m_questionText.SetText(m_definition.QuestionText);

            // show options equal to number of options
            for (int i = 0; i < m_options.Length; i++) {
                if (i >= m_definition.OptionTexts.Length) {
                    // Hide option
                    m_options[i].gameObject.SetActive(false);
                    continue;
                }

                m_options[i].gameObject.SetActive(true);
                m_options[i].SetEvaluatedState(EvalState.Pending);
                m_options[i].SetSelected(false);
            }

            // Load options in the new order
            for (int i = 0; i < order.Length; i++) {
                m_options[i].OnChoiceSelected += HandleChoiceSelected;

                // set option text and id according to order
                uint choiceOptionIndex = order[i];
                m_options[i].SetChoiceID(choiceOptionIndex);
                m_options[i].SetOptionText(m_definition.OptionTexts[choiceOptionIndex]);
            }
        }

        public void SetDefinition(MultipleChoiceDefinition def) {
            m_definition = def;

            ResetState();
        }

        #region IEvaluable

        public override bool IsCorrect() {
            for (int i = 0; i < m_definition.OptionTexts.Length; i++) {
                if (m_definition.CorrectID == m_selectedID) {
                    return true;
                }
            }

            return false;
        }

        public override void ResetState() {
            base.ResetState();

            SetOrder(m_randomOrder);
            UpdateOptions(m_order);
        }

        public override void HandleEvaluation(bool correct) {
            if (m_evaluated) {
                // no need for duplicate evaluations
                return;
            }

            for (int i = 0; i < m_order.Length; i++) {
                if (m_options[i].GetChoiceID() == m_selectedID) {
                    EvalState newState = correct ? EvalState.Correct : EvalState.Incorrect;
                    m_options[i].SetEvaluatedState(newState);
                }
            }

            m_evaluated = true;
        }

        #endregion // IEvaluable

        private void SetOrder(bool random) {
            m_order = new uint[m_definition.OptionTexts.Length];

            if (random) {
                // TODO: Shuffle
            }
            else {
                // default order
                for (uint i = 0; i < m_order.Length; i++) {
                    m_order[i] = i;
                }
            }
        }

        #region Handlers

        private void HandleChoiceSelected(object sender, IDEventArgs args) {
            if (m_evaluated) {
                // no more inputs after evaluation (until reset)
                return;
            }

            m_selectedID = args.ID;

            // Reset selected state
            for (int i = 0; i < m_order.Length; i++) {
                m_options[i].SetSelected(m_options[i].GetChoiceID() == m_selectedID);
            }
        }

        #endregion // Handlers
    }
}