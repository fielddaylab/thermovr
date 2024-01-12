using BeauUtil;
using JetBrains.Annotations;
using System;
using System.Collections;
using System.Collections.Generic;
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

    public enum EvalState
    {
        Pending,
        Correct,
        Incorrect,
        Missed
    }

    public struct MultipleChoiceDefinition
    {
        public string InitialConditionText;
        public string[] QuestionTexts;
        public string[] OptionTexts;
        public List<uint> CorrectIDs; // between 0 and number of options - 1

        public MultipleChoiceDefinition(string initText, string[] qTexts, string[] optionTexts, List<uint> ids) {
            InitialConditionText = initText;
            QuestionTexts = qTexts;
            OptionTexts = optionTexts;
            CorrectIDs = ids;
        }
    }

    public class MultipleChoiceHub : Evaluable
    {
        [SerializeField] private TMP_Text m_initText;
        [SerializeField] private TMP_Text m_questionText;
        [SerializeField] private MultipleChoiceOption[] m_options; // option "slots"; not all questions will use all slots // TODO: make pools
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
                if (m_definition.CorrectIDs[0] == m_selectedID) {
                    return true;
                }
            }

            return false;
        }

        public override bool AnswerSelected() {
            return m_selectedID != uint.MaxValue;
        }

        public override void ResetState() {
            base.ResetState();

            SetOrder(m_randomOrder);

            m_initText.SetText(m_definition.InitialConditionText);

            string questionStr = "";
            for (int i = 0; i < m_definition.QuestionTexts.Length; i++)
            {
                questionStr += m_definition.QuestionTexts[i];
                questionStr += "\n";
            }
            m_questionText.SetText(questionStr);

            UpdateOptions(m_order);
            m_selectedID = uint.MaxValue;
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