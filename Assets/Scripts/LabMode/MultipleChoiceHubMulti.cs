using BeauUtil;
using JetBrains.Annotations;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using ThermoVR.UI;
using TMPro;
using UnityEngine;
using UnityEngine.InputSystem.LowLevel;

namespace ThermoVR.Lab
{
    public class MultipleChoiceHubMulti : Evaluable
    {

        [SerializeField] private TMP_Text m_initText;
        // [SerializeField] private TMP_Text m_questionText;
        [SerializeField] private MultipleChoiceOption[] m_options; // option "slots"; not all questions will use all slots // TODO: make pools
        [SerializeField] private bool m_randomOrder = false;

        [SerializeField] private InstructionLineGenerator m_lineGenerator;

        private MultipleChoiceDefinition m_definition;

        private List<uint> m_selectedIDs;

        private uint[] m_order;

        #region Unity Callbacks

        private void OnEnable() {
            if (m_evaluated) {
                // no need to re-order or re-load
            }
            else {
                // ResetState();
            }

            if (m_selectedIDs == null) {
                m_selectedIDs = new List<uint>();
            }

            for (int i = 0; i < m_options.Length; i++)
            {
                m_options[i].FixCollider();
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
                m_options[i].OnChoiceSelected -= HandleChoiceSelected;
                m_options[i].OnChoiceSelected += HandleChoiceSelected;

                // set option text and id according to order
                uint choiceOptionIndex = order[i];
                m_options[i].SetChoiceID(choiceOptionIndex);
                m_options[i].SetOptionText(m_definition.OptionTexts[choiceOptionIndex]);

                m_options[i].FixCollider();
            }
        }

        public void SetDefinition(MultipleChoiceDefinition def) {
            m_definition = def;

            ResetState();
        }

        #region IEvaluable

        public override bool IsCorrect() {
            bool mismatchFound = false;

            // for each response id
            for (int i = 0; i < m_definition.OptionTexts.Length; i++) {
                uint optionID = m_options[i].GetChoiceID();

                // see if it was selected; if so, see if it is correct; if not, see if not correct
                if (m_definition.CorrectIDs.Contains(optionID) != m_selectedIDs.Contains(optionID)) {
                    mismatchFound = true;
                    break;
                }
            }

            return !mismatchFound;
        }

        public override bool AnswerSelected() {
            return m_selectedIDs.Count > 0;
        }

        public override void ResetState() {
            base.ResetState();

            SetOrder(m_randomOrder);

            m_initText.SetText(m_definition.InitialConditionText);

            for (int i = 0; i < m_definition.QuestionTexts.Length; i++)
            {
                var transform = m_lineGenerator.GenerateLine(m_definition.QuestionTexts[i]);
                transform.localPosition += new Vector3(0, -1.3f * i, 0);
            }

            UpdateOptions(m_order);
            m_selectedIDs.Clear();
        }

        public override void HandleEvaluation(bool correct) {
            if (m_evaluated) {
                // no need for duplicate evaluations
                return;
            }

            // for each response id
            for (int i = 0; i < m_order.Length; i++) {
                uint optionID = m_options[m_order[i]].GetChoiceID();

                // see if it was selected; if so, see if it is correct; if not, see if not correct
                if (m_definition.CorrectIDs.Contains(optionID) == m_selectedIDs.Contains(optionID)) {
                    m_options[m_order[i]].SetEvaluatedState(EvalState.Correct);
                }
                else if (m_definition.CorrectIDs.Contains(optionID) && !m_selectedIDs.Contains(optionID)) {
                    m_options[m_order[i]].SetEvaluatedState(EvalState.Missed);
                    m_options[i].SetSelected(true);

                }
                else if (m_selectedIDs.Contains(optionID)) {
                    m_options[m_order[i]].SetEvaluatedState(EvalState.Incorrect);
                }
                else {
                    // not selected, and didn't need to be selected. No change.

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

            bool selectedState = false;
            if (m_selectedIDs.Contains(args.ID)) {
                m_selectedIDs.Remove(args.ID);
                selectedState = false;
            }
            else {
                m_selectedIDs.Add(args.ID);
                selectedState = true;
            }

            for (int i = 0; i < m_order.Length; i++) {
                if (m_options[i].GetChoiceID() == args.ID) {
                    m_options[i].SetSelected(selectedState);
                }
            }
        }

        #endregion // Handlers
    }
}