using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR.UI;
using TMPro;
using UnityEngine;
using UnityEngine.UI;

namespace ThermoVR.Lab
{
    public struct WordBankDefinition
    {
        public string InitialConditionText;
        public string[] QuestionTexts;
        public string[] OptionTexts;
        public uint CorrectID; // between 0 and number of options - 1

        public WordBankDefinition(string initText, string[] qTexts, string[] optionTexts, uint id) {
            InitialConditionText = initText;
            QuestionTexts = qTexts;
            OptionTexts = optionTexts;
            CorrectID = id;
        }
    }

    public class WordBankHub : Evaluable
    {
        [SerializeField] private TMP_Text m_initText;
        // [SerializeField] private TMP_Text m_questionText;
        [SerializeField] private TMP_Text m_answerText;
        [SerializeField] private Image m_answerBG;
        [SerializeField] private ThermoButton m_chooseButton;
        [SerializeField] private bool m_randomOrder = false;

        [SerializeField] private InstructionLineGenerator m_lineGenerator;

        [SerializeField] private GameObject m_choicePanel;
        [SerializeField] private ThermoButton m_choicePanelCloseButton;
        [SerializeField] private WordBankOption[] m_options; // option "slots"; not all questions will use all slots // TODO: make pools


        private uint m_selectedID;

        private uint[] m_order;

        private WordBankDefinition m_definition;

        private void OnEnable() {
            m_chooseButton.OnButtonPressed += HandleChoosePressed;
        }

        private void OnDisable() {
            m_chooseButton.OnButtonPressed -= HandleChoosePressed;
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

        public override bool AnswerSelected() {
            return m_selectedID != uint.MaxValue;
        }

        public override void HandleEvaluation(bool correct) {
            if (m_evaluated) {
                // no need for duplicate evaluations
                return;
            }

            if (correct) {
                m_answerBG.color = Color.green;
            }
            else {
                m_answerBG.color = Color.red;
            }

            m_evaluated = true;
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
            m_selectedID = uint.MaxValue;
            m_answerBG.color = Color.white;
            m_choicePanel.SetActive(false);
        }

        #endregion // IEvaluable

        public void SetDefinition(WordBankDefinition def) {
            m_definition = def;

            ResetState();
        }

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

        private void UpdateOptions(uint[] order) {
            for (int i = 0; i < m_definition.QuestionTexts.Length; i++)
            {
                var transform = m_lineGenerator.GenerateLine(m_definition.QuestionTexts[i]);
                transform.localPosition += new Vector3(0, -1.3f * i, 0);
            }

            // show options equal to number of options
            for (int i = 0; i < m_options.Length; i++) {
                if (i >= m_definition.OptionTexts.Length) {
                    // Hide option
                    m_options[i].gameObject.SetActive(false);
                    continue;
                }

                m_options[i].gameObject.SetActive(true);
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

        #region Handlers

        private void HandleChoiceSelected(object sender, IDEventArgs args) {
            if (m_evaluated) {
                // no more inputs after evaluation (until reset)
                return;
            }

            m_choicePanel.SetActive(false);

            m_selectedID = args.ID;
            m_answerText.SetText(m_definition.OptionTexts[m_selectedID]);

            // Set answer text
        }

        private void HandleChoosePressed(object sender, EventArgs args) {
            if (m_evaluated) {
                // no more choices until reset
                return;
            }

            // display available choices
            m_choicePanel.SetActive(true);
            // show option buttons
            // add HandleCHoiceSelected
        }

        private void HandleChoicePanelClosePressed(object sender, EventArgs args) {
            
        }

        #endregion // Handlers
    }
}