using BeauUtil;
using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR.Lab;
using ThermoVR.State;
using TMPro;
using UnityEngine;
using UnityEngine.UI;

namespace ThermoVR.Lab
{
    public struct TextTaskDefinition
    {
        public string MainText;

        public TextTaskDefinition(string mainText)
        {
            MainText = mainText;
        }
    }

    public class TextHub : Evaluable
    {
        [SerializeField] private TMP_Text m_mainText;

        private TextTaskDefinition m_definition;

        public void SetDefinition(TextTaskDefinition def)
        {
            m_definition = def;

            m_mainText.SetText(def.MainText);
            
            ResetState();
        }

        #region IEvaluable

        public override void ResetState()
        {
            base.ResetState();

            m_mainText.SetText(m_definition.MainText);
        }

        public override bool AnswerSelected()
        {
            return true;
        }

        public override void HandleEvaluation(bool correct)
        {
            if (m_evaluated)
            {
                // no need for duplicate evaluations
                return;
            }

            m_evaluated = correct;
        }

        public override bool IsCorrect()
        {
            return true;
        }

        #endregion // IEvaluable
    }


}