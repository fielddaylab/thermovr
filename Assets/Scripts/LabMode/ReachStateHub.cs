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
    public struct ReachStateDefinition
    {
        public string QuestionText;
        public List<SimStateTarget> Targets;

        public ReachStateDefinition(string qText, List<SimStateTarget> targets) {
            QuestionText = qText;
            Targets = targets;
        }
    }


    public struct SimStateTarget {
        public VarID TargetID;
        public float TargetVal;
        public float TargetRange;

        public SimStateTarget(VarID id, float val, float range) {
            TargetID = id;
            TargetVal = val;
            TargetRange = range;
        }
    }

    public class ReachStateHub : Evaluable
    {
        [SerializeField] private TMP_Text m_questionText;
        [SerializeField] private Image m_completionState;

        private ReachStateDefinition m_definition;

        public void SetDefinition(ReachStateDefinition def) {
            m_definition = def;

            ResetState();
        }

        #region IEvaluable

        public override void ResetState() {
            base.ResetState();

            m_questionText.SetText(m_definition.QuestionText);
            m_completionState.sprite = GameDB.Instance.SocketEmpty;
        }

        public override bool AnswerSelected() {
            return false;
        }

        public override void HandleEvaluation(bool correct) {
            if (m_evaluated) {
                // no need for duplicate evaluations
                return;
            }

            if (correct) {
                m_completionState.sprite = GameDB.Instance.Correct;

                m_evaluated = true;
            }
        }

        public override bool IsCorrect() {
            for (int i = 0; i < m_definition.Targets.Count; i++) {
                SimStateTarget currTarget = m_definition.Targets[i];

                double varVal = ThermoPresent.Instance.get_state_var(currTarget.TargetID);

                if (varVal < currTarget.TargetVal - currTarget.TargetRange || varVal > currTarget.TargetVal + currTarget.TargetRange) {
                    // is outside of limits
                    return false;
                }
            }

            // return true if no reqs were outside of limits
            return true;
        }

        #endregion // IEvaluable
    }


}