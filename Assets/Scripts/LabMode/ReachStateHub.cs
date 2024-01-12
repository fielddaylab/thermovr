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
        public string InitialConditionText;
        public string[] QuestionTexts;
        public List<SimStateTarget> Targets;

        public ReachStateDefinition(string initText, string[] qTexts, List<SimStateTarget> targets) {
            InitialConditionText = initText;
            QuestionTexts = qTexts;
            Targets = targets;
        }
    }


    [Serializable]
    public struct SimStateTarget {
        public VarID TargetID;
        public float TargetVal;
        public float TargetRange;
        // TODO: targetComparison; ( less than, equal to, greater than, etc.)

        public SimStateTarget(VarID id, float val, float range) {
            TargetID = id;
            TargetVal = val;
            TargetRange = range;
        }
    }

    public class ReachStateHub : Evaluable
    {
        [SerializeField] private TMP_Text m_initText;
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

            m_initText.SetText(m_definition.InitialConditionText);

            string questionStr = "";
            for (int i = 0; i < m_definition.QuestionTexts.Length; i++)
            {
                questionStr += m_definition.QuestionTexts[i];
                questionStr += "\n";
            }
            m_questionText.SetText(questionStr);

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

                if (currTarget.TargetID == VarID.VolumeStop) {
                    bool eitherCorrect = false;
                    Tuple<double, double> stopVals = World.Instance.get_stop_vals();
                    bool stop1OutOfRange = (stopVals.Item1 < currTarget.TargetVal - currTarget.TargetRange || stopVals.Item1 > currTarget.TargetVal + currTarget.TargetRange);
                    bool stop2OutOfRange = (stopVals.Item2 < currTarget.TargetVal - currTarget.TargetRange || stopVals.Item2 > currTarget.TargetVal + currTarget.TargetRange);
                    if (stop1OutOfRange && stop2OutOfRange) {
                        // both stops are outside of limits
                        return false;
                    }
                }
                else {
                    double varVal = World.Instance.get_state_var(currTarget.TargetID);

                    if (varVal < currTarget.TargetVal - currTarget.TargetRange || varVal > currTarget.TargetVal + currTarget.TargetRange) {
                        // is outside of limits
                        return false;
                    }
                }
            }

            // return true if no reqs were outside of limits
            return true;
        }

        #endregion // IEvaluable
    }


}