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

    public enum ReachStateState
    {
        Incomplete,
        Countdown,
        Complete
    }

    public class ReachStateHub : Evaluable
    {
        [SerializeField] private TMP_Text m_initText;
        // [SerializeField] private TMP_Text m_questionText;
        [SerializeField] private Image m_completionStateImg;

        [SerializeField] private InstructionLineGenerator m_lineGenerator;

        private ReachStateDefinition m_definition;

        private ReachStateState m_completionState;

        private float m_completionTime;
        private float m_completionTimer;

        public void SetDefinition(ReachStateDefinition def) {
            m_definition = def;

            m_completionTime = 4f;

            ResetState();
        }

        private void Update()
        {
            if (IsWithinRange())
            {
                if (m_completionState == ReachStateState.Incomplete)
                {
                    // start timer
                    m_completionTimer = m_completionTime;
                    m_completionStateImg.fillAmount = 0;

                    m_completionState = ReachStateState.Countdown;
                }
                else if (m_completionState == ReachStateState.Countdown)
                {
                    // continue timer
                    m_completionTimer -= Time.deltaTime;
                    m_completionStateImg.fillAmount = 1 - (m_completionTimer / m_completionTime);

                    if (m_completionTimer <= 0)
                    {
                        m_completionState = ReachStateState.Complete;
                        m_completionStateImg.fillAmount = 1;
                    }
                }
                else
                {
                    // complete timer
                    m_completionStateImg.sprite = GameDB.Instance.ReachStateComplete;
                    m_completionStateImg.fillAmount = 1;

                    GameMgr.Events.Dispatch(GameEvents.TargetStateReached);
                    m_completionState = ReachStateState.Complete;
                }
            }
            else
            {
                m_completionStateImg.sprite = GameDB.Instance.ReachStateIncomplete;
                m_completionStateImg.fillAmount = 1;

                if (m_completionState != ReachStateState.Incomplete)
                {
                    GameMgr.Events.Dispatch(GameEvents.TargetStateLost, GetDiscrepancies());
                    m_completionState = ReachStateState.Incomplete;
                }
            }
        }

        #region IEvaluable

        public override void ResetState() {
            base.ResetState();

            m_initText.SetText(m_definition.InitialConditionText);

            m_lineGenerator.GenerateLines(m_definition.QuestionTexts);

            m_completionStateImg.sprite = GameDB.Instance.ReachStateIncomplete;
            m_completionState = ReachStateState.Incomplete;

            m_completionTimer = m_completionTime;
        }

        public override bool AnswerSelected() {
            return false;
        }

        public override void HandleEvaluation(bool correct) {
            if (m_evaluated) {
                // no need for duplicate evaluations
                return;
            }
        
            m_evaluated = correct;
        }

        public override bool IsCorrect() {
            return m_completionState == ReachStateState.Complete;
        }

        #endregion // IEvaluable

        private bool IsWithinRange()
        {
            for (int i = 0; i < m_definition.Targets.Count; i++)
            {
                SimStateTarget currTarget = m_definition.Targets[i];

                if (currTarget.TargetID == VarID.VolumeStop)
                {
                    Tuple<double, double> stopVals = World.Instance.get_stop_vals();
                    bool stop1OutOfRange = (stopVals.Item1 < currTarget.TargetVal - currTarget.TargetRange || stopVals.Item1 > currTarget.TargetVal + currTarget.TargetRange);
                    bool stop2OutOfRange = (stopVals.Item2 < currTarget.TargetVal - currTarget.TargetRange || stopVals.Item2 > currTarget.TargetVal + currTarget.TargetRange);
                    if (stop1OutOfRange && stop2OutOfRange)
                    {
                        // both stops are outside of limits
                        return false;
                    }
                }
                else
                {
                    double varVal = World.Instance.get_state_var(currTarget.TargetID);

                    if (varVal < currTarget.TargetVal - currTarget.TargetRange || varVal > currTarget.TargetVal + currTarget.TargetRange)
                    {
                        // is outside of limits
                        return false;
                    }
                }
            }

            // return true if no reqs were outside of limits
            return true;
        }

        private List<string> GetDiscrepancies()
        {
            List<string> discrepancies = new List<string>();

            for (int i = 0; i < m_definition.Targets.Count; i++)
            {
                SimStateTarget currTarget = m_definition.Targets[i];

                if (currTarget.TargetID == VarID.VolumeStop)
                {
                    Tuple<double, double> stopVals = World.Instance.get_stop_vals();
                    bool stop1OutOfRange = (stopVals.Item1 < currTarget.TargetVal - currTarget.TargetRange || stopVals.Item1 > currTarget.TargetVal + currTarget.TargetRange);
                    bool stop2OutOfRange = (stopVals.Item2 < currTarget.TargetVal - currTarget.TargetRange || stopVals.Item2 > currTarget.TargetVal + currTarget.TargetRange);
                    if (stop1OutOfRange && stop2OutOfRange)
                    {
                        // both stops are outside of limits
                        discrepancies.Add(currTarget.TargetID.ToString());
                    }
                }
                else
                {
                    double varVal = World.Instance.get_state_var(currTarget.TargetID);

                    if (varVal < currTarget.TargetVal - currTarget.TargetRange || varVal > currTarget.TargetVal + currTarget.TargetRange)
                    {
                        // is outside of limits
                        discrepancies.Add(currTarget.TargetID.ToString());
                    }
                }
            }

            return discrepancies;
        }
    }


}