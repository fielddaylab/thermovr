using BeauUtil;
using BeauUtil.Extensions;
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

        [SerializeField] private bool m_isLabTask = true;

        private ReachStateDefinition m_definition;

        private ReachStateState m_completionState;

        private float m_completionTime;
        private float m_completionTimer;

        private bool m_completed;

        private bool m_initialized;
        private bool m_hasDef = false;

        private bool m_gameModeTargetReached = false;

        public void SetDefinition(ReachStateDefinition def) {
            m_definition = def;
            m_hasDef = true;

            m_completionTime = 4f;

            ResetState();
        }

        private void OnEnable()
        {
            if (!m_isLabTask)
            {
                EventMgr.Events.Register(GameEvents.GameModeTargetEntered, HandleGameModeTargetEntered);
                EventMgr.Events.Register(GameEvents.GameModeTargetExited, HandleGameModeTargetExited);
                EventMgr.Events.Register(GameEvents.GameModeExited, HandleGameModeExited);
                EventMgr.Events.Register(GameEvents.GameModeBeginGenerateTarget, HandleGameModeBeginGenerateTarget);
            }

            PlaceTargetZone();
            if (m_initialized)
            {
                if (m_isLabTask)
                {
                    EventMgr.Events?.Dispatch(GameEvents.TargetStateTaskBegan);
                }
            }
            m_initialized = true;
        }

        private void OnDisable()
        {
            if (m_isLabTask)
            {
                EventMgr.Events?.Dispatch(GameEvents.ClearTargetZone);
                EventMgr.Events?.Dispatch(GameEvents.TargetStateTaskEnded);
            }
            else
            {
                EventMgr.Events.Deregister(GameEvents.GameModeTargetEntered, HandleGameModeTargetEntered);
                EventMgr.Events.Deregister(GameEvents.GameModeTargetExited, HandleGameModeTargetExited);
                EventMgr.Events.Deregister(GameEvents.GameModeExited, HandleGameModeExited);
                EventMgr.Events.Deregister(GameEvents.GameModeBeginGenerateTarget, HandleGameModeBeginGenerateTarget);
            }
        }

        private void Update()
        {
            if (!m_hasDef) { return; }

            if (IsWithinRange() || m_gameModeTargetReached)
            {
                if (m_completionState == ReachStateState.Incomplete)
                {
                    // start timer
                    m_completionTimer = m_completionTime;
                    m_completionStateImg.fillAmount = 0;

                    m_completionState = ReachStateState.Countdown;
                    if (m_isLabTask)
                    {
                        EventMgr.Events.Dispatch(GameEvents.TargetStateEntered);
                    }
                    m_completed = false;
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
                    m_completed = false;
                }
                else
                {
                    // complete timer
                    if (!m_completed)
                    {
                        m_completionStateImg.sprite = GameDB.Instance.ReachStateComplete;
                        m_completionStateImg.fillAmount = 1;

                        if (m_isLabTask)
                        {
                            EventMgr.Events.Dispatch(GameEvents.TargetStateCompleted);
                        }
                        m_completionState = ReachStateState.Complete;
                    }
                    m_completed = true;
                }
            }
            else
            {
                m_completionStateImg.sprite = GameDB.Instance.ReachStateIncomplete;
                m_completionStateImg.fillAmount = 1;

                if (m_completionState != ReachStateState.Incomplete)
                {
                    if (m_isLabTask)
                    {
                        EventMgr.Events.Dispatch(GameEvents.TargetStateLost, EvtArgs.Ref(GetDiscrepancies()));
                    }

                    m_completionState = ReachStateState.Incomplete;
                }
                m_completed = false;
            }
        }

        private void PlaceTargetZone()
        {
            // try construct p, v, t
            double p, v, t;
            double pRange, vRange, tRange;
            p = v = t = -1;
            pRange = vRange = tRange = Mathf.Infinity;

            if (m_definition.Targets == null) { return; }

            for (int i = 0; i < m_definition.Targets.Count; i++)
            {
                SimStateTarget currTarget = m_definition.Targets[i];

                if (currTarget.TargetID == VarID.Pressure)
                {
                    // convert from kPa to Pa
                    p = currTarget.TargetVal * 1000;
                    pRange = currTarget.TargetRange * 1000;
                }
                else if (currTarget.TargetID == VarID.Volume)
                {
                    v = currTarget.TargetVal;
                    vRange = currTarget.TargetRange;
                }
                else if (currTarget.TargetID == VarID.Temperature)
                {
                    t = currTarget.TargetVal;
                    tRange = currTarget.TargetRange;
                }
            }

            if (p == -1 || v == -1 || t == -1)
            {
                // p and v, calc t
                if (p != -1 && v != -1)
                {
                    // TODO: improve this estimate
                    t = ThermoMath.iterate_t_given_pv(p, v, t);
                }

                // p and t, calc v
                else if (p != -1 && t != -1)
                {
                    v = ThermoMath.v_given_pt(p, t);
                }

                // v and t, calc p
                else if (v != -1 && t != -1)
                {
                    p = ThermoMath.p_given_vt(v, t);
                }
            }

            if (p < ThermoMath.p_min || v < ThermoMath.v_min || t < ThermoMath.t_min
                || p > ThermoMath.p_max || v > ThermoMath.v_max || t > ThermoMath.t_max)
            {
                Debug.Log("[ReachStateHub] No Target Zone generated. Insufficient valid dimension points.");
                EventMgr.Events.Dispatch(GameEvents.ClearTargetZone);
                return;
            }

            Vector3 targetZoneMaxPos = ThermoPresent.Instance.plot(
                Math.Clamp(p + pRange, ThermoMath.p_min, ThermoMath.p_max),
                Math.Clamp(v + vRange, ThermoMath.v_min, ThermoMath.v_max),
                Math.Clamp(t + tRange, ThermoMath.t_min, ThermoMath.t_max)
                );
            Vector3 targetZoneMinPos = ThermoPresent.Instance.plot(
                Math.Clamp(p - pRange, ThermoMath.p_min, ThermoMath.p_max),
                Math.Clamp(v - vRange, ThermoMath.v_min, ThermoMath.v_max),
                Math.Clamp(t - tRange, ThermoMath.t_min, ThermoMath.t_max)
                );
            Vector3 targetZoneCenterPos = (targetZoneMaxPos + targetZoneMinPos) / 2.0f;
            Vector3 targetZoneDims = new Vector3(
                    (float)(targetZoneMaxPos.x - targetZoneMinPos.x),
                    (float)(targetZoneMaxPos.y - targetZoneMinPos.y),
                    (float)(targetZoneMaxPos.z - targetZoneMinPos.z)
                );
            EventMgr.Events.Dispatch(GameEvents.TargetZoneUpdated, EvtArgs.Create(new STuple<Vector3, Vector3>(targetZoneCenterPos, targetZoneDims)));
        }

        #region IEvaluable

        public override void ResetState() {
            base.ResetState();

            m_initText.SetText(m_definition.InitialConditionText);

            m_lineGenerator.GenerateLines(m_definition.QuestionTexts);

            m_completionStateImg.sprite = GameDB.Instance.ReachStateIncomplete;
            m_completionState = ReachStateState.Incomplete;

            m_completionTimer = m_completionTime;

            PlaceTargetZone();
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
                    STuple<double, double> stopVals = World.Instance.get_stop_vals();
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
                    STuple<double, double> stopVals = World.Instance.get_stop_vals();
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

        #region Handlers

        private void HandleGameModeTargetEntered()
        {
            m_gameModeTargetReached = true;
            Debug.Log("[ReachState] Game mode reached");
        }

        private void HandleGameModeTargetExited()
        {
            m_gameModeTargetReached = false;
            Debug.Log("[ReachState] Game mode exited");
        }

        private void HandleGameModeExited()
        {
            m_gameModeTargetReached = false;
        }

        private void HandleGameModeBeginGenerateTarget()
        {
            m_gameModeTargetReached = false;
        }

        #endregion // Handlers
    }


}