using BeauRoutine;
using Oculus.Interaction;
using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR.State;
using UnityEngine;
using UnityEngine.UIElements;
using static IF97.Region3Backwards;

namespace ThermoVR
{
    public class GridLineMgr : MonoBehaviour
    {
        #region Inspector

        [SerializeField] private bool m_generateOnStart;

        [SerializeField] private GameObject m_linePrefab;
        public GameObject graph;

        [SerializeField] private Material m_pLineMat;
        [SerializeField] private Material m_vLineMat;
        [SerializeField] private Material m_tLineMat;

        [SerializeField] private Transform m_pressureLinesContainer;
        [SerializeField] private Transform m_volumeLinesContainer;
        [SerializeField] private Transform m_temperatureLinesContainer;

        [Space(5)]
        [Header("Intervals")]

        [Header("Origins")]
        [SerializeField] private int m_pressureOriginSpacing;
        [SerializeField] private float m_volumeOriginSpacing;
        [SerializeField] private int m_temperatureOriginSpacing;

        [Header("Steps (High Fidelity)")]
        [SerializeField] private float m_pStepSpacing;
        [SerializeField] private float m_vStepSpacing;
        [SerializeField] private float m_vStepDomeSpacing;
        [SerializeField] private float m_tStepSpacing;

        [Header("Steps (Low Fidelity)")]
        [SerializeField] private float m_pStepSpacingLow;
        [SerializeField] private float m_vStepSpacingLow;
        [SerializeField] private float m_vStepDomeSpacingLow;
        [SerializeField] private float m_tStepSpacingLow;


        [Header("Spacing Multiplier")]
        [SerializeField] private int m_stepsBetweenSpacing;
        [SerializeField] private int m_spacingMult;

        [Header("Significance Thresholds")]
        [SerializeField] private float m_sigYThreshold;

        [Header("Processing")]
        [SerializeField] private float m_processPThreshold;
        [SerializeField] private float m_processVThreshold;
        [SerializeField] private float m_processTThreshold;

        [Header("Live Generation")]
        [SerializeField] private GameObject m_graphBall;
        [SerializeField] private LineRenderer[] m_tempLines;
        [SerializeField] private bool m_realTimeGeneration;
        [SerializeField] private float m_refreshThreshold;

        private Vector3 m_lastKnownPos;

        // [SerializeField] private float m_volumeOriginSpacing;


        #endregion // Inspector

        private List<LineRenderer> m_permanentLines;

        private Routine m_initRoutine;
        private bool m_generationInProgress;
        private Routine m_generationRoutine;

        private void Start()
        {
            if (m_generateOnStart)
            {
                m_permanentLines = new List<LineRenderer>();

                m_initRoutine.Replace(InitRoutine());
            }

            m_lastKnownPos = m_graphBall.transform.position;
        }

        private void Update()
        {
            if (m_realTimeGeneration)
            {
                var dist = Vector3.Distance(m_graphBall.transform.position, m_lastKnownPos);
                if (dist > m_refreshThreshold)
                {
                    m_tempLines[0].positionCount = 0;
                    PopulatePLine(World.Instance.get_state_var(VarID.Pressure) * 1000, ref m_tempLines[0], false);
                    m_tempLines[1].positionCount = 0;
                    PopulateVLine(World.Instance.get_state_var(VarID.Volume), ref m_tempLines[1], false);
                    m_tempLines[2].positionCount = 0;
                    PopulateTLine(World.Instance.get_state_var(VarID.Temperature), ref m_tempLines[2], false);

                    m_lastKnownPos = m_graphBall.transform.position;
                }
            }
        }

        public void GenerateLines(VarID axis, bool isPermanent)
        {
            // can only generate lines acros valid dimensions
            if (axis != VarID.Pressure && axis != VarID.Temperature && axis != VarID.Volume)
            {
                return;
            }

            if (m_generationInProgress) { return; }


            if (axis == VarID.Pressure)
            {
                m_generationRoutine.Replace(GeneratePressureLinesRoutine(isPermanent));
            }
            else if (axis == VarID.Volume)
            {
                m_generationRoutine.Replace(GenerateVolumeLinesRoutine(isPermanent));
            }
            else if (axis == VarID.Temperature)
            {
                m_generationRoutine.Replace(GenerateVolumeLinesRoutine(isPermanent));
            }
        }


        public void ShowPermanentLines()
        {
            foreach (var line in m_permanentLines)
            {
                line.enabled = true;
            }
        }

        public void HidePermanentLines()
        {
            foreach (var line in m_permanentLines)
            {
                line.enabled = false;
            }
        }

        public void ShowTempLines()
        {
            foreach (var line in m_tempLines)
            {
                line.enabled = true;
            }
        }

        public void HideTempLines()
        {
            foreach (var line in m_tempLines)
            {
                line.enabled = false;
            }
        }

        #region Routines

        private IEnumerator InitRoutine()
        {
            // generate permanent lines for p, v, and t
            yield return GeneratePressureLinesRoutine(true);
            yield return GenerateVolumeLinesRoutine(true);
            yield return GenerateTemperatureLinesRoutine(true);
        }

        private IEnumerator GeneratePressureLinesRoutine(bool isPermanent)
        {
            m_generationInProgress = true;

            // const P varies T
            int currPLine = (int)ThermoMath.p_min;
            int spacingMult = 1;
            int pStepIndex = 0;

            if (isPermanent)
            {
                // for each interval
                while (currPLine < ThermoMath.p_max)
                {
                    if (pStepIndex % m_stepsBetweenSpacing == 0 && pStepIndex != 0)
                    {
                        spacingMult *= m_spacingMult;
                    }
                    currPLine += m_pressureOriginSpacing * spacingMult * pStepIndex;
                    LineRenderer newLine = Instantiate(m_linePrefab, m_pressureLinesContainer).GetComponent<LineRenderer>();
                    PopulatePLine(currPLine, ref newLine);
                    m_permanentLines.Add(newLine);
                    pStepIndex++;
                }
            }
            else
            {
                /*
                LineRenderer newLine = Instantiate(m_linePrefab, m_pressureLinesContainer).GetComponent<LineRenderer>();
                PopulatePLine(currPLine, ref newLine);
                m_tempLines.Add(newLine);
                */
            }

            m_generationInProgress = false;
            yield return null;
        }

        private IEnumerator GenerateVolumeLinesRoutine(bool isPermanent)
        {
            m_generationInProgress = true;

            // Const V varies T
            double currVLine = ThermoMath.v_effective_min + 0.001;
            int spacingMult = 1;
            int vStepIndex = 0;

            if (isPermanent)
            {
                // for each interval
                while (currVLine < ThermoMath.v_max)
                {
                    if (vStepIndex % m_stepsBetweenSpacing == 0 && vStepIndex != 0)
                    {
                        spacingMult *= m_spacingMult;
                    }
                    currVLine += m_volumeOriginSpacing * spacingMult * vStepIndex;
                    LineRenderer newLine = Instantiate(m_linePrefab, m_volumeLinesContainer).GetComponent<LineRenderer>();
                    PopulateVLine(currVLine, ref newLine);
                    m_permanentLines.Add(newLine);
                    vStepIndex++;
                }
            }
            else
            {
                // m_tempLines.Add(PopulateVLine(currVLine));
            }

            m_generationInProgress = false;
            yield return null;
        }

        private IEnumerator GenerateTemperatureLinesRoutine(bool isPermanent)
        {
            m_generationInProgress = true;

            // Const T varies P
            int currTLine = (int)ThermoMath.t_min;
            int spacingMult = 1;
            int tStepIndex = 0;

            if (isPermanent)
            {
                // for each interval
                while (currTLine < ThermoMath.t_max)
                {
                    currTLine += m_temperatureOriginSpacing * spacingMult * tStepIndex;
                    LineRenderer newLine = Instantiate(m_linePrefab, m_temperatureLinesContainer).GetComponent<LineRenderer>();
                    PopulateTLine(currTLine, ref newLine);
                    m_permanentLines.Add(newLine);
                    tStepIndex++;
                }
            }
            else
            {
                // m_tempLines.Add(PopulateTLine(currTLine));
            }

            m_generationInProgress = false;
            yield return null;
        }

        #endregion // Routines

        #region Helpers

        private void PopulatePLine(double constP, ref LineRenderer newLine, bool highResolution = true)
        {
            newLine.material = m_pLineMat;
            newLine.transform.position = m_pressureLinesContainer.position;
            int currPosIndex = 0;
            Vector3 samplePos = Vector3.zero;
            List<Vector3> allPositions = new List<Vector3>();

            double p = constP;
            double v;
            double t = ThermoMath.t_min;
            int stepIndex = 0;
            bool firstIter = true;
            
            // run down T line
            while (t < ThermoMath.t_max)
            {
                if (!firstIter)
                {
                    t += highResolution ? m_pStepSpacing : m_pStepSpacingLow;
                }

                firstIter = false;

                try
                {
                    v = ThermoMath.v_given_pt(p, t, bypass_error: true);
                }
                catch { continue; }

                samplePos = ThermoPresent.Instance.plot(p, v, t);

                if (p < ThermoMath.p_min || v < ThermoMath.v_min || t < ThermoMath.t_min
                    || p > ThermoMath.p_max || v > ThermoMath.v_max || t > ThermoMath.t_max)
                {
                    continue;
                }

                allPositions.Add(samplePos);
                stepIndex++;
            }

            ProcessPLine(ref newLine, allPositions);
        }

        private void PopulateVLine(double constV, ref LineRenderer newLine, bool highResolution = true)
        {
            newLine.material = m_vLineMat;
            newLine.transform.position = m_volumeLinesContainer.position;
            int currPosIndex = 0;
            Vector3 samplePos = Vector3.zero;
            List<Vector3> allPositions = new List<Vector3>();

            #region First Half (Vapor)

            double p;
            double v = constV;
            double t = ThermoMath.t_max;
            int stepIndex = 0;
            bool firstIter = true;
            double prev_y = -1;

            // run down T line
            while (t > ThermoMath.t_min)
            {
                if (!firstIter)
                {
                    if (t > ThermoMath.t_crit)
                    {
                        t -= highResolution ? m_vStepSpacing : m_vStepSpacingLow;
                    }
                    else
                    {
                        t -= highResolution ? m_vStepDomeSpacing : m_vStepDomeSpacingLow;
                    }
                }

                firstIter = false;

                p = ThermoMath.p_given_vt(v, t, bypass_error: true);

                try
                {
                    var region = ThermoMath.region_given_pvt(p, v, t);
                    if (region <= ThermoMath.region_twophase)
                    {
                        continue;
                    }
                }
                catch { }

                samplePos = ThermoPresent.Instance.plot(p, v, t);

                if (Math.Abs(prev_y - samplePos.y) > m_sigYThreshold && prev_y != -1)
                {
                    continue;
                }

                if (p < ThermoMath.p_min || v < ThermoMath.v_min || t < ThermoMath.t_min
                    || p > ThermoMath.p_max || v > ThermoMath.v_max || t > ThermoMath.t_max)
                {
                    continue;
                }

                allPositions.Add(samplePos);
                stepIndex++;
                prev_y = samplePos.y;
            }

            #endregion // First Half (Vapor)

            #region Second half (Two-Phase)

            currPosIndex = 0;
            samplePos = Vector3.zero;

            p = ThermoMath.p_min;
            v = constV;
            firstIter = true;
            int iterIndex = 0;
            int spacingMult = 1;
            prev_y = -1;
            bool in_two_phase = true;

            int second_half_start_pos = stepIndex;
            List<Vector3> new_positions = new List<Vector3>();

            // run down P line
            while (p < ThermoMath.p_max)
            {
                if (iterIndex != 0)
                {
                    if (iterIndex % m_stepsBetweenSpacing == 0)
                    {
                        spacingMult *= m_spacingMult;
                    }

                    p += highResolution ? m_vStepDomeSpacing * spacingMult : m_vStepDomeSpacingLow * spacingMult;
                }
                iterIndex++;

                if (iterIndex > 100) { break; }

                if (in_two_phase)
                {
                    try
                    {
                        var x = ThermoMath.x_given_pv(p, v, bypass_error:true);
                        var h = ThermoMath.h_given_px(p, x);
                        t = ThermoMath.t_given_ph(p, h, bypass_error: true);
                    }
                    catch
                    {
                        if (p > ThermoMath.p_min + 500) {
                            in_two_phase = false;
                        }
                        continue;
                    }
                }
                else
                {
                    continue;
                }

                samplePos = ThermoPresent.Instance.plot(p, v, t);

                if (samplePos.y - prev_y > m_sigYThreshold && prev_y != -1)
                {
                    // continue;
                }

                if (p < ThermoMath.p_min || v < ThermoMath.v_min || t < ThermoMath.t_min
                    || p > ThermoMath.p_max || v > ThermoMath.v_max || t > ThermoMath.t_max)
                {
                    continue;
                }

                // newLine.positionCount++;
                new_positions.Add(samplePos);
                stepIndex++;
                prev_y = samplePos.y;
            }

            new_positions.Reverse();
            for (int i = 0; i < new_positions.Count; i++)
            {
                allPositions.Add(new_positions[i]);
            }

            #endregion // Second Half (Two-Phase)

            ProcessVLine(ref newLine, allPositions);
        }

        private void PopulateTLine(double constT, ref LineRenderer newLine, bool highResolution = true)
        {
            newLine.material = m_tLineMat;
            newLine.transform.position = m_temperatureLinesContainer.position;
            int currPosIndex = 0;
            Vector3 samplePos = Vector3.zero;
            List<Vector3> allPositions = new List<Vector3>();

            double p = ThermoMath.p_min;
            double v;
            double t = constT;
            int stepIndex = 0;
            bool firstIter = true;
            int iterIndex = 0;
            int spacingMult = 1;
            bool jumpedTwoPhase = false;

            // run down P line
            while (t < ThermoMath.t_max)
            {
                if (iterIndex != 0)
                {
                    if (iterIndex % m_stepsBetweenSpacing == 0)
                    {
                        spacingMult *= m_spacingMult;
                    }

                    p += highResolution ? m_tStepSpacing * spacingMult : m_tStepSpacingLow * spacingMult;
                }
                iterIndex++;

                if (iterIndex > 100) { break; }

                try
                {
                    v = ThermoMath.v_given_pt(p, t, bypass_error: true);
                }
                catch { continue; }

                samplePos = ThermoPresent.Instance.plot(p, v, t);

                // straighten the line across 2-phase
                if (!jumpedTwoPhase && samplePos.x < 0.15 && allPositions.Count > 0 && t < ThermoMath.t_crit)
                {
                    jumpedTwoPhase = true;
                    var old = allPositions[allPositions.Count - 1];
                    samplePos.y = old.y;
                }

                if (p < ThermoMath.p_min || v < ThermoMath.v_min || t < ThermoMath.t_min
                    || p > ThermoMath.p_max || v > ThermoMath.v_max || t > ThermoMath.t_max)
                {
                    continue;
                }

                allPositions.Add(samplePos);
                stepIndex++;
            }

            ProcessTLine(ref newLine, allPositions);
        }

        private void ProcessPLine(ref LineRenderer line, List<Vector3> allPositions)
        {
            // look for significant change in x
            float lastKnownX = 0;
            int numPruned = 0;
            for (int i = 0; i < allPositions.Count; i++)
            {
                if (i == 0)
                {
                    // first index always stays
                    lastKnownX = allPositions[i].x;
                    continue;
                }
                else if (i == allPositions.Count - 1)
                {
                    // last index always stays
                    continue;
                }
                else
                {
                    // prune if not significant change in x && next point is not across two-phase
                    if (Math.Abs(allPositions[i].x - lastKnownX) < m_processPThreshold && Math.Abs(allPositions[i + 1].x - allPositions[i].x) < 0.1)
                    {
                        allPositions.RemoveAt(i);
                        i--;
                        numPruned++;
                    }
                    else
                    {
                        lastKnownX = allPositions[i].x;
                    }
                }
            }

            Debug.Log("[Line Processing] Pruned " + numPruned + " points from P line.");

            line.positionCount = allPositions.Count;
            for (int i = 0; i < allPositions.Count; i++)
            {
                line.SetPosition(i, allPositions[i]);
            }
        }

        private void ProcessVLine(ref LineRenderer line, List<Vector3> allPositions)
        {
            // look for significant change in y
            float lastKnownY = 0;
            int numPruned = 0;
            for (int i = 0; i < allPositions.Count; i++)
            {
                if (i == 0)
                {
                    // first index always stays
                    lastKnownY = allPositions[i].y;
                    continue;
                }
                else if (i == allPositions.Count - 1)
                {
                    // last index always stays
                    continue;
                }
                else
                {
                    // prune if not significant change in y
                    if (Math.Abs(allPositions[i].y - lastKnownY) < m_processVThreshold
                        //|| allPositions[i].y < 0.001
                        )
                    {
                        allPositions.RemoveAt(i);
                        i--;
                        numPruned++;
                    }
                    else
                    {
                        lastKnownY = allPositions[i].y;
                    }
                }
            }

            Debug.Log("[Line Processing] Pruned " + numPruned + " points from V line.");

            // Add positions to line
            line.positionCount = allPositions.Count;
            for (int i = 0; i < allPositions.Count; i++)
            {
                line.SetPosition(i, allPositions[i]);
            }
        }

        private void ProcessTLine(ref LineRenderer line, List<Vector3> allPositions)
        {
            // look for significant change in x
            float lastKnownX = 0;
            int numPruned = 0;
            for (int i = 0; i < allPositions.Count; i++)
            {
                if (i == 0)
                {
                    // first index always stays
                    lastKnownX = allPositions[i].x;
                    continue;
                }
                else if (i == allPositions.Count - 1)
                {
                    // last index always stays
                    continue;
                }
                else
                {
                    // prune if not significant change in x && next point is not across two-phase
                    if (Math.Abs(allPositions[i].x - lastKnownX) < m_processTThreshold && Math.Abs(allPositions[i + 1].x - allPositions[i].x) < 0.05)
                    {
                        allPositions.RemoveAt(i);
                        i--;
                        numPruned++;
                    }
                    else
                    {
                        lastKnownX = allPositions[i].x;
                    }
                }
            }

            Debug.Log("[Line Processing] Pruned " + numPruned + " points from T line.");


            line.positionCount = allPositions.Count;
            for (int i = 0; i < allPositions.Count; i++)
            {
                line.SetPosition(i, allPositions[i]);
            }
        }

        #endregion // Helpers
    }
}