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

        [SerializeField] private GameObject m_linePrefab;
        public GameObject graph;

        [SerializeField] private int m_pressureSpacing;
        [SerializeField] private float m_volumeSpacing;
        [SerializeField] private int m_temperatureSpacing;

        [SerializeField] private Material m_pLineMat;
        [SerializeField] private Material m_vLineMat;
        [SerializeField] private Material m_tLineMat;

        [SerializeField] private Transform m_pressureLinesContainer;
        [SerializeField] private Transform m_volumeLinesContainer;
        [SerializeField] private Transform m_temperatureLinesContainer;

        #endregion // Inspector

        private List<LineRenderer> m_permanentLines;
        private List<LineRenderer> m_tempLines;

        private Routine m_initRoutine;
        private bool m_generationInProgress;
        private Routine m_generationRoutine;

        private void Start()
        {
            m_permanentLines = new List<LineRenderer>();
            m_tempLines = new List<LineRenderer>();

            m_initRoutine.Replace(InitRoutine());
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

        public void ClearTempLines()
        {
            for (int i = 0; i < m_tempLines.Count; i++)
            {
                Destroy(m_tempLines[i].gameObject);
            }
            m_tempLines.Clear();
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
                    if (pStepIndex % 6 == 0 && pStepIndex != 0)
                    {
                        spacingMult *= 4;
                    }
                    currPLine += m_pressureSpacing * spacingMult * pStepIndex;
                    m_permanentLines.Add(PopulatePLine(currPLine));
                    pStepIndex++;
                }
            }

            else
            {
                m_tempLines.Add(PopulatePLine(currPLine));
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
                    if (vStepIndex % 6 == 0 && vStepIndex != 0)
                    {
                        spacingMult *= 4;
                    }
                    currVLine += m_volumeSpacing * spacingMult * vStepIndex;
                    m_permanentLines.Add(PopulateVLine(currVLine));
                    vStepIndex++;
                }
            }

            else
            {
                m_tempLines.Add(PopulateVLine(currVLine));
            }

            m_generationInProgress = false;
            yield return null;
        }

        private IEnumerator GenerateTemperatureLinesRoutine(bool isPermanent)
        {
            m_generationInProgress = true;

            // Const T varies P
            int currTLine = 0;
            int spacingMult = 1;
            int tStepIndex = 0;

            if (isPermanent)
            {
                // for each interval
                while (currTLine < ThermoMath.t_max)
                {
                    currTLine += m_temperatureSpacing * spacingMult * tStepIndex;
                    m_permanentLines.Add(PopulateTLine(currTLine));
                    tStepIndex++;
                }
            }

            else
            {
                m_tempLines.Add(PopulateTLine(currTLine));
            }

            m_generationInProgress = false;
            yield return null;
        }

        #endregion // Routines

        #region Helpers

        private LineRenderer PopulatePLine(double constP)
        {
            LineRenderer newLine = Instantiate(m_linePrefab, m_pressureLinesContainer).GetComponent<LineRenderer>();
            newLine.material = m_pLineMat;
            newLine.transform.position = m_pressureLinesContainer.position;
            int currPosIndex = 0;
            Vector3 samplePos = Vector3.zero;

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
                    t += 2;
                }

                firstIter = false;

                v = ThermoMath.v_given_pt(p, t);

                samplePos = ThermoPresent.Instance.plot(p, v, t);

                if (p < ThermoMath.p_min || v < ThermoMath.v_min || t < ThermoMath.t_min
                    || p > ThermoMath.p_max || v > ThermoMath.v_max || t > ThermoMath.t_max)
                {
                    continue;
                }

                newLine.positionCount++;
                newLine.SetPosition(stepIndex, samplePos);
                stepIndex++;
            }

            return newLine;
        }

        private LineRenderer PopulateVLine(double constV)
        {
            LineRenderer newLine = Instantiate(m_linePrefab, m_volumeLinesContainer).GetComponent<LineRenderer>();
            newLine.material = m_vLineMat;
            newLine.transform.position = m_volumeLinesContainer.position;
            int currPosIndex = 0;
            Vector3 samplePos = Vector3.zero;

            #region First Half

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
                        t -= 2;
                    }
                    else
                    {
                        t -= 0.05;
                    }
                }

                firstIter = false;

                p = ThermoMath.p_given_vt(v, t);

                try
                {
                    var region = ThermoMath.region_given_pvt(p, v, t);
                    if (region <= ThermoMath.region_twophase)
                    {
                        continue;
                    }
                }
                catch
                {
                    /*
                    if (p < ThermoMath.psat_max)
                    {
                        continue;
                    }
                    */
                }

                samplePos = ThermoPresent.Instance.plot(p, v, t);

                if (Math.Abs(prev_y - samplePos.y) > 0.02 && prev_y != -1)
                {
                    continue;
                }

                if (p < ThermoMath.p_min || v < ThermoMath.v_min || t < ThermoMath.t_min
                    || p > ThermoMath.p_max || v > ThermoMath.v_max || t > ThermoMath.t_max)
                {
                    continue;
                }

                newLine.positionCount++;
                newLine.SetPosition(stepIndex, samplePos);
                stepIndex++;
                prev_y = samplePos.y;
            }

            #endregion // First Half

            #region Second half
            
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
                    if (iterIndex % 6 == 0)
                    {
                        spacingMult *= 4;
                    }

                    p += 1 * spacingMult;
                }
                iterIndex++;

                if (iterIndex > 100) { break; }

                if (in_two_phase)
                {
                    try
                    {
                        var x = ThermoMath.x_given_pv(p, v);
                        var h = ThermoMath.h_given_px(p, x);
                        t = ThermoMath.t_given_ph(p, h);
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

                if (samplePos.y - prev_y > 0.02 && prev_y != -1)
                {
                    // continue;
                }

                if (p < ThermoMath.p_min || v < ThermoMath.v_min || t < ThermoMath.t_min
                    || p > ThermoMath.p_max || v > ThermoMath.v_max || t > ThermoMath.t_max)
                {
                    continue;
                }

                newLine.positionCount++;
                new_positions.Add(samplePos);
                stepIndex++;
                prev_y = samplePos.y;
            }

            new_positions.Reverse();
            for (int i = 0; i < new_positions.Count; i++)
            {
                newLine.SetPosition(second_half_start_pos + i, new_positions[i]);
            }

            #endregion // Second Half

            return newLine;
        }

        private LineRenderer PopulateTLine(double constT)
        {
            LineRenderer newLine = Instantiate(m_linePrefab, m_temperatureLinesContainer).GetComponent<LineRenderer>();
            newLine.material = m_tLineMat;
            newLine.transform.position = m_temperatureLinesContainer.position;
            int currPosIndex = 0;
            Vector3 samplePos = Vector3.zero;

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
                    if (iterIndex % 6 == 0)
                    {
                        spacingMult *= 4;
                    }

                    p += 1 * spacingMult;
                }
                iterIndex++;

                if (iterIndex > 100) { break; }

                v = ThermoMath.v_given_pt(p, t);

                samplePos = ThermoPresent.Instance.plot(p, v, t);

                // straighten the line across 2-phase
                if (!jumpedTwoPhase && samplePos.x < 0.12 && newLine.positionCount > 0)
                {
                    jumpedTwoPhase = true;
                    var old = newLine.GetPosition(newLine.positionCount - 1);
                    samplePos.y = old.y;
                }

                if (p < ThermoMath.p_min || v < ThermoMath.v_min || t < ThermoMath.t_min
                    || p > ThermoMath.p_max || v > ThermoMath.v_max || t > ThermoMath.t_max)
                {
                    continue;
                }

                newLine.positionCount++;
                newLine.SetPosition(stepIndex, samplePos);
                stepIndex++;
            }

            return newLine;
        }

        #endregion // Helpers
    }
}