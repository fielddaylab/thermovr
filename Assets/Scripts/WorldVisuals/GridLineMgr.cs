using BeauRoutine;
using Oculus.Interaction;
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
        [SerializeField] private int m_volumeSpacing;
        [SerializeField] private int m_temperatureSpacing;

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
                    m_permanentLines.Add(PopulatePLine(currPLine + m_pressureSpacing * spacingMult * pStepIndex));
                    pStepIndex++;
                    currPLine += m_pressureSpacing * spacingMult * pStepIndex;
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
            int currVLine = 0;


            m_generationInProgress = false;
            yield return null;
        }

        private IEnumerator GenerateTemperatureLinesRoutine(bool isPermanent)
        {
            m_generationInProgress = true;

            // Const T varies P
            int currTLine = 0;


            m_generationInProgress = false;
            yield return null;
        }

        #endregion // Routines

        #region Helpers

        private LineRenderer PopulatePLine(double constP)
        {
            LineRenderer newLine = Instantiate(m_linePrefab, m_pressureLinesContainer).GetComponent<LineRenderer>();
            newLine.transform.position = m_pressureLinesContainer.position;
            int currPosIndex = 0;
            Vector3 samplePos = Vector3.zero;

            double p;
            double v;
            double t = ThermoMath.t_min;
            int stepIndex = 0;

            // run down T line
            while (t < ThermoMath.t_max)
            {
                p = constP;
                t += 2;
                v = ThermoMath.v_given_pt(p, t);

                samplePos = ThermoPresent.Instance.plot(p, v, t);

                newLine.positionCount++;
                newLine.SetPosition(stepIndex, samplePos);
                stepIndex++;
            }

            return newLine;
        }

        private LineRenderer PopulateVLine()
        {
            LineRenderer newLine = Instantiate(m_linePrefab, m_volumeLinesContainer).GetComponent<LineRenderer>();

            int currPosIndex = 0;
            Vector3 samplePos = Vector3.zero;

            // run down T line
            {
                newLine.SetPosition(currPosIndex, samplePos);
            }

            return newLine;
        }

        private LineRenderer PopulateTLine()
        {
            LineRenderer newLine = Instantiate(m_linePrefab, m_temperatureLinesContainer).GetComponent<LineRenderer>();

            int currPosIndex = 0;
            Vector3 samplePos = Vector3.zero;

            // run down P line
            {
                newLine.SetPosition(currPosIndex, samplePos);
            }

            return newLine;
        }

        #endregion // Helpers
    }
}