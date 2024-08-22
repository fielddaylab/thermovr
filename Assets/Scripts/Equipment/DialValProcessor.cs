using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Dials
{
    public class DialValProcessor : MonoBehaviour
    {
        [SerializeField] private Dial[] m_dials;

        private int m_numDivisions;
        private float m_step;

        private int m_mostRecentIndex;

        private void Start()
        {
            m_mostRecentIndex = 0;
            m_numDivisions = m_dials.Length - 1;
            m_step = 1.0f / m_numDivisions;
        }

        public void ProcessVal(ref float toProcess)
        {
            int nearestIndex = (int)Mathf.Round(toProcess / m_step);
            float nearestSnap = nearestIndex * m_step;

            if (Mathf.Abs(toProcess - nearestSnap) < 0.05f)
            {
                if (!ModeMgr.Instance.IsDesktop)
                {
                    toProcess = nearestSnap;
                }

                if (nearestIndex != m_mostRecentIndex)
                {
                    m_mostRecentIndex = nearestIndex;

                    // dial order is reversed
                    EventMgr.Events.Dispatch(GameEvents.SelectNewDial, m_numDivisions - nearestIndex);
                }
            }
        }
    }
}