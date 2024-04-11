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

        private void Start()
        {
            m_numDivisions = m_dials.Length - 1;
            m_step = 1.0f / m_numDivisions;
        }

        public void ProcessVal(ref float toProcess)
        {
            float nearestSnap = Mathf.Round(toProcess / m_step) * m_step;

            if (Mathf.Abs(toProcess - nearestSnap) < 0.05f)
            {
                toProcess = nearestSnap;
            }
        }
    }
}