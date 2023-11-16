using BeauRoutine;
using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR
{
    /// <summary>
    /// Coordinates state between two actors in the world
    /// </summary>
    public class CoordinatedState : MonoBehaviour
    {
        [SerializeField] private GameObject m_SharedTarget; // the shared target of anyhting coordinating through this state

        private float m_InitialVal; // the initial value of the coordinated value
        private float m_Val; // a float value to coordinate
        [SerializeField] private float m_MaxValDelta; // the most this float can change

        private Routine m_CoordinatedRoutine;

        public GameObject GetSharedTarget() {
            return m_SharedTarget;
        }

        public void SetInitialVal(float initial) {
            m_InitialVal = initial;
        }

        public float GetInitialVal() {
            return m_InitialVal;
        }

        public bool TryModifyVal(float delta) {
            if (Mathf.Abs(m_Val + delta) <= m_MaxValDelta) {
                m_Val += delta;
                return true;
            }
            return false;
        }

        public float GetVal() {
            return m_Val;
        }

        public float GetMaxValDelta() {
            return m_MaxValDelta;
        }

        public void ReplaceCoordinatedRoutine(IEnumerator routine) {
            m_CoordinatedRoutine.Replace(routine);
        }
    }
}

