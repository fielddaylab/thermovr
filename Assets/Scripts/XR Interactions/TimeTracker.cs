using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Controls
{
    public class TimeTracker : MonoBehaviour
    {
        private float m_elapsedTime;

        private void OnEnable()
        {
            m_elapsedTime = 0;
        }

        private void Update()
        {
            m_elapsedTime += Time.deltaTime;

            GameMgr.Events.Dispatch(GameEvents.ElapsedTimeUpdated, m_elapsedTime);
        }
    }
}
