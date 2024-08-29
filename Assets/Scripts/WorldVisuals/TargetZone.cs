using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR {
    public class TargetZone : MonoBehaviour
    {
        private bool m_initialized = false;
        [SerializeField] private GameObject m_checkObj; // the object that triggers collision events

        private void Start()
        {
            if (!m_initialized)
            {
                m_initialized = true;
            }
        }

        private void OnTriggerEnter(Collider col)
        {
            if (col.gameObject.name == m_checkObj.name)
            {
                EventMgr.Events.Dispatch(GameEvents.GameModeTargetEntered, this);
            }
        }

        private void OnTriggerExit(Collider col)
        {
            if (col.gameObject.name == m_checkObj.name)
            {
                EventMgr.Events.Dispatch(GameEvents.GameModeTargetExited, this);
            }
        }
    }
}
