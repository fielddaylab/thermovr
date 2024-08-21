using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

namespace ThermoVR.UI
{
    [DefaultExecutionOrder(5)]
    public class UILoading : MonoBehaviour
    {
        [SerializeField] private Image m_loadIcon;
        [SerializeField] private Vector3 m_rotation;

        [SerializeField] private float m_minTimer = 3;

        private void Update()
        {
            if (m_minTimer > 0)
            {
                m_minTimer -= Time.deltaTime;
            }

            m_loadIcon.transform.Rotate(-m_rotation * Time.deltaTime, Space.Self);
        }

        public bool MinLoadTimeCompleted()
        {
            return m_minTimer <= 0;
        }
    }
}