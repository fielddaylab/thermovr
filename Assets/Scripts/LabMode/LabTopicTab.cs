using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR.UI;
using UnityEngine;
using UnityEngine.UI;

namespace ThermoVR.Lab
{
    [RequireComponent(typeof(ThermoButton))]
    public class LabTopicTab : MonoBehaviour
    {
        public List<LabTab> TaskTabs;

        [SerializeField] private ThermoButton m_thermoButton;
        public Image ButtonImage;

        [SerializeField] private BoxCollider m_collider;

        public void EnableCollider()
        {
            m_collider.enabled = true;
        }

        public void DisableCollider()
        {
            m_collider.enabled = false;
        }

        public LabTopicTab()
        {
            TaskTabs = new List<LabTab>();
        }

        public ThermoButton Button
        {
            get { return m_thermoButton; }
        }
    }
}
