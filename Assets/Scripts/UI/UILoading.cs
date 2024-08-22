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

        [Header("Elements")]
        [SerializeField] private GameObject m_desktopGroup;
        [SerializeField] private GameObject m_vrGroup;
        [SerializeField] private RectTransform m_loadGroup;

        [Header("VR layout")]
        [SerializeField] private Vector3 m_vrLoadPos;

        private void Awake()
        {
            if (ModeMgr.Instance.IsDesktop)
            {
                m_desktopGroup.SetActive(true);
                m_vrGroup.SetActive(false);
            }
            else
            {
                m_desktopGroup.SetActive(false);
                m_vrGroup.SetActive(true);

                m_loadGroup.anchoredPosition = m_vrLoadPos;
            }
        }

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