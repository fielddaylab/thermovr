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
        private bool m_loadComplete = false;

        private void Awake()
        {
            GameMgr.Events?.Register(GameEvents.InitialLoadComplete, HandleLoadComplete);
        }

        private void Update()
        {
            if (m_loadComplete && m_minTimer <= 0)
            {
                this.gameObject.SetActive(false);

                GameMgr.Events.Dispatch(GameEvents.TryNewName);
            }
            else
            {
                m_minTimer -= Time.deltaTime;

                m_loadIcon.transform.Rotate(-m_rotation * Time.deltaTime, Space.Self);
            }
        }

        #region Handlers

        private void HandleLoadComplete()
        {
            m_loadComplete = true;
        }

        #endregion // Handlers

    }
}