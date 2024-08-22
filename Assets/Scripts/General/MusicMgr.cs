using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR
{
    public class MusicMgr : MonoBehaviour
    {
        [SerializeField] private AudioSource m_src;

        private bool m_initialized = false;
        private float m_pauseTime = 0;

        private void OnEnable()
        {
            if (PersistentState.Instance.Bools[PersistentVars.MusicOnStart])
            {
                if (m_initialized)
                {
                    m_src.time = m_pauseTime;
                    m_src.Play();
                }
                else
                {
                    m_src.Play();
                    m_initialized = true;
                }
            }
            else if (m_initialized)
            {
                m_src.Play();
            }
            else
            {
                m_initialized = true;
            }
        }

        private void Update()
        {
            m_pauseTime = m_src.time;
        }

        private void OnDisable()
        {
            m_src.Pause();
        }
    }

}