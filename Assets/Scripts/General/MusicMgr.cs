using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR
{
    public class MusicMgr : MonoBehaviour
    {
        [SerializeField] private AudioSource m_src;

        private float m_pauseTime = 0;

        private void OnEnable()
        {
            m_src.time = m_pauseTime;
            m_src.Play();
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