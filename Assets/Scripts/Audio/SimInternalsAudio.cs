using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Audio
{
    public class SimInternalsAudio : MonoBehaviour
    {
        private enum PressurePlaybackState
        {
            Stopped,
            Playing
        }

        [SerializeField] private AudioClip m_pressureDropClip;
        [SerializeField] private AudioSource m_pressureAudioSrc;

        [SerializeField] private float m_pressureDropThreshold;

        [SerializeField] private float m_minPlayTime;

        private float m_minPlayTimer;

        private PressurePlaybackState m_pressurePlaybackState;


        private void Update()
        {
            if (!m_pressureAudioSrc.isPlaying && m_pressurePlaybackState == PressurePlaybackState.Playing)
            {
                m_pressurePlaybackState = PressurePlaybackState.Stopped;
            }

            if (m_minPlayTimer > 0)
            {
                m_minPlayTimer -= Time.deltaTime;
            }
        }

        public void ProcessPressureAudio(float delta_p)
        {
            if (Mathf.Abs(delta_p) >= m_pressureDropThreshold && Mathf.Sign(delta_p) < 0)
            {
                if (m_pressurePlaybackState == PressurePlaybackState.Playing)
                {
                    // continue playing audio
                    return;
                }

                // turn on audio
                m_pressureAudioSrc.Stop();
                m_pressureAudioSrc.clip = m_pressureDropClip;
                m_pressureAudioSrc.Play();
                // m_pressureAudioSrc.PlayOneShot(m_pressureDropClip);

                m_minPlayTimer = m_minPlayTime;

                m_pressurePlaybackState = PressurePlaybackState.Playing;
            }
            else if (m_minPlayTimer <= 0)
            {
                // turn off audio if playing
                m_pressureAudioSrc.Stop();
                m_pressurePlaybackState = PressurePlaybackState.Stopped;
            }
        }
    }
}