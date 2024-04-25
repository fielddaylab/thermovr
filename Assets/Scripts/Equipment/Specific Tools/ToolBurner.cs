using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Tools
{
    public class ToolBurner : Tool
    {
        [Space(5)]
        [Header("Audio")]
        [SerializeField] private AudioSource m_audioSrc;

        #region Tool

        protected override void InitializeRoutines_Impl() {

        }

        protected override IEnumerator ActivationRoutine() {
            Debug.Log("[Triggers] Burner activated!");

            gameObject.SetActive(true);
            yield return null;
        }

        protected override IEnumerator DeactivationRoutine() {
            Debug.Log("[Triggers] Burner deactivated!");

            m_audioSrc.Stop();

            gameObject.SetActive(false);
            yield return null;
        }

        protected override IEnumerator BeginAdjustRoutine() {
            Debug.Log("[Triggers] Burner begin adjust!");
            yield return null;
        }

        protected override IEnumerator EndAdjustRoutine() {
            Debug.Log("[Triggers] Burner end adjust!");
            yield return null;
        }

        protected override IEnumerator EngageRoutine() {
            Debug.Log("[Triggers] Burner engaged!");
            if (!m_audioSrc.isPlaying)
            {
                m_audioSrc.Play();
            }

            yield return null;
        }

        protected override IEnumerator DisengageRoutine() {
            Debug.Log("[Triggers] Burner disengaged!");

            m_audioSrc.Stop();

            yield return null;
        }


        #endregion // Tool
    }
}