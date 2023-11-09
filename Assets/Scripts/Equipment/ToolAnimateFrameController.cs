using System.Collections;
using System.Collections.Generic;
using ThermoVR.Dials;
using UnityEngine;

namespace ThermoVR.Tools
{
    [RequireComponent(typeof(Dial))]
    public class ToolAnimateFrameController : MonoBehaviour
    {
        #region Inspector

        [SerializeField] private Animator m_Animator; // the animator with the animation we're modifying
        [SerializeField] private int m_AnimationIndex; // the index of the specific animation clip we're working with

        [SerializeField] private float test_val;

        private AnimationClip m_Clip;

        #endregion // Inspector

        private Dial m_ToolDial;

        #region Unity Callbacks

        private void OnEnable() {
            m_ToolDial = this.GetComponent<Dial>();

            if (m_ToolDial) {
                m_ToolDial.DialMoved.AddListener(HandleToolValUpdated);
            }
            if (m_Animator) {
                m_Clip = m_Animator.runtimeAnimatorController.animationClips[m_AnimationIndex];
                // m_ClipFrames = m_Clip.length * m_Clip.frameRate;
            }
        }

        private void OnDisable() {
            if (m_ToolDial) {
                m_ToolDial.DialMoved.RemoveListener(HandleToolValUpdated);
            }
        }

        #endregion // Unity Callbacks

        #region Handlers

        private void HandleToolValUpdated() {
            if (!m_Animator || !m_ToolDial) {
                return;
            }

            float newVal = m_ToolDial.get_val();

            m_Animator.speed = 0;
            m_Animator.Play(m_Clip.name, 0, newVal);
        }

        #endregion // Handlers
    }
}