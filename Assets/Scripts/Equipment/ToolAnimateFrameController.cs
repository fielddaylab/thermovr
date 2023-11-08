using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Tools
{
    [RequireComponent(typeof(Tool))]
    public class ToolAnimateFrameController : MonoBehaviour
    {
        #region Inspector

        [SerializeField] private Animator m_Animator; // the animator with the animation we're modifying
        [SerializeField] private int m_AnimationIndex; // the index of the specific animation clip we're working with

        [SerializeField] private float test_val;

        private AnimationClip m_Clip;
        private int m_ClipFrames;

        #endregion // Inspector

        private Tool m_Tool;

        #region Unity Callbacks

        private void OnEnable() {
            m_Tool = this.GetComponent<Tool>();

            if (m_Tool) {
                m_Tool.ValUpdated.AddListener(HandleToolValUpdated);
            }
            if (m_Animator) {
                m_Clip = m_Animator.runtimeAnimatorController.animationClips[m_AnimationIndex];
                // m_ClipFrames = m_Clip.length * m_Clip.frameRate;
            }
        }

        private void OnDisable() {
            if (m_Tool) {
                m_Tool.ValUpdated.RemoveListener(HandleToolValUpdated);
            }
        }

        #endregion // Unity Callbacks

        #region Handlers

        private void HandleToolValUpdated() {
            if (!m_Animator || !m_Tool) {
                return;
            }

            float newVal = m_Tool.GetVal();

            m_Animator.speed = 0;
            m_Animator.Play(m_Clip.name, 0, newVal);
        }

        #endregion // Handlers
    }
}