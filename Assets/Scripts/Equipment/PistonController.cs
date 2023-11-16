using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR {
    public class PistonController : MonoBehaviour
    {
        #region Editor

        /// <summary>
        /// Position the piston in the scene at the desired min position, then run this function.
        /// </summary>
        [ContextMenu("Set Min Position")]
        private void SetMinPos() {
            m_MinPos = this.transform.localPosition;
        }

        /// <summary>
        /// Position the piston in the scene at the desired min position, then run this function.
        /// </summary>
        [ContextMenu("Set Max Position")]
        private void SetMaxPos() {
            m_MaxPos = this.transform.localPosition;
        }

        /// <summary>
        /// Shortcut to preview max position
        /// </summary>
        [ContextMenu("Preview Max Position")]
        private void PreviewMaxPos() {
            this.transform.localPosition = m_MaxPos;
        }

        /// <summary>
        /// Shortcut to preview min position
        /// </summary>
        [ContextMenu("Preview Min Position")]
        private void PreviewMinPos() {
            this.transform.localPosition = m_MinPos;
        }

        #endregion // Editor

        #region Inspector

        [SerializeField] private Vector3 m_MaxPos;
        [SerializeField] private Vector3 m_MinPos;

        #endregion // Inspector

        #region Accessors

        public Vector3 GetMinPos() {
            return m_MinPos;
        }

        public Vector3 GetMaxPos() {
            return m_MaxPos;
        }

        #endregion // Accessors

    }
}