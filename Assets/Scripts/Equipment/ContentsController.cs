using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR {
    public class ContentsController : MonoBehaviour
    {
        #region Editor

        /// <summary>
        /// Position the piston in the scene at the desired min position, then run this function.
        /// </summary>
        [ContextMenu("Set Min Scale")]
        private void SetMinScale() {
            m_MinScale = this.transform.localScale;
        }

        /// <summary>
        /// Position the piston in the scene at the desired min position, then run this function.
        /// </summary>
        [ContextMenu("Set Max Scale")]
        private void SetMaxSa() {
            m_MaxScale = this.transform.localScale;
        }

        /// <summary>
        /// Shortcut to preview max scale
        /// </summary>
        [ContextMenu("Preview Max Scale")]
        private void PreviewMaxScale() {
            this.transform.localScale = m_MaxScale;
        }

        /// <summary>
        /// Shortcut to preview min scale
        /// </summary>
        [ContextMenu("Preview Min Scale")]
        private void PreviewMinScale() {
            this.transform.localScale = m_MinScale;
        }

        #endregion // Editor

        #region Inspector

        public GameObject Water;
        public GameObject Steam;

        [SerializeField] private GameObject m_ContentsContainer;
        [SerializeField] private Vector3 m_MaxScale;
        [SerializeField] private Vector3 m_MinScale;

        #endregion // Inspector

        #region Accessors

        public Vector3 GetMinScale() {
            return m_MinScale;
        }

        public Vector3 GetMaxScale() {
            return m_MaxScale;
        }

        #endregion // Accessors

    }
}