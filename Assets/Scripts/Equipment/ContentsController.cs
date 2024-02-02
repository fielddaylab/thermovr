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
        private void SetMaxScale() {
            m_MaxScale = this.transform.localScale;
        }
        /// <summary>
        /// Position the piston in the scene at the desired min position, then run this function.
        /// </summary>
        [ContextMenu("Set Default Scale")]
        private void SetDefaultScale()
        {
            m_DefaultScale = this.transform.localScale;
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

        /// <summary>
        /// Shortcut to preview min scale
        /// </summary>
        [ContextMenu("Preview Default Scale")]
        private void PreviewDefaultScale()
        {
            this.transform.localScale = m_DefaultScale;
        }

        #endregion // Editor

        #region Inspector

        public GameObject Water;
        public GameObject Steam;

        [SerializeField] private GameObject m_ContentsContainer;
        [SerializeField] private Vector3 m_MaxScale;
        [SerializeField] private Vector3 m_MinScale;
        [SerializeField] private Vector3 m_DefaultScale;

        private float m_Span; // difference between min and max positions

        #endregion // Inspector

        #region Unity Callbacks

        private void Awake() {
            m_Span = m_MaxScale.y - m_MinScale.y;
        }

        #endregion // Unity Callbacks

        #region Accessors

        public Vector3 GetMinScale() {
            return m_MinScale;
        }

        public Vector3 GetMaxScale() {
            return m_MaxScale;
        }

        public float GetSpan() {
            return m_Span;
        }

        #endregion // Accessors

    }
}