using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR
{
    public class PositionController : MonoBehaviour
    {
        #region Editor

        /// <summary>
        /// Position the piston in the scene at the desired min position, then run this function.
        /// </summary>
        [ContextMenu("Set Min Position")]
        protected void SetMinPos()
        {
            m_MinPos = this.transform.localPosition;
        }

        /// <summary>
        /// Position the piston in the scene at the desired min position, then run this function.
        /// </summary>
        [ContextMenu("Set Max Position")]
        protected void SetMaxPos()
        {
            m_MaxPos = this.transform.localPosition;
        }

        /// <summary>
        /// Shortcut to preview max position
        /// </summary>
        [ContextMenu("Preview Max Position")]
        protected void PreviewMaxPos()
        {
            this.transform.localPosition = m_MaxPos;
        }

        /// <summary>
        /// Shortcut to preview min position
        /// </summary>
        [ContextMenu("Preview Min Position")]
        protected void PreviewMinPos()
        {
            this.transform.localPosition = m_MinPos;
        }

        #endregion // Editor

        #region Inspector

        [SerializeField] private Vector3 m_MaxPos;
        [SerializeField] private Vector3 m_MinPos;

        protected float m_Span; // difference between min and max positions

        #endregion // Inspector

        #region Unity Callbacks

        private void Awake()
        {
            m_Span = m_MaxPos.y - m_MinPos.y;
        }

        #endregion // Unity Callbacks

        #region Accessors

        public Vector3 GetMinPos()
        {
            return m_MinPos;
        }

        public Vector3 GetMaxPos()
        {
            return m_MaxPos;
        }

        public float GetSpan()
        {
            return m_Span;
        }

        #endregion // Accessors

    }
}