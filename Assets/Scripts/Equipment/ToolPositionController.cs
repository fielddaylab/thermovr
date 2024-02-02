using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR.Dials;
using ThermoVR.State;
using UnityEngine;

namespace ThermoVR.Tools {
    /// <summary>
    /// Controls the position of a tool based on volume values
    /// </summary>
    [RequireComponent(typeof(Tool))]
    public class ToolPositionController : MonoBehaviour
    {
        private enum MapType
        {
            Linear,
            Ln // natural log
        }

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
        [SerializeField] private float m_PistonCapHeight;
        [SerializeField] private float m_minOffsetVal; // offset from volume stop position at 0 volume to position at ThemoMath.v_min volume, which is the acting min position

        protected float m_Span; // difference between min and max positions

        [SerializeField] private MapType m_MappingType; // Note: when attached to volume stops, this type MUST match the volume calculation log base. Currently it is natural log.

        #endregion // Inspector

        private Tool m_Tool;

        #region Unity Callbacks 

        private void OnEnable()
        {
            m_Tool = this.GetComponent<Tool>();

            if (m_Tool)
            {
                m_Tool.ValUpdated.AddListener(HandleToolValUpdated);
            }

            m_Span = m_MaxPos.y - m_MinPos.y - m_PistonCapHeight;
        }

        private void OnDisable()
        {
            if (m_Tool)
            {
                m_Tool.ValUpdated.RemoveListener(HandleToolValUpdated);
            }
        }

        #endregion // Unity Callbacks

        #region Handlers

        private void HandleToolValUpdated()
        {
            if (!m_Tool)
            {
                return;
            }
            if (!m_Tool.engaged)
            {
                return;
            }

            float newToolVal = 0;

            if (m_MappingType == MapType.Ln)
            {
                // NOTE: Specifically for volume implementation
                newToolVal = m_Tool.GetVal();

                double totalHeight = Math.Log(newToolVal / ThermoState.piston_area) + ThermoState.log_offset_volume;

                float log_map = (float)(totalHeight / ThermoPresent.max_height_log); // map log height to range from 0 to 1

                float toolSpan = GetSpan(); // distance between min and max positions
                Vector3 new_tool_pos = m_Tool.transform.localPosition;
                new_tool_pos.y = GetMinPos().y + log_map * toolSpan;

                // if above sim volume (and thus above piston), add piston height
                if (newToolVal > ThermoState.Instance.volume)
                {
                    new_tool_pos.y += m_PistonCapHeight;
                }

                m_Tool.transform.localPosition = new_tool_pos;

                // Add dial val offset

            }
        }

        #endregion // Handlers

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