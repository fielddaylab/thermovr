using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR.Lab;
using ThermoVR.Tools;
using UnityEngine;

namespace ThermoVR {
    public class WorldModMgr : MonoBehaviour
    {
        private List<ToolType> m_activeTools = new List<ToolType>();
        private LimitsGroup m_activeLimits = new LimitsGroup();

        private List<ToolType> m_allTools = new List<ToolType>();

        private void Awake() {
            foreach (ToolType tool in Enum.GetValues(typeof(ToolType))) {
                m_allTools.Add(tool);
            }
        }

        #region Tools

        public void SetAllowedTools(List<ToolType> allowedTools) {
            if (allowedTools == null) {
                // No restrictions specified; default is all
                m_activeTools = new List<ToolType>(m_allTools);
            }
            else {
                m_activeTools = allowedTools;
            }
            GameMgr.Events.Dispatch(GameEvents.UpdateAllowedTools, m_activeTools);
        }

        public void ResetToolRestrictions() {
            m_activeTools.Clear();
            m_activeTools = new List<ToolType>(m_allTools);
            GameMgr.Events.Dispatch(GameEvents.ResetToolRestrictions);
        }

        #endregion // Tools

        #region Limits

        public void SetLimits(LimitsGroup limits) {
            m_activeLimits = limits;
        }

        public void ResetLimits() {
            m_activeLimits = new LimitsGroup();
        }

        #endregion // Limits
    }
}