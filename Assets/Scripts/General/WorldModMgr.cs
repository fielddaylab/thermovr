using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR.Lab;
using ThermoVR.State;
using ThermoVR.Tools;
using UnityEngine;

namespace ThermoVR
{
    public class WorldModMgr : MonoBehaviour
    {
        [SerializeField] private World m_world;

        private List<ToolType> m_activeTools = new List<ToolType>();
        private LimitsGroup m_activeLimits = new LimitsGroup(true);
        private bool m_limitsEnabled; // whether sim needs to check for limits
        private bool m_graphBallInteractable;

        private List<ToolType> m_allTools = new List<ToolType>();

        private void Awake() {
            foreach (ToolType tool in Enum.GetValues(typeof(ToolType))) {
                m_allTools.Add(tool);
            }

            m_limitsEnabled = false;

            m_graphBallInteractable = true;
        }

        #region Tools

        public void SetAllowedTools(List<ToolType> allowedTools) {
            if (allowedTools == null) {
                // No restrictions specified; default is all
                m_activeTools = new List<ToolType>(m_allTools);
            }
            else {
                m_activeTools = new List<ToolType>(allowedTools);
            }
            GameMgr.Events.Dispatch(GameEvents.UpdateAllowedTools, m_activeTools);
        }

        public void ResetToolRestrictions() {
            m_activeTools.Clear();
            m_activeTools = new List<ToolType>(m_allTools);
            GameMgr.Events?.Dispatch(GameEvents.ResetToolRestrictions);
        }

        public void EnableGraphBallInteractions()
        {
            m_graphBallInteractable = true;
        }

        public void DisableGraphBallInteractions()
        {
            m_graphBallInteractable = false;
        }

        public bool GraphBallInteractable()
        {
            return m_graphBallInteractable;
        }

        #endregion // Tools

        #region Limits

        public void SetLimits(LimitsGroup limits) {
            m_activeLimits = limits;
            m_limitsEnabled = true;
        }

        public void ResetLimits() {
            m_activeLimits = new LimitsGroup(true);
            m_limitsEnabled = false;
        }

        public bool LimitsEnabled() {
            return m_limitsEnabled;
        }

        public bool CrossedLimit(VarID simVar, double val, out bool ceiling) {
            double upper = -1;
            double lower = -1;

            switch (simVar) {
                case VarID.Pressure:
                    upper = m_activeLimits.Pressure.Ceiling;
                    lower = m_activeLimits.Pressure.Floor;
                    break;
                case VarID.Temperature:
                    upper = m_activeLimits.Temperature.Ceiling;
                    lower = m_activeLimits.Temperature.Floor;
                    break;
                case VarID.Volume:
                    upper = m_activeLimits.Volume.Ceiling;
                    lower = m_activeLimits.Volume.Floor;
                    break;
                case VarID.InternalEnergy:
                    upper = m_activeLimits.InternalEnergy.Ceiling;
                    lower = m_activeLimits.InternalEnergy.Floor;
                    break;
                case VarID.Entropy:
                    upper = m_activeLimits.Entropy.Ceiling;
                    lower = m_activeLimits.Entropy.Floor;
                    break;
                case VarID.Enthalpy:
                    upper = m_activeLimits.Enthalpy.Ceiling;
                    lower = m_activeLimits.Enthalpy.Floor;
                    break;
                case VarID.Quality:
                    upper = m_activeLimits.Quality.Ceiling;
                    lower = m_activeLimits.Quality.Floor;
                    break;
                default:
                    break;
            }

            if (upper != -1 && val > upper) {
                ceiling = true;
                return true;
            }
            else if (lower != -1 && val < lower) {
                ceiling = false;
                return true;
            }

            ceiling = false;
            return false;
        }

        #endregion // Limits
    }
}