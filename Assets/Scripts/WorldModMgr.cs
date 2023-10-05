using System.Collections;
using System.Collections.Generic;
using ThermoVR.Lab;
using ThermoVR.Tools;
using UnityEngine;

namespace ThermoVR {
    public class WorldModMgr : MonoBehaviour
    {
        private List<ToolType> ActiveTools = new List<ToolType>();
        private LimitsGroup ActiveLimits = new LimitsGroup();

        #region Tools

        public void SetAllowedTools(List<ToolType> allowedTools) {
            ActiveTools = allowedTools;
            GameMgr.Events.Dispatch(GameEvents.UpdateAllowedTools, allowedTools);
        }

        public void ResetToolRestrictions() {
            ActiveTools.Clear();
            GameMgr.Events.Dispatch(GameEvents.ResetToolRestrictions);
        }

        #endregion // Tools

        #region Limits

        public void SetLimits(LimitsGroup limits) {
            ActiveLimits = limits;
        }

        public void ResetLimits() {
            ActiveLimits = new LimitsGroup();
        }

        #endregion // Limits
    }
}