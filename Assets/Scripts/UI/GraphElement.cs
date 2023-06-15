using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.UI.GraphElements
{
    public enum GraphElementID : byte
    {
        AxisNumbers,
        RegionLabels,
        GridLines,
        ConstantLines
    }

    public struct GraphSettingUpdate
    {
        public GraphElementID GraphElementID;
        public bool ToggleVal;

        public GraphSettingUpdate(GraphElementID id, bool toggleVal) {
            GraphElementID = id;
            ToggleVal = toggleVal;
        }
    }


    public class GraphElement : MonoBehaviour
    {
        [SerializeField] private GraphElementID m_elementID;
        [SerializeField] private bool m_startVisible = false;

        private void Awake() {
            // start invisible
            this.gameObject.SetActive(m_startVisible);
            GameMgr.Events?.Register<GraphSettingUpdate>(GameEvents.UpdateGraphSetting, HandleUpdateGraphSetting);
        }

        #region Handlers

        private void HandleUpdateGraphSetting(GraphSettingUpdate update) {
            if (update.GraphElementID == m_elementID) {
                SetActive(update.ToggleVal);
            }
        }

        #endregion // Handlers

        #region Helpers

        private void SetActive(bool active) {
            this.gameObject.SetActive(active);
        }

        #endregion // Helpers
    }
}