using System.Collections;
using System.Collections.Generic;
using TMPro;
using UnityEngine;

namespace ThermoVR.Tools
{
    public enum AccumulationType { 
        HEAT,
        WEIGHT
    }

    public class AccumulationTracker : MonoBehaviour
    {
        [SerializeField] private TextMeshPro m_accumulatedEnergyText;
        [SerializeField] private AccumulationType m_accumulationType;

        private void Start()
        {
            switch (m_accumulationType) {
                case AccumulationType.HEAT:
                    EventMgr.Events.Register(GameEvents.AccumHeatEnergyUpdated, HandleAccumHeatEnergyUpdated);
                    break;
                case AccumulationType.WEIGHT:
                    EventMgr.Events.Register(GameEvents.AccumWeightEnergyUpdated, HandleAccumWeightEnergyUpdated);
                    break;
                default:
                    break;
            }

        }

        private void HandleAccumHeatEnergyUpdated()
        {
            m_accumulatedEnergyText.SetText(string.Format(DigitFormat.Heat, ToolMgr.Instance.GetAccumulatedHeatEnergy()));
        }

        private void HandleAccumWeightEnergyUpdated()
        {
            m_accumulatedEnergyText.SetText(string.Format(DigitFormat.Weight, ToolMgr.Instance.GetAccumulatedWeightEnergy()));
        }
    }
}