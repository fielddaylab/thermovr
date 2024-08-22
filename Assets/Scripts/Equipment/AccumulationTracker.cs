using System.Collections;
using System.Collections.Generic;
using TMPro;
using UnityEngine;

namespace ThermoVR.Tools
{
    public class AccumulationTracker : MonoBehaviour
    {
        [SerializeField] private TextMeshPro m_accumulatedHeatEnergyText;

        private void Start()
        {
            EventMgr.Events.Register(GameEvents.AccumHeatEnergyUpdated, HandleAccumHeatEnergyUpdated);
        }

        private void HandleAccumHeatEnergyUpdated()
        {
            m_accumulatedHeatEnergyText.SetText(string.Format(DigitFormat.Heat, ToolMgr.Instance.GetAccumulatedHeatEnergy()));
        }
    }
}