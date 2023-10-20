using BeauUtil;
using System.Collections;
using System.Collections.Generic;
using ThermoVR.State;
using TMPro;
using UnityEngine;

namespace ThermoVR.Tools
{
    public class UnitText : MonoBehaviour
    {
        [SerializeField] private ToolType m_ToolType;
        [SerializeField] private TMP_Text m_Text;


        private void Awake() {
            string toSet = "";
            switch (m_ToolType) {
                case ToolType.Burner:
                case ToolType.Coil:
                    toSet = Units.Heat;
                    break;
                case ToolType.Weight:
                case ToolType.NegativeWeight:
                    toSet = Units.Weight;
                    break;
                case ToolType.Insulator:
                    toSet = Units.Percent;
                    break;
                case ToolType.SurroundingTemperature:
                    toSet = Units.TemperatureK;
                    break;
                case ToolType.SurroundingPressure:
                    toSet = Units.AmbientPressure;
                    break;
                case ToolType.Stops:
                    toSet = Units.Volume;
                    break;
                default:
                    // No matching id found!
                    return;
            }

            m_Text.SetText(toSet);
        }
    }
}