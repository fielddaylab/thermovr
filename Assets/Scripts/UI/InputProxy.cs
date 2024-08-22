using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR.Dials;
using ThermoVR.Tools;
using TMPro;
using UnityEngine;

namespace ThermoVR
{
    public class InputProxy : MonoBehaviour
    {
        [SerializeField] private TMP_Text m_textToChange;
        [SerializeField] private Dial m_dialToSet;
        [SerializeField] private float m_mult = 1;

        public string InputUnits;

        private void Start()
        {
            if (m_dialToSet)
            {
                InputUnits = m_dialToSet.get_relevant_tools()[0].display_unit;
            }
        }

        public ToolType ToolType()
        {
            return m_dialToSet.get_relevant_tools()[0].tool_type;
        }

        public void SetValue(string newValStr)
        {
            float newVal;

            try
            {
                // simple input
                newVal = float.Parse(newValStr);
            }
            catch 
            { 
                try
                {
                    decimal d = Decimal.Parse(newValStr, System.Globalization.NumberStyles.Float);
                    newVal = (float)(d);
                }
                catch
                {
                    EventMgr.Events.Dispatch(GameEvents.SetInvalidToolVal, newValStr);
                    return;
                }
            }

            newVal *= m_mult;

            if (m_dialToSet)
            {
                if (m_dialToSet.val_within_range(newVal))
                {
                    m_dialToSet.convert_and_set_map(newVal);
                    EventMgr.Events.Dispatch(GameEvents.ProxyInputSubmitted, newVal);
                }
                else
                {
                    EventMgr.Events.Dispatch(GameEvents.SetInvalidToolVal, newValStr);
                }
            }
            else if (m_textToChange)
            {
                m_textToChange.SetText(newValStr);
            }
        }
    }
}