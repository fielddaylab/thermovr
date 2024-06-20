using System.Collections;
using System.Collections.Generic;
using ThermoVR.Dials;
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

        public void SetValue(string newValStr)
        {
            float newVal;

            try
            {
                newVal = float.Parse(newValStr);
            }
            catch { return; }

            newVal *= m_mult;

            if (m_dialToSet)
            {
                if (m_dialToSet.val_within_range(newVal))
                {
                    m_dialToSet.convert_and_set_map(newVal);
                }
            }
            else if (m_textToChange)
            {
                m_textToChange.SetText(newValStr);
            }
        }
    }
}