using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Lab
{
    public class InstructionLineGenerator : MonoBehaviour
    {
        [SerializeField] private GameObject m_linePrefab;
        [SerializeField] private Transform m_lineContainer;

        public Transform GenerateLine(string text)
        {
            InstructionLine instruction = Instantiate(m_linePrefab, m_lineContainer).GetComponent<InstructionLine>();
            instruction.Text.SetText(text);
            return instruction.transform;
        }
    }
}
