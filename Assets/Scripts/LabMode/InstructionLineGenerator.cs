using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Lab
{
    public class InstructionLineGenerator : MonoBehaviour
    {
        [SerializeField] private GameObject m_linePrefab;
        [SerializeField] private Transform m_lineContainer;

        [SerializeField] private float m_spacing = 1.3f;

        public void GenerateLines(string[] lines)
        {
            float adjustedSpacing = m_spacing;
            float overrideTextHeight = 0;
            float overrideImgHeight = 0;
            if (lines.Length > 2)
            {
                adjustedSpacing -= 0.6f;
                overrideTextHeight -= 4;
                overrideImgHeight -= 6.6f;
            }

            for (int i = 0; i < lines.Length; i++)
            {
                var transform = GenerateLine(lines[i], overrideTextHeight, overrideImgHeight);
                transform.localPosition += new Vector3(0, -adjustedSpacing * i, 0);
            }
        }

        private Transform GenerateLine(string text, float overrideTextHeight, float overrideImgHeight)
        {
            InstructionLine instruction = Instantiate(m_linePrefab, m_lineContainer).GetComponent<InstructionLine>();
            instruction.Text.SetText(text);
            instruction.Rect.sizeDelta = new Vector2(instruction.Rect.sizeDelta.x, instruction.Rect.sizeDelta.y + overrideTextHeight);
            instruction.Img.rectTransform.sizeDelta = new Vector2(instruction.Img.rectTransform.sizeDelta.x, instruction.Img.rectTransform.sizeDelta.y + overrideImgHeight);
            return instruction.transform;
        }
    }
}
