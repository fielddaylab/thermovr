using System.Text;
using BeauUtil;
using TMPro;
using UnityEngine;

namespace ThermoVR
{
    static public class UIUtility
    {
        static public bool TryPopulate(this TMP_Text text, string textString)
        {
            if (string.IsNullOrEmpty(textString))
            {
                text.gameObject.SetActive(false);
                return false;
            }

            text.SetText(textString);
            text.gameObject.SetActive(true);
            return true;
        }

        static public bool TryPopulate(this TMP_Text text, StringSlice textString)
        {
            if (textString.IsEmpty)
            {
                text.gameObject.SetActive(false);
                return false;
            }

            text.SetText(textString);
            text.gameObject.SetActive(true);
            return true;
        }

        static public bool TryPopulate(this TMP_Text text, StringBuilder textString)
        {
            if (textString.Length == 0)
            {
                text.gameObject.SetActive(false);
                return false;
            }

            text.SetText(textString);
            text.gameObject.SetActive(true);
            return true;
        }

        static public void ApplyOffsets(RectTransform rect, RectOffset offset)
        {
            rect.offsetMin = new Vector2(offset.left, offset.bottom);
            rect.offsetMax = new Vector2(-offset.right, -offset.top);
        }

        static public void SetAnchorsPreservePosition(this RectTransform rect, Vector2 min, Vector2 max)
        {
            Vector3 pos = rect.position;
            rect.anchorMin = min;
            rect.anchorMax = max;
            rect.position = pos;
        }
    }
}