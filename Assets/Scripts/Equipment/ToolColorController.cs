using BeauUtil;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using ThermoVR.Dials;
using UnityEngine;

namespace ThermoVR.Tools
{
    [RequireComponent(typeof(Dial))]
    public class ToolColorController : MonoBehaviour
    {
        #region Inspector

        [SerializeField] private Color m_LowestColor;
        [SerializeField] private Color m_HighestColor;

        [SerializeField] private MeshRenderer[] m_Renderers;

        #endregion // Inspector

        private Dial m_ToolDial;

        #region Unity Callbacks

        private void OnEnable()
        {
            m_ToolDial = this.GetComponent<Dial>();

            if (m_ToolDial)
            {
                m_ToolDial.DialMoved.AddListener(HandleToolValUpdated);
            }
        }

        private void OnDisable()
        {
            if (m_ToolDial)
            {
                m_ToolDial.DialMoved.RemoveListener(HandleToolValUpdated);
            }
        }

        #endregion // Unity Callbacks

        #region Handlers

        private void HandleToolValUpdated()
        {
            if (m_Renderers.Length == 0 || !m_ToolDial)
            {
                return;
            }

            float newVal = m_ToolDial.get_val();

            for (int i = 0; i < m_Renderers.Length; i++)
            {
                var mats = m_Renderers[i].sharedMaterials;
                for (int m = 0; m < mats.Length; m++)
                {
                    mats[m].color = Color.Lerp(m_LowestColor, m_HighestColor, newVal);
                }

                m_Renderers[i].sharedMaterials = mats;
            }
        }

        #endregion // Handlers

    }
}
