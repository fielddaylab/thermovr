using BeauUtil;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using ThermoVR.Dials;
using UnityEngine;

namespace ThermoVR.Tools
{
    [RequireComponent(typeof(Dial))]
    public class ToolTransparencyController : MonoBehaviour
    {
        #region Inspector

        [SerializeField] private float m_lowestAlpha; // alpha at lowest val of slider
        [SerializeField] private float m_highestAlpha; // alpha at highest val of slider

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
                    var newColor = mats[m].color;
                    newColor.a = m_lowestAlpha + newVal * (m_highestAlpha - m_lowestAlpha);
                    mats[m].color = newColor;
                }

                m_Renderers[i].sharedMaterials = mats;
            }
        }

        #endregion // Handlers

    }
}
