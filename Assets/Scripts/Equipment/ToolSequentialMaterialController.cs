using BeauUtil;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using ThermoVR.Dials;
using UnityEngine;

namespace ThermoVR.Tools
{
    [RequireComponent(typeof(Dial))]
    public class ToolSequentialMaterialController : MonoBehaviour
    {
        #region Inspector

        [SerializeField] private MeshRenderer[] m_Renderers; // from bottom to top
        [SerializeField] private Material m_LitMat;
        [SerializeField] private Material m_UnlitMat;

        #endregion // Inspector

        private Dial m_ToolDial;
        private float m_LightStep; // value difference for each light to trigger

        #region Unity Callbacks

        private void OnEnable() {
            m_ToolDial = this.GetComponent<Dial>();

            if (m_ToolDial) {
                m_ToolDial.DialMoved.AddListener(HandleToolValUpdated);
            }

            if (m_Renderers.Length > 0) {
                m_LightStep = 1.0f / m_Renderers.Length;
            }
            else {
                m_LightStep = 1.0f;
            }
        }

        private void OnDisable() {
            if (m_ToolDial) {
                m_ToolDial.DialMoved.RemoveListener(HandleToolValUpdated);
            }
        }

        #endregion // Unity Callbacks

        #region Handlers

        private void HandleToolValUpdated() {
            if (m_Renderers.Length == 0 || !m_ToolDial) {
                return;
            }

            float newVal = m_ToolDial.get_val();
            int numLit = (int)((newVal + 0.05f) / m_LightStep);

            for (int i = 0; i < numLit; i++) {
                m_Renderers[i].material = m_LitMat;
            }
            for (int i = numLit; i < m_Renderers.Length; i++) {
                m_Renderers[i].material = m_UnlitMat;
            }
        }

        #endregion // Handlers

    }
}
