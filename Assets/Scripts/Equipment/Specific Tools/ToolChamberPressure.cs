using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Tools
{
    public class ToolChamberPressure : Tool
    {
        private static int GLASS_INDEX = 1;
        [SerializeField] private MeshRenderer m_Mesh;
        [SerializeField] private Material m_ActiveGlassMat;
        [SerializeField] private Material m_InactiveGlassMat;

        #region Tool

        protected override void InitializeRoutines_Impl() {

        }

        protected override IEnumerator ActivationRoutine() {
            // ALWAYS ENGAGED
            /*
            Material[] materials = m_Mesh.materials;
            materials[GLASS_INDEX] = m_ActiveGlassMat;
            m_Mesh.materials = materials;
            */

            // gameObject.SetActive(true);
            yield return null;
        }

        protected override IEnumerator DeactivationRoutine() {
            // ALWAYS ENGAGED
            /*
            Material[] materials = m_Mesh.materials;
            materials[GLASS_INDEX] = m_InactiveGlassMat;
            m_Mesh.materials = materials;
            */

            // gameObject.SetActive(false);
            yield return null;
        }

        protected override IEnumerator BeginAdjustRoutine() {
            yield return null;
        }

        protected override IEnumerator EndAdjustRoutine() {
            yield return null;
        }

        protected override IEnumerator EngageRoutine() {
            yield return null;
        }

        protected override IEnumerator DisengageRoutine() {
            yield return null;
        }

        #endregion // Tool
    }
}