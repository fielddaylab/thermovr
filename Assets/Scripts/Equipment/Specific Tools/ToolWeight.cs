using BeauRoutine;
using System.Collections;
using System.Collections.Generic;
using ThermoVR.Physics;
using UnityEngine;

namespace ThermoVR.Tools
{
    public class ToolWeight : Tool
    {
        private static float ENTRY_TIME = 0.2f;

        #region Inspector

        [SerializeField] private FixedAnchor m_HydraulicPressAnchor;
        [SerializeField] private MeshRenderer m_ArrowRenderer;
        [SerializeField] private Material m_LitMat;
        [SerializeField] private Material m_UnlitMat;

        #endregion // Inspector

        #region Tool

        protected override void InitializeRoutines_Impl() {
            m_DeactivatedBasePos = new Vector3(0, 0.6f, 0);
        }

        protected override IEnumerator ActivationRoutine() {
            // TODO: on first activation, object occasionally dips too low

            transform.localPosition = m_DeactivatedBasePos;
            m_ActivatedBasePos = transform.InverseTransformPoint(m_HydraulicPressAnchor.GetAnchorPoint());

            // disable anchor until activation completed
            m_HydraulicPressAnchor.enabled = false;

            gameObject.SetActive(true);

            yield return transform.MoveTo(m_ActivatedBasePos, ENTRY_TIME / m_RoutineSpeed, Axis.Y, Space.Self);

            m_HydraulicPressAnchor.enabled = true;
            yield return null;
        }

        protected override IEnumerator DeactivationRoutine() {
            // disable anchor until deactivation completed
            m_HydraulicPressAnchor.enabled = false;

            yield return transform.MoveTo(m_DeactivatedBasePos, ENTRY_TIME / m_RoutineSpeed, Axis.Y, Space.Self);

            gameObject.SetActive(false);

            m_HydraulicPressAnchor.enabled = true;
            yield return null;
        }

        protected override IEnumerator EngageRoutine() {
            Material[] materials = m_ArrowRenderer.materials;
            materials[0] = m_LitMat;
            m_ArrowRenderer.materials = materials;

            yield return null;
        }

        protected override IEnumerator DisengageRoutine() {
            Material[] materials = m_ArrowRenderer.materials;
            materials[0] = m_UnlitMat;
            m_ArrowRenderer.materials = materials;

            yield return null;
        }

        protected override IEnumerator BeginAdjustRoutine() {
            yield return null;
        }

        protected override IEnumerator EndAdjustRoutine() {
            yield return null;
        }

        #endregion // Tool
    }
}