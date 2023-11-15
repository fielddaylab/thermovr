using BeauRoutine;
using System.Collections;
using System.Collections.Generic;
using ThermoVR.Physics;
using UnityEngine;

namespace ThermoVR.Tools
{
    public class ToolNegativeWeight : Tool
    {
        private static float ENTRY_TIME = 0.2f;


        #region Inspector

        [SerializeField] private FixedAnchor m_ElasticAnchor;

        #endregion // Inspector

        #region Tool

        protected override void InitializeRoutines_Impl() {
            m_DeactivatedBasePos = new Vector3(0, 0.4f, 0);
        }
        protected override IEnumerator ActivationRoutine() {
            // TODO: on first activation, object occasionally dips too low

            transform.localPosition = m_DeactivatedBasePos;
            m_ActivatedBasePos = transform.InverseTransformPoint(m_ElasticAnchor.GetAnchorPoint());

            // disable anchor until activation completed
            m_ElasticAnchor.enabled = false;

            gameObject.SetActive(true);

            yield return transform.MoveTo(m_ActivatedBasePos, ENTRY_TIME / m_RoutineSpeed, Axis.Y, Space.Self);

            m_ElasticAnchor.enabled = true;
            yield return null;
        }

        protected override IEnumerator DeactivationRoutine() {
            // disable anchor for deactivation
            m_ElasticAnchor.enabled = false;

            yield return transform.MoveTo(m_DeactivatedBasePos, ENTRY_TIME / m_RoutineSpeed, Axis.Y, Space.Self);

            gameObject.SetActive(false);
            m_ElasticAnchor.enabled = true;
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