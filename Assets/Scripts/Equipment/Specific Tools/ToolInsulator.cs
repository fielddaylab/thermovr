using BeauRoutine;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Tools
{
    public class ToolInsulator : Tool
    {
        private static float ENTRY_TIME = 0.2f;

        #region Tool

        protected override void InitializeRoutines_Impl() {
            m_DeactivatedBasePos = new Vector3(0, 1, 0);
            m_ActivatedBasePos = new Vector3(0, 0, 0);
        }

        protected override IEnumerator ActivationRoutine() {
            transform.localPosition = m_DeactivatedBasePos;
            gameObject.SetActive(true);

            yield return transform.MoveTo(m_ActivatedBasePos, ENTRY_TIME / m_RoutineSpeed, Axis.XYZ, Space.Self);
        }

        protected override IEnumerator DeactivationRoutine() {
            yield return transform.MoveTo(m_DeactivatedBasePos, ENTRY_TIME / m_RoutineSpeed, Axis.XYZ, Space.Self);

            gameObject.SetActive(false);
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
