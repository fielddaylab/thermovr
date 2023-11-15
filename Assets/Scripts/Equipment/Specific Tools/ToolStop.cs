using BeauRoutine;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Tools
{
    public class ToolStop : Tool
    {
        [SerializeField] private GameObject m_Left;
        [SerializeField] private GameObject m_Right;

        private static float ADJUST_ANGLE = 90f;
        private static float ROTATE_SPEED = 0.075f;

        #region Tool

        protected override void InitializeRoutines_Impl() {

        }

        protected override IEnumerator ActivationRoutine() {
            gameObject.SetActive(true);
            yield return null;
        }

        protected override IEnumerator DeactivationRoutine() {
            gameObject.SetActive(false);
            yield return null;
        }

        protected override IEnumerator BeginAdjustRoutine() {
            yield return SnapOff();
        }

        protected override IEnumerator EndAdjustRoutine() {
            yield return SnapOn();
        }

        protected override IEnumerator EngageRoutine() {
            yield return null;
        }

        protected override IEnumerator DisengageRoutine() {
            yield return null;
        }

        #endregion // Tool

        #region Helpers

        private IEnumerator SnapOn() {
            Vector3 targetLeftRotation = m_Left.transform.localEulerAngles;
            targetLeftRotation.y = targetLeftRotation.y - ADJUST_ANGLE;
            Vector3 targetRightRotation = m_Right.transform.localEulerAngles;
            targetRightRotation.y = targetRightRotation.y + ADJUST_ANGLE;

            yield return Routine.Combine(
                m_Left.transform.RotateTo(targetLeftRotation, ROTATE_SPEED / m_RoutineSpeed, Axis.Y, Space.Self),
                m_Right.transform.RotateTo(targetRightRotation, ROTATE_SPEED / m_RoutineSpeed, Axis.Y, Space.Self)
                );
        }

        private IEnumerator SnapOff() {
            Vector3 targetLeftRotation = m_Left.transform.localEulerAngles;
            targetLeftRotation.y = targetLeftRotation.y + ADJUST_ANGLE;
            Vector3 targetRightRotation = m_Right.transform.localEulerAngles;
            targetRightRotation.y = targetRightRotation.y - ADJUST_ANGLE;

            yield return Routine.Combine(
                m_Left.transform.RotateTo(targetLeftRotation, ROTATE_SPEED / m_RoutineSpeed, Axis.Y, Space.Self),
                m_Right.transform.RotateTo(targetRightRotation, ROTATE_SPEED / m_RoutineSpeed, Axis.Y, Space.Self)
                );
        }


        #endregion // Helpers
    }
}