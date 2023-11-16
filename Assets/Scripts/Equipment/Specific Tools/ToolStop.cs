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

        private float m_LeftStartAngle, m_RightStartAngle;
        private float m_StartY;

        private static float ADJUST_ANGLE = 90f;
        private static float ROTATE_SPEED = 0.075f;

        #region Tool

        protected override void InitializeRoutines_Impl() {
            m_LeftStartAngle = m_Left.transform.eulerAngles.y;
            m_RightStartAngle = m_Right.transform.eulerAngles.y;

            m_StartY = m_Left.transform.position.y;

            m_Elements.Add(m_Left);
            m_Elements.Add(m_Right);
        }

        protected override IEnumerator ActivationRoutine() {
            Vector3 currLeft = m_Left.transform.position;
            m_Left.transform.position = new Vector3(currLeft.x, m_StartY, currLeft.z);
            Vector3 currRight = m_Right.transform.position;
            m_Right.transform.position = new Vector3(currRight.x, m_StartY, currRight.z);

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
            yield return Routine.Combine(
                m_Left.transform.RotateTo(m_LeftStartAngle, ROTATE_SPEED / m_RoutineSpeed, Axis.Y, Space.World),
                m_Right.transform.RotateTo(m_RightStartAngle, ROTATE_SPEED / m_RoutineSpeed, Axis.Y, Space.World)
                );
        }

        private IEnumerator SnapOff() {
            yield return Routine.Combine(
                m_Left.transform.RotateTo(m_LeftStartAngle + ADJUST_ANGLE, ROTATE_SPEED / m_RoutineSpeed, Axis.Y, Space.World),
                m_Right.transform.RotateTo(m_RightStartAngle - ADJUST_ANGLE, ROTATE_SPEED / m_RoutineSpeed, Axis.Y, Space.World)
                );
        }


        #endregion // Helpers
    }
}