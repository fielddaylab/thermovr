using BeauRoutine;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

namespace ThermoVR.Controls
{
    public class DesktopTabletController : MonoBehaviour
    {
        private enum State
        {
            Visible,
            Hidden
        }

        [SerializeField] private GameObject m_tablet;
        [SerializeField] private float m_transitionTime;
        [SerializeField] private Button m_toggleButton;
        [SerializeField] private Image m_toggleArrow;

        [SerializeField] private Vector3 m_visiblePos;
        [SerializeField] private Quaternion m_visibleRot;
        [SerializeField] private Vector3 m_hiddenPos;
        [SerializeField] private Quaternion m_hiddenRot;


        private State m_state;
        private Routine m_transitionRoutine;
        private bool m_inTransition;

        private void OnEnable()
        {
            m_inTransition = false;
            m_state = State.Visible;

            m_toggleButton.onClick.AddListener(HandleToggleClicked);
        }

        public void ToggleTabletVisibility()
        {
            if (m_inTransition) { return; }

            switch(m_state)
            {
                case State.Visible:
                    HideTablet();
                    break;
                case State.Hidden:
                    ShowTablet();
                    break;
                default:
                    break;
            }
        }

        public void ShowTablet()
        {
            m_transitionRoutine.Replace(ShowRoutine());
        }

        public void HideTablet()
        {
            m_transitionRoutine.Replace(HideRoutine());
        }

        #region Routines

        private IEnumerator ShowRoutine()
        {
            m_inTransition = true;

            yield return Routine.Combine(
                m_tablet.transform.MoveTo(m_visiblePos, m_transitionTime),
                m_tablet.transform.RotateQuaternionTo(m_visibleRot, m_transitionTime)
            );

            m_toggleArrow.transform.localScale = new Vector3(
                Mathf.Abs(m_toggleArrow.transform.localScale.x),
                m_toggleArrow.transform.localScale.y,
                m_toggleArrow.transform.localScale.z
                );

            m_state = State.Visible;
            m_inTransition = false;
        }

        private IEnumerator HideRoutine()
        {
            m_inTransition = true;

            yield return Routine.Combine(
                m_tablet.transform.MoveTo(m_hiddenPos, m_transitionTime),
                m_tablet.transform.RotateQuaternionTo(m_hiddenRot, m_transitionTime)
            );

            m_toggleArrow.transform.localScale = new Vector3(
                -Mathf.Abs(m_toggleArrow.transform.localScale.x),
                m_toggleArrow.transform.localScale.y,
                m_toggleArrow.transform.localScale.z
                );

            m_state = State.Hidden;
            m_inTransition = false;
        }

        #endregion // Routines

        #region Handlers

        private void HandleToggleClicked()
        {
            ToggleTabletVisibility();
        }

        #endregion // Handlers

        #region Editor

#if UNITY_EDITOR

        [ContextMenu("Set Tablet Visible Pos")]
        private void SetTabletVisiblePos()
        {
            m_visiblePos = m_tablet.transform.position;
            m_visibleRot = m_tablet.transform.rotation;
        }

        [ContextMenu("Set Tablet Hidden Pos")]
        private void SetTabletHiddenPos()
        {
            m_hiddenPos = m_tablet.transform.position;
            m_hiddenRot = m_tablet.transform.rotation;
        }

        [ContextMenu("Apply Tablet Visible Pos")]
        private void ApplyTabletVisiblePos()
        {
            m_tablet.transform.position = m_visiblePos;
            m_tablet.transform.rotation = m_visibleRot;
        }

        [ContextMenu("Apply Tablet Hidden Pos")]
        private void ApplyTabletHiddenPos()
        {
            m_tablet.transform.position = m_hiddenPos;
            m_tablet.transform.rotation = m_hiddenRot;
        }


#endif // UNITY_EDITOR

        #endregion // Editor
    }
}