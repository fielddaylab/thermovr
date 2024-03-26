using BeauRoutine;
using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR;
using ThermoVR.Tools;
using UnityEngine;

public class RotatorButton : MonoBehaviour
{
    [SerializeField] private CoordinatedState m_rotationState;

    [SerializeField] private Pressable m_button;
    [SerializeField] private float m_rotationStep;

    private void OnEnable() {
        m_button.OnPress += HandleButtonPressed;

        m_rotationState.SetInitialVal(m_rotationState.GetSharedTarget().transform.localEulerAngles.y);
    }

    private void OnDisable() {
        m_button.OnPress -= HandleButtonPressed;
    }

    #region Handlers

    private void HandleButtonPressed(object sender, EventArgs args) {
        if (!m_rotationState.GetSharedTarget()) {
            return;
        }

        // Try to apply rotation
        if (m_rotationState.TryModifyVal(m_rotationStep)) {
            // within bounds; apply rotation
            m_rotationState.ReplaceCoordinatedRoutine(RotationRoutine());

            if (GameMgr.I.AudioEnabled) { m_button.ClickAudio.Play(); }
        }
    }

    #endregion // Handlers

    #region Routines

    private IEnumerator RotationRoutine() {
        // new target rotation is initial rotation + new rotation amount
        float targetYRotation = m_rotationState.GetInitialVal() + m_rotationState.GetVal();

        Transform targetTransform = m_rotationState.GetSharedTarget().transform;
        Vector3 targetRotation = targetTransform.localEulerAngles;
        targetRotation.y = targetYRotation;

        if (m_rotationStep > 0) { GameMgr.Events.Dispatch(GameEvents.RotateGraphClickedCW, new Tuple<float, float>(targetTransform.eulerAngles.y, targetRotation.y)); }
        else { GameMgr.Events.Dispatch(GameEvents.RotateGraphClickedCCW, new Tuple<float, float>(targetTransform.eulerAngles.y, targetRotation.y)); }
        
        yield return targetTransform.RotateTo(targetRotation, 0.5f, Axis.Y);

        yield return null;
    }

    #endregion // Routines
}
