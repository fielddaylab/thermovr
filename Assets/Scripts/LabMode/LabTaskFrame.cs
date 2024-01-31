using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR.UI;
using TMPro;
using UnityEngine;

namespace ThermoVR.Lab
{
    public class LabTaskFrame : MonoBehaviour
    {
        public AnswerEvaluator AnswerEvaluator;
        public ThermoButton TaskResetButton;

        [SerializeField] AudioClip m_taskResetClip;

        [SerializeField] private Evaluable[] m_evaluables;

        private void OnEnable() {
            TaskResetButton.OnButtonPressed += HandleResetPressed;

            UpdateResetButtonState();
        }

        private void OnDisable() {
            TaskResetButton.OnButtonPressed -= HandleResetPressed;
        }

        public void LoadCompleted(bool completed)
        {
            AnswerEvaluator.LoadCompleted(completed);
        }

        private void Update()
        {
            UpdateResetButtonState();
        }

        private void UpdateResetButtonState()
        {
            bool anyEvaluated = false;
            foreach (var evaluable in m_evaluables)
            {
                if (evaluable.HasBeenEvaluated())
                {
                    anyEvaluated = true;
                }
            }

            TaskResetButton.gameObject.SetActive(anyEvaluated);
        }

        public Evaluable[] GetEvaluables() {
            return m_evaluables;
        }

        private void HandleResetPressed(object sender, EventArgs args) {
            for (int i = 0; i < m_evaluables.Length; i++) {
                m_evaluables[i].ResetState();
            }

            if (GameMgr.I.AudioEnabled) { Tablet.Instance.PlayUIAudio(m_taskResetClip); }

            GameMgr.Events.Dispatch(GameEvents.TaskResetPressed);
        }
    }
}

