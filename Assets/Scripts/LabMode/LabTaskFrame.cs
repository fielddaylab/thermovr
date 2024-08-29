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
        public ThermoButton NextButton;

        [SerializeField] AudioClip m_taskResetClip;

        [SerializeField] private Evaluable[] m_evaluables;
        [SerializeField] private bool m_usesSubmit = true;

        private void OnEnable() {
            if (TaskResetButton)
            {
                TaskResetButton.OnButtonPressed += HandleResetPressed;
            }
            if (NextButton)
            {
                NextButton.OnButtonPressed += HandleNextPressed;
                NextButton.SetInteractable(false);
            }

            bool anyEvaluated = AnyEvaluated();

            UpdateResetButtonState(anyEvaluated);
        }

        private void OnDisable() {
            if (TaskResetButton)
            {
                TaskResetButton.OnButtonPressed -= HandleResetPressed;
            }
            if (NextButton)
            {
                NextButton.OnButtonPressed -= HandleNextPressed;
            }
        }

        public void LoadCompleted(bool completed)
        {
            AnswerEvaluator.LoadCompleted(completed);
        }

        private bool AnyEvaluated()
        {
            bool anyEvaluated = false;
            foreach (var evaluable in m_evaluables)
            {
                if (evaluable.HasBeenEvaluated())
                {
                    anyEvaluated = true;
                }
            }

            return anyEvaluated;
        }

        private void Update()
        {
            bool anyEvaluated = AnyEvaluated();

            UpdateResetButtonState(anyEvaluated);
            UpdateNextButtonState(anyEvaluated);
        }

        private void UpdateResetButtonState(bool anyEvaluated)
        {
            if (TaskResetButton)
            {
                TaskResetButton.gameObject.SetActive(anyEvaluated);
            }
        }

        private void UpdateNextButtonState(bool anyEvaluated)
        {
            if (NextButton)
            {
                NextButton.SetInteractable(anyEvaluated);
                NextButton.gameObject.SetActive(anyEvaluated);
            }
            
            if (m_usesSubmit)
            {
                AnswerEvaluator.SubmitButton.gameObject.SetActive(!anyEvaluated);
            }
        }

        public Evaluable[] GetEvaluables() {
            return m_evaluables;
        }

        public void ConfigureAsLastTask()
        {
            if (NextButton)
            {
                NextButton.gameObject.SetActive(false);
            }
        }

        private void HandleResetPressed(object sender, EventArgs args) {
            for (int i = 0; i < m_evaluables.Length; i++) {
                m_evaluables[i].ResetState();
            }

            if (GameMgr.I.AudioEnabled) { Tablet.Instance.PlayUIAudio(m_taskResetClip); }

            EventMgr.Events.Dispatch(GameEvents.TaskResetPressed);
        }

        private void HandleNextPressed(object sender, EventArgs args)
        {
            EventMgr.Events.Dispatch(GameEvents.TaskNextPressed);
        }
    }
}

