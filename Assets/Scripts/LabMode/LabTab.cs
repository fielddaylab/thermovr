using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR.UI;
using UnityEngine;
using UnityEngine.UI;

namespace ThermoVR.Lab
{
    [RequireComponent(typeof(ThermoButton))]
    public class LabTab : MonoBehaviour
    {
        [SerializeField] private ThermoButton m_thermoButton;
        [SerializeField] private Image m_completionSocket;

        public Image ButtonImage;
        public RectTransform ButtonRect;

        private LabTaskFrame m_taskFrame;

        public ThermoButton Button {
            get { return m_thermoButton; }
        }

        public void OnDisable() {
            if (m_taskFrame != null) {
                m_taskFrame.AnswerEvaluator.OnEvaluationUpdated -= HandleEvalUpdate;
                m_taskFrame.TaskResetButton.OnButtonPressed -= HandleFrameReset;
            }
        }

        public void OnEnable()
        {
            if (m_taskFrame != null)
            {
                m_taskFrame.AnswerEvaluator.OnEvaluationUpdated -= HandleEvalUpdate;
                m_taskFrame.TaskResetButton.OnButtonPressed -= HandleFrameReset;

                m_taskFrame.AnswerEvaluator.OnEvaluationUpdated += HandleEvalUpdate;
                m_taskFrame.TaskResetButton.OnButtonPressed += HandleFrameReset;
            }
        }

        public void RegisterFrame(LabTaskFrame frame) {
            m_taskFrame = frame;

            m_taskFrame.AnswerEvaluator.OnEvaluationUpdated += HandleEvalUpdate;
            m_taskFrame.TaskResetButton.OnButtonPressed += HandleFrameReset;
        }

        public void ShowCompletionSprite()
        {
            m_completionSocket.enabled = true;
        }

        public void HideCompletionSprite()
        {
            m_completionSocket.enabled = false;
        }

        public void SetSocketSprite(Sprite newSprite) {
            m_completionSocket.sprite = newSprite;
        }


        private void HandleEvalUpdate(object sender, BoolEventArgs args) {
            if (args.Value) {
                ShowCompletionSprite();
            }
            else {
                HideCompletionSprite();
            }
        }

        private void HandleFrameReset(object sender, EventArgs args) {
            HideCompletionSprite();
        }
    }
}
