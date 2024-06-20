using System.Collections;
using System.Collections.Generic;
using System.Text;
using TMPro;
using UnityEngine;
using UnityEngine.UI;

namespace ThermoVR
{
    public class UIInputProxy : MonoBehaviour
    {
        [SerializeField] private CanvasGroup m_mainGroup;
        [SerializeField] private TMP_InputField m_inputField;
        [SerializeField] private TMP_Text m_inputUnitsText;
        [SerializeField] private Button m_closeButton;
        [SerializeField] private Button m_enterButton;

        private InputProxy m_activeProxy;

        private void Start()
        {
            m_mainGroup.alpha = 0;
            m_activeProxy = null;
            GameMgr.Events.Register<InputProxy>(GameEvents.InputProxySelected, HandleInputProxySelected);

            m_closeButton.onClick.AddListener(HandleCloseClicked);
            m_enterButton.onClick.AddListener(HandleEnterClicked);

            m_inputField.onValueChanged.AddListener(HandleValueChanged);
        }

        private void HandleInputProxySelected(InputProxy proxy)
        {
            m_mainGroup.alpha = 1;
            m_activeProxy = proxy;

            m_inputField.text = string.Empty;
            m_inputField.Select();
            m_inputUnitsText.text = proxy.InputUnits;
        }

        private void HandleCloseClicked()
        {
            m_mainGroup.alpha = 0;
            m_activeProxy = null;
        }

        private void HandleEnterClicked()
        {
            m_activeProxy.SetValue(m_inputField.text);
        }

        private void HandleValueChanged(string newVal)
        {

        }
    }
}