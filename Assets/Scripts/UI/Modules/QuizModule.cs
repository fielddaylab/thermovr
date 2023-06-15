using System.Collections;
using System.Collections.Generic;
using ThermoVR.UI;
using ThermoVR.UI.Interfaces;
using TMPro;
using UnityEngine;

public class QuizModule : UIModule
{
    [SerializeField] private TMP_InputField m_inputField;

    #region IUIModule

    public override void Open() {
        base.Open();

        m_inputField.Select();
        m_inputField.ActivateInputField();
    }

    public override void Close() {
        base.Close();
    }

    #endregion // IUIModule
}
