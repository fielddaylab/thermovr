using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR;
using ThermoVR.Lab;
using ThermoVR.State;
using ThermoVR.UI;
using ThermoVR.UI.Interfaces;
using TMPro;
using UnityEngine;

public class QuizModule : UIModule
{
    #region Inspector

    [SerializeField] private TMP_Text m_headerText;
    [SerializeField] private TMP_Text m_activeText;

    #endregion // Inspector

    private bool m_labIsActive;


    #region IUIModule

    public override void Init() {
        base.Init();

        m_labIsActive = false;

        GameMgr.Events?.Register<LabInfo>(GameEvents.ActivateLab, HandleActivateLab);
        GameMgr.Events?.Register(GameEvents.DeactivateLab, HandleDeactivateLab);

        GameMgr.Events?.Register(GameEvents.BeginLab, HandleBeginLab);

        m_hub.InitializeRegistered();
    }

    public override void Open() {
        base.Open();

        if (m_labIsActive) {
            m_hub.OpenUI(UIID.QuizLabTasks);
        }
        else {
            m_hub.OpenUI(UIID.QuizSelect);
        }
    }

    public override void Close() {
        base.Close();
    }

    #endregion // IUIModule

    #region Handlers

    
    private void HandleActivateLab(LabInfo labInfo) {
        Debug.Log("[Lab Activation] Lab " + labInfo.Name + " activated!");

        m_labIsActive = true;
        m_hub.OpenUI(UIID.QuizLabTasks);
    }

    private void HandleDeactivateLab()
    {
        m_labIsActive = false;
        m_hub.OpenUI(UIID.QuizSelect);
    }

    private void HandleBeginLab() {
        m_hub.OpenUI(UIID.QuizLabTasks);
    }

    #endregion // Handlers
}
