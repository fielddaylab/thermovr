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

        EventMgr.Events?.Register<LabInfo>(GameEvents.ActivateLab, HandleActivateLab);
        EventMgr.Events?.Register(GameEvents.DeactivateLab, HandleDeactivateLab);

        EventMgr.Events?.Register(GameEvents.BeginLab, HandleBeginLab);

        m_hub.InitializeRegistered();
    }

    public override void Open() {
        this.gameObject.SetActive(true);

        if (m_labIsActive) {
            m_hub.OpenUI(UIID.QuizLabTasks);
        }
        else {
            m_hub.OpenUI(UIID.QuizSelect);
        }
    }

    public override void Close() {
        this.gameObject.SetActive(false);
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
