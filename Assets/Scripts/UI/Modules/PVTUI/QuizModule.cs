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

    private LabInfo m_activeLab;


    #region IUIModule

    public override void Init() {
        base.Init();

        // m_activeLab = null;

        GameMgr.Events?.Register<LabInfo>(GameEvents.ActivateLab, HandleActivateLab);

        GameMgr.Events?.Register(GameEvents.BeginLab, HandleBeginLab);

        m_hub.InitializeRegistered();
    }

    public override void Open() {
        base.Open();

        /*
        if (m_activeCartridge) {
            m_hub.OpenUI(UIID.QuizLoaded);
        }
        else {
            m_hub.OpenUI(UIID.QuizDefault);
        }
        */

        m_hub.OpenUI(UIID.QuizSelect);
    }

    public override void Close() {
        base.Close();
    }

    #endregion // IUIModule

    #region Handlers

    
    private void HandleActivateLab(LabInfo labInfo) {
        Debug.Log("[Lab Activation] Lab " + labInfo.Name + " activated!");

        m_hub.OpenUI(UIID.QuizLabTasks);
    }

    private void HandleBeginLab() {
        m_hub.OpenUI(UIID.QuizLabTasks);

        // Pass in cartridge data
    }

    #endregion // Handlers
}
