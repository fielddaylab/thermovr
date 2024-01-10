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

    private Cartridge m_activeCartridge;


    #region IUIModule

    public override void Init() {
        base.Init();

        m_activeCartridge = null;

        GameMgr.Events?.Register<Cartridge>(GameEvents.ActivateCartridge, HandleActivateCartridge);
        GameMgr.Events?.Register<Cartridge>(GameEvents.DeactivateCartridge, HandleDeactivateCartridge);

        GameMgr.Events?.Register(GameEvents.BeginLab, HandleBeginLab);

        m_hub.InitializeRegistered(); //TODO: move away from Init(). This is a bit wonky
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

    private void HandleActivateCartridge(Cartridge cartridge) {
        switch (cartridge.GetCartridgeType()) {
            case Cartridge.CartridgeType.Lab:
                Debug.Log("[Cartridge] Lab " + cartridge.GetInfo().Name + " activated!");
                m_activeText.SetText("LOADED: \n" + cartridge.GetInfo().Name);
                m_activeCartridge = cartridge;
                // m_hub.OpenUI(UIID.QuizLoaded);
                break;
            case Cartridge.CartridgeType.Sandbox:
                /*
                Debug.Log("[Cartridge] Sandbox activated!");
                m_activeText.SetText("LOADED: \n" + "Sandbox");
                m_activeCartridge = cartridge;
                m_hub.OpenUI(UIID.QuizLoaded);
                */
                break;
            default:
                break;
        }
    }

    private void HandleDeactivateCartridge(Cartridge cartridge) {
        Debug.Log("[Cartridge] Cartridge Deactivated.");

        if (m_activeCartridge == cartridge) {
            m_activeCartridge = null;
        }

        // m_hub.OpenUI(UIID.QuizDefault);
    }

    private void HandleBeginLab() {
        m_hub.OpenUI(UIID.QuizLabTasks);

        // Pass in cartridge data
    }

    #endregion // Handlers
}
