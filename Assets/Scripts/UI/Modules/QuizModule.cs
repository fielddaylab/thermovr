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
    [SerializeField] private ThermoButton m_beginButton;

    private Cartridge m_activeCartridge;

    #endregion // Inspector

    #region IUIModule

    public override void Init() {
        base.Init();

        m_activeCartridge = null;

        GameMgr.Events?.Register<Cartridge>(GameEvents.ActivateCartridge, HandleActivateCartridge);
        GameMgr.Events?.Register<Cartridge>(GameEvents.DeactivateCartridge, HandleDeactivateCartridge);
    }

    public override void Open() {
        base.Open();

        if (m_activeCartridge) {
            DisplayLoadedScreen();
        }
        else {
            HideLoadedScreen();
        }
    }

    public override void Close() {
        base.Close();

        m_beginButton.OnButtonPressed -= HandleBeginPressed;
    }

    #endregion // IUIModule

    #region Handlers

    private void HandleActivateCartridge(Cartridge cartridge) {
        switch (cartridge.GetCartridgeType()) {
            case Cartridge.CartridgeType.Lab:
                Debug.Log("[Cartridge] Lab " + cartridge.GetInfo().Name + " activated!");
                m_activeText.SetText("LOADED: \n" + cartridge.GetInfo().Name);
                m_activeCartridge = cartridge;
                DisplayLoadedScreen();
                break;
            case Cartridge.CartridgeType.Sandbox:
                Debug.Log("[Cartridge] Sandbox activated!");
                m_activeText.SetText("LOADED: \n" + "Sandbox");
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

        m_activeText.SetText("");
    }

    private void HandleBeginPressed(object sender, EventArgs args) {
        Debug.Log("[Press] Begin Pressed!");
        HideLoadedScreen();
    }

    #endregion // Handlers

    private void DisplayLoadedScreen() {
        m_beginButton.OnButtonPressed += HandleBeginPressed;

        m_beginButton.gameObject.SetActive(true);
        m_activeText.gameObject.SetActive(true);
    }

    private void HideLoadedScreen() {
        m_beginButton.OnButtonPressed -= HandleBeginPressed;

        m_beginButton.gameObject.SetActive(false);
        m_activeText.gameObject.SetActive(false);
    }
}
