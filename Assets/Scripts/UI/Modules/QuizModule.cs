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

    #region IUIModule

    public override void Init() {
        base.Init();

        GameMgr.Events?.Register<Cartridge>(GameEvents.ActivateCartridge, HandleActivateCartridge);
        GameMgr.Events?.Register<Cartridge>(GameEvents.DeactivateCartridge, HandleDeactivateCartridge);
    }

    public override void Open() {
        base.Open();
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

        m_activeText.SetText("");
    }

    #endregion // Handlers
}
