using System.Collections;
using System.Collections.Generic;
using ThermoVR.State;
using ThermoVR;
using ThermoVR.UI;
using ThermoVR.UI.Interfaces;
using UnityEngine;
using TMPro;

public class ReadoutModule : UIModule
{
    #region Inspector

    [SerializeField] private TextMeshPro text_pressure;
    [SerializeField] private TextMeshPro text_temperature;
    [SerializeField] private TextMeshPro text_volume;
    [SerializeField] private TextMeshPro text_internalenergy;
    [SerializeField] private TextMeshPro text_entropy;
    [SerializeField] private TextMeshPro text_enthalpy;
    [SerializeField] private TextMeshPro text_quality;
    [SerializeField] private TextMeshPro text_region;

    #endregion // Inspector

    #region IUIModule

    public override void Init() {
        base.Init();

        GameMgr.Events?.Register<VarUpdate>(GameEvents.UpdateVarText, HandleUpdateVarText, this);
    }

    public override void Open() {
        base.Open();
    }

    public override void Close() {
        base.Close();
    }

    #endregion // IUIModule

    #region Handlers

    private void HandleUpdateVarText(VarUpdate update) {
        TextMeshPro toUpdate = null;
        switch (update.ID) {
            case VarID.Pressure:
                toUpdate = text_pressure;
                break;
            case VarID.Temperature:
                toUpdate = text_temperature;
                break;
            case VarID.Volume:
                toUpdate = text_volume;
                break;
            case VarID.InternalEnergy:
                toUpdate = text_internalenergy;
                break;
            case VarID.Entropy:
                toUpdate = text_entropy;
                break;
            case VarID.Enthalpy:
                toUpdate = text_enthalpy;
                break;
            case VarID.Quality:
                toUpdate = text_quality;
                break;
            case VarID.Region:
                toUpdate = text_region;
                break;
            default:
                // No matching id found!
                return;
        }

        if (toUpdate!= null) {
            toUpdate.SetText(update.NewText);
        }
    }

    #endregion // Handlers
}
