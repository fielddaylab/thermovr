using System.Collections;
using System.Collections.Generic;
using ThermoVR.State;
using ThermoVR;
using ThermoVR.UI;
using ThermoVR.UI.Interfaces;
using UnityEngine;
using TMPro;
using ThermoVR.Lab;

public class ReadoutModule : UIModule
{
    #region Inspector

    [SerializeField] private SensorReadout readout_pressure;
    [SerializeField] private SensorReadout readout_temperature;
    [SerializeField] private SensorReadout readout_volume;
    [SerializeField] private SensorReadout readout_internalenergy;
    [SerializeField] private SensorReadout readout_entropy;
    [SerializeField] private SensorReadout readout_enthalpy;
    [SerializeField] private SensorReadout readout_quality;
    [SerializeField] private SensorReadout readout_region;

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
        SensorReadout toUpdate = null;
        switch (update.ID) {
            case VarID.Pressure:
                toUpdate = readout_pressure;
                break;
            case VarID.Temperature:
                toUpdate = readout_temperature;
                break;
            case VarID.Volume:
                toUpdate = readout_volume;
                break;
            case VarID.InternalEnergy:
                toUpdate = readout_internalenergy;
                break;
            case VarID.Entropy:
                toUpdate = readout_entropy;
                break;
            case VarID.Enthalpy:
                toUpdate = readout_enthalpy;
                break;
            case VarID.Quality:
                toUpdate = readout_quality;
                break;
            case VarID.Region:
                toUpdate = readout_region;
                break;
            default:
                // No matching id found!
                return;
        }

        if (toUpdate!= null) {
            toUpdate.TMP.SetText(update.NewText);

            // TextMeshPro thing;
            // thing.SetText(update.NewText);
        }
    }

    #endregion // Handlers
}
