using System.Collections;
using System.Collections.Generic;
using ThermoVR.Dials;
using UnityEngine;

public interface ITool
{

    public abstract void InitializeRoutines();

    public abstract void TriggerActivation();

    public abstract void TriggerDeactivation();

    public abstract void TriggerEngage();

    public abstract void TriggerDisengage();

    public abstract void TriggerBeginAdjust();

    public abstract void TriggerEndAdjust();
}
