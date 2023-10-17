using System.Collections;
using System.Collections.Generic;
using ThermoVR.Dials;
using UnityEngine;

public interface ITool
{
    public abstract void TriggerActivation();

    public abstract void TriggerDeactivation();
}
