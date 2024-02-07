using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Analytics
{
    [System.Serializable]
    public class TargetTaskLogData : TaskLogData
    {
        public bool IsComplete;
        public Dictionary<string, float> TargetStateTarget;
        public List<string> ConstantVariableTarget;
    }
}