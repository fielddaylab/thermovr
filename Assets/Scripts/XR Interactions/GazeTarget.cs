using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Controls
{
    public enum GazeTargetType
    {
        NONE,
        TABLET,
        PISTON,
        GRAPH,
        CONTROLS,
    }

    public class GazeTarget : MonoBehaviour
    {
        public GazeTargetType TargetType;
    }
}
