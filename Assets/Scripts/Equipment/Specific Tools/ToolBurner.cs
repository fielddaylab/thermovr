using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Tools
{
    public class ToolBurner : Tool
    {
        #region Tool

        protected override IEnumerator ActivationRoutine() {
            gameObject.SetActive(true);
            yield return null;
        }

        protected override IEnumerator DeactivationRoutine() {
            gameObject.SetActive(false);
            yield return null;
        }

        #endregion // Tool
    }
}