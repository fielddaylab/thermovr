using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Tools
{
    public class ToolChamberTemp : Tool
    {
        #region Tool

        protected override void InitializeRoutines_Impl() {

        }
        protected override IEnumerator ActivationRoutine() {
            gameObject.SetActive(true);
            yield return null;
        }

        protected override IEnumerator DeactivationRoutine() {
            gameObject.SetActive(false);
            yield return null;
        }

        protected override IEnumerator BeginAdjustRoutine() {
            yield return null;
        }

        protected override IEnumerator EndAdjustRoutine() {
            yield return null;
        }

        protected override IEnumerator EngageRoutine() {
            yield return null;
        }

        protected override IEnumerator DisengageRoutine() {
            yield return null;
        }

        #endregion // Tool
    }
}