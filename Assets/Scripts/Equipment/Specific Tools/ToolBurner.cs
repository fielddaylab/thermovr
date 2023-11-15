using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Tools
{
    public class ToolBurner : Tool
    {
        #region Tool

        protected override IEnumerator ActivationRoutine() {
            Debug.Log("[Triggers] Burner activated!");

            gameObject.SetActive(true);
            yield return null;
        }

        protected override IEnumerator DeactivationRoutine() {
            Debug.Log("[Triggers] Burner deactivated!");

            gameObject.SetActive(false);
            yield return null;
        }

        protected override IEnumerator BeginAdjustRoutine() {
            Debug.Log("[Triggers] Burner begin adjust!");
            yield return null;
        }

        protected override IEnumerator EndAdjustRoutine() {
            Debug.Log("[Triggers] Burner end adjust!");
            yield return null;
        }

        protected override IEnumerator EngageRoutine() {
            Debug.Log("[Triggers] Burner engaged!");
            yield return null;
        }

        protected override IEnumerator DisengageRoutine() {
            Debug.Log("[Triggers] Burner disengaged!");
            yield return null;
        }


        #endregion // Tool
    }
}