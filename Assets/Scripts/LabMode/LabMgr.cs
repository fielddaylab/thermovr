using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Lab
{
    public class LabMgr : MonoBehaviour
    {
        public static LabMgr Instance;

        public List<LabInfo> AvailableLabs;
        // dict to list of completed questions

        private void Awake() {
            if (Instance == null)
            {
                Instance = this;
            }
            else if (Instance != this)
            {
                Destroy(this.gameObject);
                return;
            }

            AvailableLabs = new List<LabInfo>();
            GameMgr.Events.Register<LabInfo>(GameEvents.LabLoaded, HandleLabLoaded);
        }

        #region Handlers

        private void HandleLabLoaded(LabInfo labInfo) {
            AvailableLabs.Add(labInfo);
        }

        #endregion // Handlers
    }
}
