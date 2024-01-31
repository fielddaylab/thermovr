using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Lab
{
    public class LabMgr : MonoBehaviour
    {
        public static LabMgr Instance;

        public LabStatsMgr Stats;

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

            LabStats newStats = new LabStats();
            newStats.CompletionState = new bool[labInfo.Topics.Count][];
            for (int i = 0; i < labInfo.Topics.Count; i++)
            {
                newStats.CompletionState[i] = new bool[labInfo.Topics[i].Tasks.Count];
            }
            newStats.Progress = 0;

            Stats.LabMap.Add(labInfo.ID, newStats);
        }

        #endregion // Handlers
    }
}
