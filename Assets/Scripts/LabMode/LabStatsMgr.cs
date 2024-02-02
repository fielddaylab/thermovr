using BeauUtil;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Lab
{
    public struct LabStats
    {
        public bool[][] CompletionState; // an entry for each task, separated into topics. True if completed, false otherwise.
        public float Progress;

        public void RefreshProgress()
        {
            int tally = 0;
            int total = 0;

            for (int i = 0; i < CompletionState.Length; i++)
            {
                for (int j = 0; j < CompletionState[i].Length; j++)
                {
                    total++;

                    if (CompletionState[i][j])
                    {
                        tally++;
                    }
                }
            }

            if (total == 0)
            {
                Progress = 1;
            }
            else
            {
                Progress = tally / (float)total;
            }
        }
    }

    public class LabStatsMgr : MonoBehaviour
    {
        public Dictionary<StringHash32, LabStats> LabMap;

        private void Awake()
        {
            LabMap = new Dictionary<StringHash32, LabStats>();
        }
    }
}