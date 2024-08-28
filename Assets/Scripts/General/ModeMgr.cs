using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR
{   [DefaultExecutionOrder(-5)]
    public class ModeMgr : MonoBehaviour
    {
        public static ModeMgr Instance;

        public bool IsAlphaRelease = true; // temp solution to managing alpha release channel
        public bool IsDesktop = false;

        private void Awake()
        {
            if (Instance == null)
            {
                Instance = this;
            }
            else if (this != Instance)
            {
                Destroy(this.gameObject);
                return;
            }
        }
    }
}