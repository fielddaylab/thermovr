using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR
{
    public enum PlatformType { 
        VR,
        DESKTOP
    }

    [DefaultExecutionOrder(-5)]
    public class ModeMgr : MonoBehaviour
    {
        public static ModeMgr Instance;

        public bool IsAlphaRelease = true; // temp solution to managing alpha release channel
        public bool IsDesktop = false;

        public PlatformType Platform;

        private void Awake()
        {
            if (Instance == null)
            {
                Instance = this;
                Platform = IsDesktop ? PlatformType.DESKTOP : PlatformType.VR;
            }
            else if (this != Instance)
            {
                Destroy(this.gameObject);
                return;
            }
        }
    }
}