using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR
{
    public class PersistentState : MonoBehaviour
    {
        public static PersistentState Instance;

        public Dictionary<string, bool> Bools = new Dictionary<string, bool>();

        private void Awake()
        {
            if (Instance == null )
            {
                Instance = this;
                DontDestroyOnLoad(this.gameObject);
            }
            else if (this != Instance)
            {
                Destroy(this.gameObject);
                return;
            }
        }
    }
}