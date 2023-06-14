using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using ThermoVR.UI.Interfaces;

namespace ThermoVR.UI
{
    public enum UIID : byte
    {
        Readout,
        Sandbox,
        Quiz,
        Graph
    }


    /// <summary>
    /// Coordinates opening/closing between various UI systems
    /// </summary>
    public class UIHub : MonoBehaviour
    {
        [SerializeField] private UIModule[] m_registered;

        [SerializeField] private UIID m_initialUI;

        public void Init() {
            for (int i = 0; i < m_registered.Length; i++) {
                m_registered[i].Init();
            }

            OpenUI(m_initialUI);
        }

        public void OpenUI(UIID id) {
            UIModule toOpen = null;

            for (int i = 0; i < m_registered.Length; i++) {
                UIModule currModule = m_registered[i];
                if (currModule.ID == id) {
                    toOpen = currModule;
                }
                else {
                    currModule.Close();
                }
            }

            toOpen.Open();
        }
    }
}

