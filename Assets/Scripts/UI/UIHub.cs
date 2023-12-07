using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using ThermoVR.UI.Interfaces;

namespace ThermoVR.UI
{
    public enum UIID : byte
    {
        // PVT Tablet
        Readout,
        Sandbox,
        Lab,
        Graph,

        // Quiz UI
        QuizDefault,
        QuizLoaded,
        QuizLabTasks
    }


    /// <summary>
    /// Coordinates opening/closing between various UI systems
    /// </summary>
    public class UIHub : MonoBehaviour
    {
        [SerializeField] private UIModule[] m_registered;

        [SerializeField] private UIModule m_initialUI;

        private bool m_initialized = false;

        public void Awake() {
            InitializeRegistered();

            if (m_initialUI) {
                OpenUI(m_initialUI.ID);
            }
        }

        public void InitializeRegistered() {
            if (m_initialized) {
                return;
            }

            for (int i = 0; i < m_registered.Length; i++) {
                m_registered[i].Init();
            }

            m_initialized = true;
        }

        public void OpenUI(UIID id) {
            UIModule toOpen = null;

            for (int i = 0; i < m_registered.Length; i++) {
                UIModule currModule = m_registered[i];
                if (currModule.ID == id) {
                    toOpen = currModule;
                }
                else {
                    if (!currModule.AlwaysAvailable) {
                        currModule.Close();
                    }
                }
            }

            toOpen.Open();
        }
    }
}

