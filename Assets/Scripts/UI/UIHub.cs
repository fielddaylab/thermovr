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
        Quiz,
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

        public void Awake() {
            for (int i = 0; i < m_registered.Length; i++) {
                m_registered[i].Init();
            }

            if (m_initialUI) {
                OpenUI(m_initialUI.ID);
            }
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

