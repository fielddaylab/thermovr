using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Controls
{
    [DefaultExecutionOrder(-5)]
    public class ListenerController : MonoBehaviour
    {
        [SerializeField] private ResonanceAudioListener m_listener;

        private void Awake()
        {
            m_listener.enabled = false;
            GameMgr.Events?.Register(GameEvents.InitialLoadComplete, HandleInitialLoadComplete);
        }

        #region Handlers

        private void HandleInitialLoadComplete()
        {
            m_listener.enabled = true;
        }

        #endregion // Handlers
    }
}