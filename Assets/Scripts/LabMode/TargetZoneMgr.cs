using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Tools
{
    public class TargetZoneMgr : MonoBehaviour
    {
        private Vector3 DEFAULT_ZONE_DIMS = new Vector3(0.1f, 0.1f, 0.1f);

        [SerializeField] private GameObject m_targetZone;

        private bool m_isGameModeZone = false;

        private void Start()
        {
            EventMgr.Events?.Register<Tuple<Vector3, Vector3>>(GameEvents.TargetZoneUpdated, HandleTargetZoneUpdated);
            EventMgr.Events?.Register(GameEvents.ClearTargetZone, HandleClearTargetZone);
            EventMgr.Events?.Register(GameEvents.GameModeStarted, HandleGameModeStarted);
            EventMgr.Events?.Register(GameEvents.GameModeExited, HandleGameModeExited);
        }

        #region Handlers

        private void HandleTargetZoneUpdated(Tuple<Vector3, Vector3> args)
        {
            m_targetZone.transform.localPosition = args.Item1;

            m_targetZone.transform.localScale = m_isGameModeZone ? DEFAULT_ZONE_DIMS : args.Item2;

            m_targetZone.SetActive(true);
        }

        private void HandleClearTargetZone()
        {
            m_targetZone.SetActive(false);
        }

        private void HandleGameModeStarted()
        {
            m_isGameModeZone = true;
        }

        private void HandleGameModeExited()
        {
            m_isGameModeZone = false;
        }

        #endregion // Handlers
    }
}
