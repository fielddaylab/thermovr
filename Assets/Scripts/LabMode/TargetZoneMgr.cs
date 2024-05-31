using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Tools
{
    public class TargetZoneMgr : MonoBehaviour
    {
        [SerializeField] private GameObject m_targetZone;

        private void Start()
        {
            GameMgr.Events?.Register<Tuple<Vector3, Vector3>>(GameEvents.TargetZoneUpdated, HandleTargetZoneUpdated);
            GameMgr.Events?.Register(GameEvents.ClearTargetZone, HandleClearTargetZone);
        }

        #region Handlers

        private void HandleTargetZoneUpdated(Tuple<Vector3, Vector3> args)
        {
            m_targetZone.transform.localPosition = args.Item1;
            m_targetZone.transform.localScale = args.Item2;

            m_targetZone.SetActive(true);
        }

        private void HandleClearTargetZone()
        {
            m_targetZone.SetActive(false);
        }

        #endregion // Handlers
    }
}
