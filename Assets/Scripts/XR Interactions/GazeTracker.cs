using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Controls
{
    public class GazeTracker : MonoBehaviour
    {
        #region Inspector

        public Transform GazeOrigin;
        [SerializeField] private LayerMask GazeLayer;
        
        #endregion // Inspector

        private GazeTargetType m_CurrentGazeTarget;
        private float m_GazeTime;

        #region Unity Callbacks

        private void Update()
        {
            UpdateGaze();
        }

        #endregion // Unity Callbacks

        #region Helpers

        private void UpdateGaze()
        {
            // Raycast for target
            UnityEngine.Physics.Raycast(GazeOrigin.position, GazeOrigin.forward, out RaycastHit hit, Mathf.Infinity, GazeLayer);

            if (hit.collider)
            {
                var target = hit.collider.GetComponent<GazeTarget>();

                if (m_CurrentGazeTarget == GazeTargetType.NONE && target.TargetType != GazeTargetType.NONE)
                {
                    // Begin gaze on new target
                    StartGaze(target.TargetType);
                }
                else if (m_CurrentGazeTarget == target.TargetType)
                {
                    // continue gazing at same object
                    ContinueGaze();
                }
                else
                {
                    // stop gazing at previous object
                    EndGaze();

                    // gaze at new object
                    StartGaze(target.TargetType);
                }
            }
            else if (m_CurrentGazeTarget != GazeTargetType.NONE)
            {
                // stop gazing at previous target
                EndGaze();
            }
        }

        private void StartGaze(GazeTargetType targetType)
        {
            Debug.Log("[Gaze Tracker] Beginning gaze on type " + targetType);

            // set new target and timer
            m_GazeTime = 0;
            m_CurrentGazeTarget = targetType;
        }

        private void ContinueGaze()
        {
            Debug.Log("[Gaze Tracker] continuing gaze on type " + m_CurrentGazeTarget);

            m_GazeTime += Time.deltaTime;
        }

        private void EndGaze()
        {
            Debug.Log("[Gaze Tracker] ending gaze on type " + m_CurrentGazeTarget + ". Time was " + m_GazeTime);

            // dispatch event
            GameMgr.Events.Dispatch(GameEvents.GazeEnd, new Tuple<GazeTargetType, float>(m_CurrentGazeTarget, m_GazeTime));

            // reset target
            m_CurrentGazeTarget = GazeTargetType.NONE;
        }

        #endregion // Helpers
    }
}
