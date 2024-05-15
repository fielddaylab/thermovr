using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Linq;
using ThermoVR.Controls;
using ThermoVR.Lab;
using System;

namespace ThermoVR
{
    public class TracerManager : MonoBehaviour
    {
        public static TracerManager Instance;

        #region Structs & Enums

        private enum State
        {
            Stopped,
            Tracing
        }

        #endregion // Structs & Enums

        #region Inspector

        [SerializeField] private GameObject m_graphBall;
        [SerializeField] private TrailRenderer m_trail;
        [SerializeField] private int m_maxPositions;

        [SerializeField] private Material m_visibleMat;
        [SerializeField] private Material m_hiddenMat;

        #endregion // Inspector

        private State m_state;
        private bool m_isEnabled;
        private bool m_isVisible;
        private bool m_clearExisting;

        #region Unity Callbacks

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

            EndTrace();
            HideTrace();

            GameMgr.Events.Register<Hand>(GameEvents.GraphBallGrabbed, HandleGraphBallGrabbed, this);
            GameMgr.Events?.Register<Tuple<double, double, double>>(GameEvents.WarpPVT, HandleWarpPVT);
        }

        private void Update()
        {
            if (m_state == State.Tracing)
            {
                if (m_maxPositions == -1) { return; }

                // enforce max length
                if (m_trail.positionCount > m_maxPositions)
                {
                    var positions = new Vector3[m_trail.positionCount];
                    m_trail.GetPositions(positions);
                    var positionsList = positions.ToList();

                    positionsList.RemoveRange(0, m_trail.positionCount - m_maxPositions);
                    m_trail.Clear();
                    m_trail.AddPositions(positionsList.ToArray());
                }
            }
        }

        #endregion // Unity Callbacks

        public void ShowTrace()
        {
            m_isVisible = true;

            var mats = m_trail.sharedMaterials;
            mats[0] = m_visibleMat;
            m_trail.sharedMaterials = mats;
        }

        public void HideTrace()
        {
            m_isVisible = false;

            var mats = m_trail.sharedMaterials;
            mats[0] = m_hiddenMat;
            m_trail.sharedMaterials = mats;
        }

        public void StartTrace()
        {
            m_state = State.Tracing;
        }

        /*
        public void ResumeTrace()
        {

        }

        public void PauseTrace()
        {

        }
        */

        public void EndTrace()
        {
            m_state = State.Stopped;
            m_trail.Clear();
        }

        public void LoadSettings(TrailGroup settings)
        {
            m_isEnabled = settings.IsEnabled;
            m_clearExisting = settings.ClearExisting;
            m_maxPositions = settings.MaxLength;
        }

        #region Handlers

        private void HandleGraphBallGrabbed(Hand hand)
        {
            EndTrace();
        }

        private void HandleWarpPVT(Tuple<double, double, double> args)
        {
            EndTrace();

            // enable / disable line
            if (m_isEnabled)
            {
                StartTrace();
                ShowTrace();
            }
            else
            {
                HideTrace();
                EndTrace();
            }
        }

        #endregion // Handlers
    }
}