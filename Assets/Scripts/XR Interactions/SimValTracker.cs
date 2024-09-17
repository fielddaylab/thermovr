using System;
using System.Collections;
using System.Collections.Generic;
using BeauUtil.Extensions;
using UnityEngine;

namespace ThermoVR.State
{
    [Serializable]
    public struct SimStateDataFrame
    {
        public double P;
        public double V;
        public double T;
        public double u;
        public double s;
        public double h;
        public double x;
    }

    public class SimValTracker : MonoBehaviour
    {
        private const int SAMPLE_SIZE = 30;

        [SerializeField] private ThermoPresent m_thermoPresent;

        private SimStateDataFrame[] m_stateBuffer;

        private int m_frameCounter;

        #region Unity Callbacks

        private void OnEnable()
        {
            m_stateBuffer = new SimStateDataFrame[SAMPLE_SIZE];

            m_frameCounter = 0;
        }

        private void Update()
        {
            UpdateBuffers(m_frameCounter);

            if (m_frameCounter == SAMPLE_SIZE - 1)
            {
                // dispatch frames
                EventMgr.Events.Dispatch(GameEvents.SimStateData, EvtArgs.Ref(m_stateBuffer));

                // reset (old samples will be overriden frame by frame)
                m_frameCounter = 0;
            }
            else
            {
                m_frameCounter++;
            }
        }

        #endregion // Unity Callbacks

        #region Helpers

        private void UpdateBuffers(int frameCount)
        {
            // add current frame data
            LoadFrameToBuffer(m_thermoPresent, ref m_stateBuffer, frameCount);
        }

        private void LoadFrameToBuffer(ThermoPresent toLoad, ref SimStateDataFrame[] buffer, int frameIndex)
        {
            ref SimStateDataFrame newDataFrame = ref buffer[frameIndex];

            newDataFrame.P = toLoad.get_state_var(VarID.Pressure);
            newDataFrame.V = toLoad.get_state_var(VarID.Volume);
            newDataFrame.T = toLoad.get_state_var(VarID.Temperature);
            newDataFrame.u = toLoad.get_state_var(VarID.InternalEnergy);
            newDataFrame.s = toLoad.get_state_var(VarID.Entropy);
            newDataFrame.h = toLoad.get_state_var(VarID.Enthalpy);
            newDataFrame.x = toLoad.get_state_var(VarID.Quality);
        }

        #endregion // Helpers
    }
}
