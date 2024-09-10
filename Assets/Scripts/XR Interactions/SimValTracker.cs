using System;
using System.Collections;
using System.Collections.Generic;
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

        public void Init()
        {
            P = 0;
            V = 0;
            T = 0;
            u = 0;
            s = 0;
            h = 0;
            x = 0;
        }
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
            InitBuffer(ref m_stateBuffer);

            m_frameCounter = 0;
        }

        private void Update()
        {
            UpdateBuffers(m_frameCounter);

            if (m_frameCounter == SAMPLE_SIZE - 1)
            {
                // dispatch frames
                EventMgr.Events.Dispatch(GameEvents.SimStateData, m_stateBuffer);

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

        private void InitBuffer(ref SimStateDataFrame[] buffer)
        {
            for (int i = 0; i < buffer.Length; i++)
            {
                buffer[i].Init();
            }
        }

        private void UpdateBuffers(int frameCount)
        {
            // add current frame data
            LoadFrameToBuffer(m_thermoPresent, ref m_stateBuffer, frameCount);
        }

        private void LoadFrameToBuffer(ThermoPresent toLoad, ref SimStateDataFrame[] buffer, int frameIndex)
        {
            SimStateDataFrame newDataFrame = buffer[frameIndex];

            newDataFrame.P = toLoad.get_state_var(VarID.Pressure);
            newDataFrame.V = toLoad.get_state_var(VarID.Volume);
            newDataFrame.T = toLoad.get_state_var(VarID.Temperature);
            newDataFrame.u = toLoad.get_state_var(VarID.InternalEnergy);
            newDataFrame.s = toLoad.get_state_var(VarID.Entropy);
            newDataFrame.h = toLoad.get_state_var(VarID.Enthalpy);
            newDataFrame.x = toLoad.get_state_var(VarID.Quality);

            buffer[frameIndex] = newDataFrame;
        }

        #endregion // Helpers
    }
}
