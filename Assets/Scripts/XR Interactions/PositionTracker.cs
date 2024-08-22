using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Controls
{
    [Serializable]
    public struct PositionDataFrame
    {
        public float[] pos;
        public float[] rot;

        public void Init()
        {
            pos = new float[3];
            rot = new float[4];
        }
    }

    public class PositionTracker : MonoBehaviour
    {
        private const int SAMPLE_SIZE = 30;

        [SerializeField] private Transform m_viewport;
        [SerializeField] private Transform m_leftHand;
        [SerializeField] private Transform m_rightHand;

        private PositionDataFrame[] m_viewportBuffer;
        private PositionDataFrame[] m_leftHandBuffer;
        private PositionDataFrame[] m_rightHandBuffer;

        private int m_frameCounter;

        #region Unity Callbacks

        private void OnEnable()
        {
            m_viewportBuffer = new PositionDataFrame[SAMPLE_SIZE];
            m_leftHandBuffer = new PositionDataFrame[SAMPLE_SIZE];
            m_rightHandBuffer = new PositionDataFrame[SAMPLE_SIZE];
            InitBuffer(ref m_viewportBuffer);
            InitBuffer(ref m_leftHandBuffer);
            InitBuffer(ref m_rightHandBuffer);

            m_frameCounter = 0;
        }

        private void Update()
        {
            UpdateBuffers(m_frameCounter);

            if (m_frameCounter == SAMPLE_SIZE - 1)
            {
                // dispatch frames
                EventMgr.Events.Dispatch(GameEvents.ViewportData, m_viewportBuffer);
                EventMgr.Events.Dispatch(GameEvents.LeftHandData, m_rightHandBuffer);
                EventMgr.Events.Dispatch(GameEvents.RightHandData, m_leftHandBuffer);

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

        private void InitBuffer(ref PositionDataFrame[] buffer)
        {
            for (int i = 0; i < buffer.Length; i++)
            {
                buffer[i].Init();
            }
        }

        private void UpdateBuffers(int frameCount)
        {
            // add current frame data
            LoadFrameToBuffer(m_viewport, ref m_viewportBuffer, frameCount);
            LoadFrameToBuffer(m_leftHand, ref m_leftHandBuffer, frameCount);
            LoadFrameToBuffer(m_rightHand, ref m_rightHandBuffer, frameCount);

            EventMgr.Events.Dispatch(GameEvents.HeadsetPosUpdated, m_viewportBuffer[frameCount]);
        }

        private void LoadFrameToBuffer(Transform toLoad, ref PositionDataFrame[] buffer, int frameIndex)
        {
            PositionDataFrame newDataFrame = buffer[frameIndex];

            newDataFrame.pos[0] = toLoad.position.x;
            newDataFrame.pos[1] = toLoad.position.y;
            newDataFrame.pos[2] = toLoad.position.z;
            newDataFrame.rot[0] = toLoad.rotation.x;
            newDataFrame.rot[1] = toLoad.rotation.y;
            newDataFrame.rot[2] = toLoad.rotation.z;
            newDataFrame.rot[3] = toLoad.rotation.w;

            buffer[frameIndex] = newDataFrame;
        }

        #endregion // Helpers
    }
}
