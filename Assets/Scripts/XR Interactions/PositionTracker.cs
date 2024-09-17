using System;
using System.Collections;
using System.Collections.Generic;
using BeauUtil;
using BeauUtil.Extensions;
using UnityEngine;

namespace ThermoVR.Controls
{
    [Serializable]
    public struct PositionDataFrame
    {
        public unsafe fixed float pos[3];
        public unsafe fixed float rot[4];

        public unsafe Vector3 posVector {
            get {
                fixed (float* p = pos) {
                    return Unsafe.FastReinterpret<float, Vector3>(p);
                }
            }
            set {
                fixed(float* p = pos) {
                    *(Vector3*) p = value;
                }
            }
        }

        public unsafe Quaternion rotQuat {
            get {
                fixed (float* r = rot) {
                    return Unsafe.FastReinterpret<float, Quaternion>(r);
                }
            }
            set {
                fixed (float* r = rot) {
                    *(Quaternion*) r = value;
                }
            }
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
            //InitBuffer(ref m_viewportBuffer);
            //InitBuffer(ref m_leftHandBuffer);
            //InitBuffer(ref m_rightHandBuffer);

            m_frameCounter = 0;
        }

        private void Update()
        {
            UpdateBuffers(m_frameCounter);

            if (m_frameCounter == SAMPLE_SIZE - 1)
            {
                // dispatch frames
                EventMgr.Events.Dispatch(GameEvents.ViewportData, EvtArgs.Ref(m_viewportBuffer));
                EventMgr.Events.Dispatch(GameEvents.LeftHandData, EvtArgs.Ref(m_rightHandBuffer));
                EventMgr.Events.Dispatch(GameEvents.RightHandData, EvtArgs.Ref(m_leftHandBuffer));

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

        //private void InitBuffer(ref PositionDataFrame[] buffer)
        //{
        //    for (int i = 0; i < buffer.Length; i++)
        //    {
        //        buffer[i].Init();
        //    }
        //}

        private void UpdateBuffers(int frameCount)
        {
            // add current frame data
            LoadFrameToBuffer(m_viewport, ref m_viewportBuffer, frameCount);
            LoadFrameToBuffer(m_leftHand, ref m_leftHandBuffer, frameCount);
            LoadFrameToBuffer(m_rightHand, ref m_rightHandBuffer, frameCount);

            EventMgr.Events.Dispatch(GameEvents.HeadsetPosUpdated, EvtArgs.Create(m_viewportBuffer[frameCount]));
        }

        private unsafe void LoadFrameToBuffer(Transform toLoad, ref PositionDataFrame[] buffer, int frameIndex)
        {
            PositionDataFrame newDataFrame = buffer[frameIndex];

            toLoad.GetPositionAndRotation(out Vector3 pos, out Quaternion rot);

            newDataFrame.pos[0] = pos.x;
            newDataFrame.pos[1] = pos.y;
            newDataFrame.pos[2] = pos.z;
            newDataFrame.rot[0] = rot.x;
            newDataFrame.rot[1] = rot.y;
            newDataFrame.rot[2] = rot.z;
            newDataFrame.rot[3] = rot.w;

            buffer[frameIndex] = newDataFrame;
        }

        #endregion // Helpers
    }
}
