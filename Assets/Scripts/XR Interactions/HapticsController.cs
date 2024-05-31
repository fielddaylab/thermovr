using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Controls
{
    public class HapticsController : MonoBehaviour
    {
        private OVRHapticsClip m_hapticsClip;

        private void Start()
        {
            GameMgr.Events.Register<Hand>(GameEvents.HandStartPress, OnHandStartPress, this);

            byte[] samples = new byte[50];
            for (int i = 0; i < samples.Length; i++)
            {
                samples[i] = 50; // out of 255
            }
            m_hapticsClip = new OVRHapticsClip(samples, samples.Length);
        }

        private void OnHandStartPress(Hand inHand)
        {
            if (inHand == Hand.LEFT)
            {
                OVRHaptics.Channels[0].Mix(m_hapticsClip);
            }
            else if (inHand == Hand.RIGHT)
            {
                OVRHaptics.Channels[1].Mix(m_hapticsClip);
            }
        }
    }
}
