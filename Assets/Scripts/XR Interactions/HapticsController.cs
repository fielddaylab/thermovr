using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Controls
{
    public class HapticsController : MonoBehaviour
    {
        private OVRHapticsClip m_pressHapticsClip;
        private OVRHapticsClip m_detentHapticsClip;

        private void Start()
        {
            if (GameMgr.I.IsDesktop) { return; }
            
            GameMgr.Events.Register<Hand>(GameEvents.HandStartPress, OnHandStartPress, this);
            GameMgr.Events.Register<Hand>(GameEvents.DetentHit, OnDetentHit, this);

            // Press Haptics
            byte[] samples = new byte[50];
            for (int i = 0; i < samples.Length; i++)
            {
                samples[i] = 50; // out of 255
            }

            if (OVRHaptics.Config.SampleSizeInBytes != 0) {
                m_pressHapticsClip = new OVRHapticsClip(samples, samples.Length);
            }
            else
            {
                m_pressHapticsClip = null;
            }

            // Detent Haptics
            samples = new byte[25];
            for (int i = 0; i < samples.Length; i++)
            {
                samples[i] = 75; // out of 255
            }

            if (OVRHaptics.Config.SampleSizeInBytes != 0)
            {
                m_detentHapticsClip = new OVRHapticsClip(samples, samples.Length);
            }
            else
            {
                m_detentHapticsClip = null;
            }
        }

        private void OnHandStartPress(Hand inHand)
        {
            PlayHaptics(inHand, m_pressHapticsClip);
        }

        private void OnDetentHit(Hand inHand)
        {
            PlayHaptics(inHand, m_detentHapticsClip);
        }

        private void PlayHaptics(Hand inHand, OVRHapticsClip clip)
        {
            if (inHand == Hand.LEFT)
            {
                OVRHaptics.Channels[0].Mix(clip);
            }
            else if (inHand == Hand.RIGHT)
            {
                OVRHaptics.Channels[1].Mix(clip);
            }
        }
    }
}
