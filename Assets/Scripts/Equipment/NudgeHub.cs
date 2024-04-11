using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR.UI;
using UnityEngine;

namespace ThermoVR.Dials
{
    public class NudgeHub : MonoBehaviour
    {
        [SerializeField] private Material[] m_dialIconMats;

        [SerializeField] private MeshRenderer m_dialIconRenderer;

        [SerializeField] private Pressable m_NudgeUpButton;
        [SerializeField] private Pressable m_NudgeDownButton;

        private int m_currDialIndex;

        private void Start()
        {
            m_currDialIndex = 0;

            GameMgr.Events.Register<int>(GameEvents.SelectNewDial, HandleSelectNewDial);

            m_NudgeUpButton.OnPress += HandleNudgeUpPressed;
            m_NudgeDownButton.OnPress += HandleNudgeDownPressed;
        }

        #region Handlers


        private void HandleSelectNewDial(int dialIndex)
        {
            m_currDialIndex = dialIndex;
            var mats = m_dialIconRenderer.sharedMaterials;
            mats[0] = m_dialIconMats[dialIndex];
            m_dialIconRenderer.sharedMaterials = mats;
        }


        private void HandleNudgeUpPressed(object sender, EventArgs args)
        {
            GameMgr.Events.Dispatch(GameEvents.NudgeUpClicked, m_currDialIndex);

        }

        private void HandleNudgeDownPressed(object sender, EventArgs args)
        {
            GameMgr.Events.Dispatch(GameEvents.NudgeDownClicked, m_currDialIndex);
        }

        #endregion // Handlers
    }
}
