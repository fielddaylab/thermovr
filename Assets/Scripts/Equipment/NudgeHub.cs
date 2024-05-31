using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR.Tools;
using ThermoVR.UI;
using TMPro;
using UnityEngine;

namespace ThermoVR.Dials
{
    public class NudgeHub : MonoBehaviour
    {
        [SerializeField] private TextMeshPro m_currDialValText;
        [SerializeField] private UnitText m_currDialUnits;
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

            GameMgr.Events.Register(GameEvents.DialTextUpdated, HandleDialTextUpdated);
        }

        #region Handlers


        private void HandleSelectNewDial(int dialIndex)
        {
            m_currDialIndex = dialIndex;
            var mats = m_dialIconRenderer.sharedMaterials;
            mats[0] = m_dialIconMats[dialIndex];
            m_dialIconRenderer.sharedMaterials = mats;

            var tools = ToolMgr.Instance.Dials[m_currDialIndex].GetRelevantTools();
            m_currDialValText.SetText(ToolMgr.Instance.Dials[m_currDialIndex].textv.text);
            if (tools.Count > 0)
            {
                m_currDialUnits.UpdateToolType(tools[0].tool_type);
            }
        }

        private void HandleDialTextUpdated()
        {
            m_currDialValText.SetText(ToolMgr.Instance.Dials[m_currDialIndex].textv.text);
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
