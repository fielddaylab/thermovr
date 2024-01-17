using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR.UI;
using TMPro;
using UnityEngine;
using UnityEngine.UI;

namespace ThermoVR.Lab
{
    public class LabOption : MonoBehaviour
    {
        public TMP_Text LabTitle;
        public TMP_Text LabAuthor;
        public ThermoButton LoadButton;

        public BoxCollider Collider;

        [SerializeField] private Slider m_Slider;

        private LabInfo m_Lab;

        public void SetLabInfo(LabInfo labInfo)
        {
            m_Lab = labInfo;
            LoadButton.OnButtonPressed += OnLoadButtonPressed;
        }

        public void FixCollider(BoxCollider compare)
        {
            Collider.center = compare.center;
            Collider.size = compare.size;
        }

        private void OnLoadButtonPressed(object sender, EventArgs args)
        {
            GameMgr.Events?.Dispatch(GameEvents.PreActivateLab, m_Lab);
            GameMgr.Events?.Dispatch(GameEvents.ActivateLab, m_Lab);
        }
    }
}