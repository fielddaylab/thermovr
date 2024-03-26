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
        [SerializeField] private AudioClip m_audioClip;

        private LabInfo m_Lab;
        private int m_LabIndex;

        public void SetLabInfo(LabInfo labInfo, int labIndex)
        {
            m_Lab = labInfo;
            m_LabIndex = labIndex;
            LoadButton.OnButtonPressed += OnLoadButtonPressed;
        }

        public LabInfo GetLabInfo()
        {
            return m_Lab;
        }

        public void SetSlider(float val)
        {
            m_Slider.value = val;
        }

        public void FixCollider(BoxCollider compare)
        {
            Collider.center = compare.center;
            Collider.size = compare.size;
        }

        public void EnableCollider()
        {
            Collider.enabled = true;
        }

        public void DisableCollider()
        {
            Collider.enabled = false;
        }

        private void OnLoadButtonPressed(object sender, EventArgs args)
        {
            if (GameMgr.I.AudioEnabled)
            {
                Tablet.Instance.PlayUIAudio(m_audioClip);
            }

            GameMgr.Events?.Dispatch(GameEvents.PreActivateLab, new Tuple<LabInfo, int>(m_Lab, m_LabIndex));
            GameMgr.Events?.Dispatch(GameEvents.SelectLab);
            GameMgr.Events?.Dispatch(GameEvents.ActivateLab, m_Lab);
        }
    }
}