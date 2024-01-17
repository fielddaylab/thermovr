using System;
using ThermoVR;
using ThermoVR.UI;
using UnityEngine;
using UnityEngine.UI;
using BeauPools;
using BeauUtil;
using BeauRoutine;
using System.Collections.Generic;

namespace ThermoVR.Lab
{
    public class QuizLabSelectModule : UIModule
    {
        [Serializable] public class LabOptionPool : SerializablePool<LabOption> { }

        #region Inspector

        [SerializeField] private ThermoButton m_ScrollUpBtn;
        [SerializeField] private ThermoButton m_ScrollDownBtn;
        [SerializeField] private Transform m_ScrollContainer;
        [SerializeField] private float m_OptionSpacing = 1.9f;
        public LabOptionPool m_OptionPool;

        #endregion // Inspector

        private List<LabOption> m_ActiveOptions = new List<LabOption>();

        #region IUIModule

        public override void Open()
        {
            this.gameObject.SetActive(true);

            FreeLabs();

            // Add button listeners
            m_ScrollUpBtn.OnButtonPressed -= HandleScrollUp;
            m_ScrollDownBtn.OnButtonPressed -= HandleScrollDown;
            m_ScrollUpBtn.OnButtonPressed += HandleScrollUp;
            m_ScrollDownBtn.OnButtonPressed += HandleScrollDown;

            // Spawn and populate available labs
            for (int i = 0; i < LabMgr.Instance.AvailableLabs.Count; i++)
            {
                var option = m_OptionPool.Alloc(m_ScrollContainer);
                option.transform.localPosition = new Vector3(option.transform.localPosition.x, option.transform.localPosition.y - (m_OptionSpacing * i), option.transform.localPosition.z);
                option.LabTitle.text = LabMgr.Instance.AvailableLabs[i].Name;
                option.LabAuthor.text = LabMgr.Instance.AvailableLabs[i].Author;
                option.LoadButton.ClearListeners();
                option.SetLabInfo(LabMgr.Instance.AvailableLabs[i]);
                option.FixCollider(m_OptionPool.Prefab.Collider);

                // TODO: set Slider progess
                m_ActiveOptions.Add(option);
            }
        }

        public override void Close()
        {
            this.gameObject.SetActive(false);

            // Remove button listeners
            m_ScrollUpBtn.OnButtonPressed -= HandleScrollUp;
            m_ScrollDownBtn.OnButtonPressed -= HandleScrollDown;

            // Free labs
            FreeLabs();
        }

        #endregion // IUIModule

        #region Handlers

        private void HandleScrollUp(object sender, EventArgs args)
        {
            if (m_ScrollContainer.transform.localPosition.y <= 0) { return; }
            m_ScrollContainer.transform.localPosition = m_ScrollContainer.transform.localPosition - new Vector3(0, m_OptionSpacing, 0);
        }

        private void HandleScrollDown(object sender, EventArgs args)
        {
            if (m_ScrollContainer.transform.localPosition.y >= ((m_ActiveOptions.Count - 3) * m_OptionSpacing)) { return; }
            m_ScrollContainer.transform.localPosition = m_ScrollContainer.transform.localPosition + new Vector3(0, m_OptionSpacing, 0);
        }

        #endregion // Handlers

        private void FreeLabs()
        {
            foreach (LabOption option in m_ActiveOptions)
            {
                m_OptionPool.Free(option);
            }
            m_ActiveOptions.Clear();
        }
    }

}