using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

namespace ThermoVR
{
    [DefaultExecutionOrder(10)]
    public class TutorialMgr : MonoBehaviour
    {
        [SerializeField] private GameObject m_nudgeTutorial;
        [SerializeField] private Image m_nudgeImage;
        [SerializeField] private int m_numTimesDisplayNudge = 4;

        private int m_timesSeenNudge = 0;
        private bool m_showingNudge = false;

        private void OnEnable()
        {
            EventMgr.Events.Register<bool>(GameEvents.ShowNudgeTutorial, HandleShowNudgeTutorial);
            EventMgr.Events.Register<bool>(GameEvents.HideNudgeTutorial, HandleHideNudgeTutorial);

            if (ModeMgr.Instance.IsDesktop)
            {
                m_nudgeImage.sprite = GameDB.Instance.TutorialNudgeDesktop;
            }
            else
            {
                m_nudgeImage.sprite = GameDB.Instance.TutorialNudgeVR;
            }

            m_nudgeTutorial.SetActive(false);
        }

        #region Handlers

        private void HandleShowNudgeTutorial(bool overrideShow)
        {
            if (overrideShow)
            {
                m_nudgeTutorial.SetActive(true);
                return;
            }

            if (m_timesSeenNudge >= m_numTimesDisplayNudge)
            {
                return;
            }

            m_showingNudge = true;
            m_nudgeTutorial.SetActive(true);
        }

        private void HandleHideNudgeTutorial(bool overrideHide)
        {
            if (overrideHide)
            {
                m_nudgeTutorial.SetActive(false);
                return;
            }


            if (m_showingNudge)
            {
                m_timesSeenNudge++;
            }
            m_showingNudge = false;

            m_nudgeTutorial.SetActive(false);
        }

        #endregion // Handlers
    }
}