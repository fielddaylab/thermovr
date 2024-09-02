using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

namespace ThermoVR
{
    /// <summary>
    /// Unique identifiers for tutorial types
    /// </summary>
    public enum TutorialID : byte
    {
        Nudge
    }

    [DefaultExecutionOrder(10)]
    public class TutorialMgr : MonoBehaviour
    {
        public static TutorialMgr Instance;

        [SerializeField] private GameObject m_nudgeTutorial;
        [SerializeField] private Image m_nudgeImage;
        [SerializeField] private int m_numTimesDisplayNudge = 4;

        private int m_timesSeenNudge = 0;
        private bool m_showingNudge = false;

        public bool ForceNudge = false;

        private void Awake()
        {
            if (Instance == null)
            {
                Instance = this;
            }
            else if (this != Instance)
            {
                Destroy(this.gameObject);
                return;
            }
        }

        private void OnEnable()
        {
            EventMgr.Events.Register(GameEvents.ShowNudgeTutorial, HandleShowNudgeTutorial);
            EventMgr.Events.Register(GameEvents.HideNudgeTutorial, HandleHideNudgeTutorial);

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

        private void HandleShowNudgeTutorial()
        {
            if (ForceNudge)
            {
                m_nudgeTutorial.SetActive(true);

                m_showingNudge = true;

                EventMgr.Events.Dispatch(GameEvents.NudgeHintDisplayed);
            }
            else
            {
                if (m_timesSeenNudge >= m_numTimesDisplayNudge)
                {
                    return;
                }

                m_showingNudge = true;
                m_nudgeTutorial.SetActive(true);

                EventMgr.Events.Dispatch(GameEvents.NudgeHintDisplayed);
            }
        }

        private void HandleHideNudgeTutorial()
        {
            if (ForceNudge)
            {
                m_nudgeTutorial.SetActive(false);
            }
            else
            {
                if (m_showingNudge)
                {
                    m_timesSeenNudge++;
                }

                m_nudgeTutorial.SetActive(false);
            }

            m_showingNudge = false;

            EventMgr.Events.Dispatch(GameEvents.NudgeHintHidden);
        }

        #endregion // Handlers
    }
}