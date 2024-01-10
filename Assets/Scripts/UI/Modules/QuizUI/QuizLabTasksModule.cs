using System;
using System.Collections.Generic;
using ThermoVR.UI;
using UnityEngine;

namespace ThermoVR.Lab
{
    // Manages loading of task tabs
    public class QuizLabTasksModule : UIModule
    {
        [SerializeField] private RectTransform m_tabContainer;
        [SerializeField] private GameObject m_tabPrefab;
        [SerializeField] private GameObject m_taskFramePrefabMC;
        [SerializeField] private GameObject m_taskFramePrefabMCMulti;
        [SerializeField] private GameObject m_taskFramePrefabWordBank;
        [SerializeField] private GameObject m_taskFramePrefabReachState;

        [SerializeField] private float m_xSpacing;

        private List<LabTab> m_tabs; // TODO: make these pools
        private List<LabTaskFrame> m_frames;
        private List<TaskInfo> m_tasks;

        private int m_activeTabIndex;

        #region IUIModule

        public override void Init() {
            base.Init();

            GameMgr.Events?.Register<LabInfo>(GameEvents.ActivateLab, HandleActivateLab);
            GameMgr.Events?.Register<LabInfo>(GameEvents.DeactivateLab, HandleDeactivateLab);

            GameMgr.Events?.Register(GameEvents.TaskResetPressed, HandleTaskResetPressed);


            m_tabs = new List<LabTab>();
            m_frames = new List<LabTaskFrame>();
            m_activeTabIndex = -1;
        }

        public override void Open() {
            base.Open();

            if (m_tasks == null) {
                return;
            }

            for (int i = 0; i < m_tasks.Count; i++) {
                GameObject newTabObj = Instantiate(m_tabPrefab, m_tabContainer.transform);
                newTabObj.transform.localPosition += new Vector3(m_xSpacing * i, 0, 0);

                LabTab newTab = newTabObj.GetComponent<LabTab>();
                m_tabs.Add(newTab);

                newTab.Button.SetText("Task " + (i + 1));
                int tabIndex = i;
                newTab.Button.OnButtonPressed += delegate { HandleLabTabPressed(tabIndex); };

                LabTaskFrame newFrame = PopulateLabTaskFrame(m_tasks[tabIndex]);

                float startingOffset = 0.1f;
                RectTransform tabRect = newTabObj.GetComponent<RectTransform>();
                float tabBuffer = m_xSpacing - tabRect.sizeDelta.x * tabRect.localScale.x;
                m_tabContainer.sizeDelta = new Vector2((startingOffset + tabRect.sizeDelta.x * tabRect.localScale.x + tabBuffer) * m_tasks.Count, m_tabContainer.sizeDelta.y);

                if (newFrame != null) {
                    newTab.RegisterFrame(newFrame);
                }
            }

            // Open first tab
            if (m_tabs.Count > 0) {
                HandleLabTabPressed(0);
            }

            // Disable placement ball
        }

        public override void Close() {
            base.Close();

            ResetWorldMods();

            for (int i = 0; i < m_tabs.Count; i++) {
                Destroy(m_tabs[i].gameObject);
            }
            m_tabs.Clear();


            for (int i = 0; i < m_frames.Count; i++) {
                Destroy(m_frames[i].gameObject);
            }
            m_frames.Clear();

            m_activeTabIndex = -1;
        }

        #endregion // IUIModule

        private LabTaskFrame PopulateLabTaskFrame(TaskInfo taskInfo) {
            GameObject newFrameObj = null;
            LabTaskFrame newFrame = null;

            bool framePopulated = true;
            Evaluable[] evaluables = null;

            switch (taskInfo.TaskType) {
                case TaskType.MultipleChoice:
                    // populate multiple choice
                    newFrameObj = Instantiate(m_taskFramePrefabMC, this.transform);
                    newFrame = newFrameObj.GetComponent<LabTaskFrame>();
                    evaluables = newFrame.GetEvaluables();
                    for (int i = 0; i < evaluables.Length; i++) {
                        MultipleChoiceHub mcHub = evaluables[i].GetComponent<MultipleChoiceHub>();

                        MultipleChoiceDefinition newDef = new MultipleChoiceDefinition();
                        newDef.QuestionText = taskInfo.TaskQuestion;
                        newDef.OptionTexts = taskInfo.SecondaryTexts.ToArray();
                        newDef.CorrectIDs = taskInfo.CorrectIDs;

                        mcHub.SetDefinition(newDef);
                    }

                    break;
                case TaskType.MultipleChoiceMulti:
                    // populate multiple choice
                    newFrameObj = Instantiate(m_taskFramePrefabMCMulti, this.transform);
                    newFrame = newFrameObj.GetComponent<LabTaskFrame>();
                    evaluables = newFrame.GetEvaluables();
                    for (int i = 0; i < evaluables.Length; i++) {
                        MultipleChoiceHubMulti mcHub = evaluables[i].GetComponent<MultipleChoiceHubMulti>();

                        MultipleChoiceDefinition newDef = new MultipleChoiceDefinition();
                        newDef.QuestionText = taskInfo.TaskQuestion;
                        newDef.OptionTexts = taskInfo.SecondaryTexts.ToArray();
                        newDef.CorrectIDs = taskInfo.CorrectIDs;

                        mcHub.SetDefinition(newDef);
                    }

                    break;
                case TaskType.WordBank:
                    // populate word bank
                    newFrameObj = Instantiate(m_taskFramePrefabWordBank, this.transform);
                    newFrame = newFrameObj.GetComponent<LabTaskFrame>();
                    evaluables = newFrame.GetEvaluables();
                    for (int i = 0; i < evaluables.Length; i++) {
                        WordBankHub wordBankHub = evaluables[i].GetComponent<WordBankHub>();

                        WordBankDefinition newDef = new WordBankDefinition();
                        newDef.QuestionText = taskInfo.TaskQuestion;
                        newDef.OptionTexts = taskInfo.SecondaryTexts.ToArray();
                        newDef.CorrectID = taskInfo.CorrectIDs[0];

                        wordBankHub.SetDefinition(newDef);
                    }

                    break;
                case TaskType.ReachState:
                    // populate reach state instructions
                    newFrameObj = Instantiate(m_taskFramePrefabReachState, this.transform);
                    newFrame = newFrameObj.GetComponent<LabTaskFrame>();
                    evaluables = newFrame.GetEvaluables();
                    for (int i = 0; i < evaluables.Length; i++) {
                        ReachStateHub reachStateHub = evaluables[i].GetComponent<ReachStateHub>();

                        ReachStateDefinition newDef = new ReachStateDefinition();
                        newDef.QuestionText = taskInfo.TaskQuestion;
                        newDef.Targets = taskInfo.Targets;

                        reachStateHub.SetDefinition(newDef);
                    }


                    break;
                default:
                    framePopulated = false;
                    break;
            }



            if (framePopulated) {
                newFrameObj.SetActive(false);
                m_frames.Add(newFrame);
            }

            return newFrame;
        }

        private void ApplyWorldMods(TaskInfo mods) {
            // Tools
            World.Instance.ModMgr.SetAllowedTools(mods.AllowedTools);

            // Sets
            int validCount = 0;
            if (mods.Sets.P != -1) { validCount++; }
            if (mods.Sets.V != -1) { validCount++; }
            if (mods.Sets.T != -1) { validCount++; }
            if (validCount >= 2) {
                World.Instance.WarpPVTPartial(mods.Sets.P, mods.Sets.V, mods.Sets.T);
            }

            // Limits
            World.Instance.ModMgr.SetLimits(mods.Limits);
        }

        private void ResetWorldMods() {
            // Tools
            World.Instance?.ModMgr.ResetToolRestrictions();

            // Sets
            // N/A

            // Limits
            World.Instance?.ModMgr.ResetLimits();

            // Restore move ball functionality
            World.Instance?.ModMgr.EnableGraphBallInteractions();
        }

        #region Handlers

        private void HandleActivateLab(LabInfo info) {
            m_tasks = info.Tasks;
        }

        private void HandleDeactivateLab(LabInfo info) {
            m_tasks = null;
        }

        private void HandleLabTabPressed(int newTabIndex) {
            if (m_activeTabIndex == -1) {
                // no tab activated yet; activate new tab
                ActivateTab(newTabIndex);
            }
            else if (m_activeTabIndex == newTabIndex) {
                // same as current tab; no change needed
            }
            else {
                // hide existing tab
                DeactivateTab(m_activeTabIndex);

                // activate new tab
                ActivateTab(newTabIndex);
            }

            m_activeTabIndex = newTabIndex;
        }

        private void HandleTaskResetPressed() {
            ApplyWorldMods(m_tasks[m_activeTabIndex]);
        }

        #endregion // Handlers

        private void ActivateTab(int index) {
            m_frames[index].gameObject.SetActive(true);
            m_tabs[index].Button.SetColor(GameDB.Instance.TabSelectedColor);
            ApplyWorldMods(m_tasks[index]);
        }

        private void DeactivateTab(int index) {
            m_frames[index].gameObject.SetActive(false);
            m_tabs[index].Button.SetColor(GameDB.Instance.TabDefaultColor);
        }
    }

}