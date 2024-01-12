using System;
using System.Collections.Generic;
using ThermoVR.UI;
using TMPro;
using UnityEngine;

namespace ThermoVR.Lab
{
    // Manages loading of task tabs
    public class QuizLabTasksModule : UIModule
    {
        private const float TAB_HEIGHT_INACTIVE = 60f;
        private const float TAB_HEIGHT_ACTIVE = 70.66f;

        [SerializeField] private RectTransform m_tabContainer;
        [SerializeField] private GameObject m_tabPrefab;
        [SerializeField] private RectTransform m_topicTabContainer;
        [SerializeField] private GameObject m_topicTabPrefab;
        [SerializeField] private GameObject m_taskFramePrefabMC;
        [SerializeField] private GameObject m_taskFramePrefabMCMulti;
        [SerializeField] private GameObject m_taskFramePrefabWordBank;
        [SerializeField] private GameObject m_taskFramePrefabReachState;

        [SerializeField] private float m_xSpacing;
        [SerializeField] private float m_ySpacing;

        [SerializeField] private ThermoButton m_homeButton;

        [SerializeField] private TMP_Text m_topicHeader;

        private List<LabTopicTab> m_tabs; // TODO: make these pools
        private List<List<LabTaskFrame>> m_frames;
        // private List<TaskInfo> m_tasks;

        private LabInfo m_currLab;
        private bool m_labIsActive;

        private int m_activeTabIndex;
        private int m_activeTopicIndex;

        #region IUIModule

        public override void Init() {
            base.Init();

            GameMgr.Events?.Register<LabInfo>(GameEvents.ActivateLab, HandleActivateLab);
            GameMgr.Events?.Register(GameEvents.DeactivateLab, HandleDeactivateLab);

            GameMgr.Events?.Register(GameEvents.TaskResetPressed, HandleTaskResetPressed);

            m_homeButton.OnButtonPressed += HandleHomeButtonPressed;

            m_tabs = new List<LabTopicTab>();
            m_frames = new List<List<LabTaskFrame>>();
            m_activeTabIndex = -1;
            m_activeTopicIndex = -1;
        }

        public override void Open() {
            base.Open();

            if (!m_labIsActive) {
                return;
            }
            m_tabs.Clear();
            m_frames.Clear();

            for (int topicIndex = 0; topicIndex < m_currLab.Topics.Count; topicIndex++) {
                m_frames.Add(new List<LabTaskFrame>());

                // create a LabTopicButton, HandleLabTopicTabPressed
                GameObject newTopicTabObj = Instantiate(m_topicTabPrefab, m_topicTabContainer.transform);
                newTopicTabObj.transform.localPosition -= new Vector3(0, m_ySpacing * topicIndex, 0);

                LabTopicTab newTopicTab = newTopicTabObj.GetComponent<LabTopicTab>();
                m_tabs.Add(newTopicTab);

                newTopicTab.Button.SetText("" + (topicIndex + 1));
                int topicTabIndex = topicIndex;
                newTopicTab.Button.OnButtonPressed += delegate { HandleLabTopicTabPressed(topicTabIndex); };

                float topicStartingOffset = 0.1f;
                RectTransform topicTabRect = newTopicTab.GetComponent<RectTransform>();
                float topicTabBuffer = m_ySpacing - topicTabRect.sizeDelta.y * topicTabRect.localScale.y;
                // m_topicTabContainer.sizeDelta = new Vector2((topicStartingOffset + topicTabRect.sizeDelta.x * topicTabRect.localScale.x + topicTabBuffer) * m_currLab.Topics.Count, m_topicTabContainer.sizeDelta.y);

                // Create lab task tabs
                for (int taskIndex = 0; taskIndex < m_currLab.Topics[topicIndex].Tasks.Count; taskIndex++)
                {
                    float startingOffset = 0.1f;

                    GameObject newTabObj = Instantiate(m_tabPrefab, m_tabContainer.transform);
                    newTabObj.transform.localPosition += new Vector3(startingOffset + m_xSpacing * taskIndex, 0, 0);

                    LabTab newTab = newTabObj.GetComponent<LabTab>();
                    m_tabs[topicIndex].TaskTabs.Add(newTab);

                    newTab.Button.SetText("Task " + (taskIndex + 1));
                    int currTabIndex = taskIndex;
                    int currTopicIndex = topicIndex;
                    newTab.Button.OnButtonPressed += delegate { HandleLabTabPressed(currTopicIndex, currTabIndex); };

                    LabTaskFrame newFrame = PopulateLabTaskFrame(m_currLab.Topics[topicIndex].Tasks[currTabIndex], currTopicIndex);

                    RectTransform tabRect = newTabObj.GetComponent<RectTransform>();
                    float tabBuffer = m_xSpacing - tabRect.sizeDelta.x * tabRect.localScale.x;
                    // m_tabContainer.sizeDelta = new Vector2((startingOffset + tabRect.sizeDelta.x * tabRect.localScale.x + tabBuffer) * m_currLab.Topics[topicIndex].Tasks.Count, m_tabContainer.sizeDelta.y);

                    if (newFrame != null)
                    {
                        newTab.RegisterFrame(newFrame);
                    }
                }
            }

            // Open first tab
            if (m_tabs.Count > 0) {
                HandleLabTopicTabPressed(0);
            }

            // Disable placement ball
        }

        public override void Close() {
            base.Close();

            ResetWorldMods();

            for (int i = 0; i < m_tabs.Count; i++) {
                for (int j = 0; j < m_tabs[i].TaskTabs.Count; j++) {
                    Destroy(m_tabs[i].TaskTabs[j].gameObject);
                }
            }
            m_tabs.Clear();


            for (int i = 0; i < m_frames.Count; i++) {
                for (int j = 0; j < m_frames[i].Count; j++) {
                    Destroy(m_frames[i][j].gameObject);
                }
            }
            m_frames.Clear();

            m_activeTabIndex = -1;
            m_activeTopicIndex = -1;
        }

        #endregion // IUIModule

        private LabTaskFrame PopulateLabTaskFrame(TaskInfo taskInfo, int currTopicIndex) {
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
                        newDef.InitialConditionText = taskInfo.InitialConditions;
                        newDef.QuestionTexts = taskInfo.TaskQuestions.ToArray();
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
                        newDef.InitialConditionText = taskInfo.InitialConditions;
                        newDef.QuestionTexts = taskInfo.TaskQuestions.ToArray();
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
                        newDef.InitialConditionText = taskInfo.InitialConditions;
                        newDef.QuestionTexts = taskInfo.TaskQuestions.ToArray();
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
                        newDef.InitialConditionText = taskInfo.InitialConditions;
                        newDef.QuestionTexts = taskInfo.TaskQuestions.ToArray();
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
                m_frames[currTopicIndex].Add(newFrame);
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
            m_currLab = info;
            m_labIsActive = true;
        }

        private void HandleDeactivateLab() {
            m_labIsActive = false;
        }

        private void HandleLabTabPressed(int newTopicIndex, int newTabIndex)
        {
            if (m_activeTabIndex == -1)
            {
                // no tab activated yet; activate new tab
                ActivateTab(newTopicIndex, newTabIndex);
            }
            else if (m_activeTopicIndex == newTopicIndex && m_activeTabIndex == newTabIndex)
            {
                // same as current tab; no change needed
            }
            else
            {
                // hide existing tab
                DeactivateTab(m_activeTopicIndex, m_activeTabIndex);

                // activate new tab
                ActivateTab(newTopicIndex, newTabIndex);
            }

            m_activeTabIndex = newTabIndex;
            m_activeTopicIndex = newTopicIndex;
        }

        private void HandleLabTopicTabPressed(int newTopicIndex) {
            if (m_activeTopicIndex == -1) {
                // no tab activated yet; activate new tab
                ActivateTopicTab(newTopicIndex);
            }
            else if (m_activeTopicIndex == newTopicIndex) {
                // same as current tab; no change needed
            }
            else {
                // hide existing tab
                DeactivateTopicTab(m_activeTopicIndex);

                // activate new tab
                ActivateTopicTab(newTopicIndex);
            }

            m_activeTabIndex = 0;
            m_activeTopicIndex = newTopicIndex;
        }

        private void HandleTaskResetPressed() {
            ApplyWorldMods(m_currLab.Topics[m_activeTopicIndex].Tasks[m_activeTabIndex]);
        }


        private void HandleHomeButtonPressed(object sender, EventArgs args)
        {
            GameMgr.Events?.Dispatch(GameEvents.DeactivateLab);
        }

        #endregion // Handlers

        private void ActivateTab(int topicIndex, int taskIndex)
        {
            // TODO: highlight active tab

            m_frames[topicIndex][taskIndex].gameObject.SetActive(true);
            m_tabs[topicIndex].TaskTabs[taskIndex].Button.SetColor(GameDB.Instance.TabSelectedColor);
            ApplyWorldMods(m_currLab.Topics[topicIndex].Tasks[taskIndex]);
        }

        private void DeactivateTab(int topicIndex, int taskIndex)
        {
            // TODO: un-highlight prev active tab

            m_frames[topicIndex][taskIndex].gameObject.SetActive(false);
            m_tabs[topicIndex].TaskTabs[taskIndex].Button.SetColor(GameDB.Instance.TabDefaultColor);
        }

        private void ActivateTopicTab(int topicIndex)
        {
            // TODO: highlight active topic tab
            m_topicHeader.SetText(m_currLab.Topics[topicIndex].TopicHeader);

            for (int i = 0; i < m_tabs[topicIndex].TaskTabs.Count; i++)
            {
                m_tabs[topicIndex].TaskTabs[i].gameObject.SetActive(true);
            }

            m_frames[topicIndex][0].gameObject.SetActive(true);
            m_tabs[topicIndex].TaskTabs[0].Button.SetColor(GameDB.Instance.TabSelectedColor);
            ApplyWorldMods(m_currLab.Topics[topicIndex].Tasks[0]);
        }

        private void DeactivateTopicTab(int topicIndex)
        {
            // TODO: un-highlight prev active topic tab
            m_topicHeader.SetText("");

            for (int i = 0; i < m_tabs[topicIndex].TaskTabs.Count; i++)
            {
                m_frames[topicIndex][i].gameObject.SetActive(false);
                m_tabs[topicIndex].TaskTabs[i].gameObject.SetActive(false);
                m_tabs[topicIndex].TaskTabs[i].Button.SetColor(GameDB.Instance.TabDefaultColor);
            }
        }
    }

}