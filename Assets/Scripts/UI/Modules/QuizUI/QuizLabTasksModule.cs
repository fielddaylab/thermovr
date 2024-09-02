using BeauRoutine;
using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR.Analytics;
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

        #region Inspector

        [SerializeField] private RectTransform m_tabContainer;
        [SerializeField] private GameObject m_tabPrefab;
        [SerializeField] private RectTransform m_topicTabContainer;
        [SerializeField] private GameObject m_topicTabPrefab;
        [SerializeField] private GameObject m_taskFramePrefabMC;
        [SerializeField] private GameObject m_taskFramePrefabText;
        [SerializeField] private GameObject m_taskFramePrefabMCMulti;
        [SerializeField] private GameObject m_taskFramePrefabWordBank;
        [SerializeField] private GameObject m_taskFramePrefabReachState;

        [SerializeField] private float m_xSpacing;
        [SerializeField] private float m_ySpacing;

        [SerializeField] private ThermoButton m_homeButton;
        [SerializeField] private AudioClip m_homeAudioClip;

        [SerializeField] private TMP_Text m_topicHeader;

        [Space(5)]
        [Header("Scroll")]
        [SerializeField] private Transform m_ScrollVerticalContainer;
        [SerializeField] private float m_VerticalScrollSpacing;
        [SerializeField] private ThermoButton m_ScrollUpBtn;
        [SerializeField] private ThermoButton m_ScrollDownBtn;
        private int m_ScrollVerticalValidVisibleIndex;
        private const int SCROLL_VERTICAL_NUM = 5;

        [SerializeField] private Transform m_ScrollHorizontalContainer;
        [SerializeField] private float m_HorizontalScrollSpacing;
        [SerializeField] private ThermoButton m_ScrollLeftBtn;
        [SerializeField] private ThermoButton m_ScrollRightBtn;
        private int m_ScrollHorizontalValidVisibleIndex;
        private const int SCROLL_HORIZONTAL_NUM = 4;


        #endregion // Inspector

        private float m_VerticalScrollOrigin;
        private float m_HorizontalScrollOrigin;

        private List<LabTopicTab> m_tabs; // TODO: make these pools
        private List<List<LabTaskFrame>> m_frames;
        // private List<TaskInfo> m_tasks;

        private LabInfo m_currLab;
        private bool m_labIsActive;

        private int m_activeTabIndex;
        private int m_activeTopicIndex;

        private List<IndexedTaskInfo> m_visibleTasks;
        private List<IndexedTopicInfo> m_visibleSections;

        private Routine m_wordBankRoutine;

        #region IUIModule

        public override void Init() {
            base.Init();

            EventMgr.Events?.Register<Tuple<LabInfo, int>>(GameEvents.PreActivateLab, HandlePreActivateLab);
            EventMgr.Events?.Register(GameEvents.DeactivateLab, HandleDeactivateLab);

            EventMgr.Events?.Register(GameEvents.TaskResetPressed, HandleTaskResetPressed);
            EventMgr.Events?.Register(GameEvents.TaskNextPressed, HandleTaskNextPressed);

            EventMgr.Events?.Register(GameEvents.ClickOpenWordBank, HandleWordBankOpened);
            EventMgr.Events?.Register<string>(GameEvents.WordBankClosed, HandleWordBankClosed);


            m_homeButton.OnButtonPressed += HandleHomeButtonPressed;

            m_tabs = new List<LabTopicTab>();
            m_frames = new List<List<LabTaskFrame>>();
            m_activeTabIndex = -1;
            m_activeTopicIndex = -1;
            m_visibleTasks = new List<IndexedTaskInfo>();
            m_visibleSections = new List<IndexedTopicInfo>();

            m_VerticalScrollOrigin = m_ScrollVerticalContainer.localPosition.y;
            m_HorizontalScrollOrigin = m_ScrollHorizontalContainer.localPosition.x;
        }

        private void OnDisable()
        {
            Close();
        }

        public override void Open() {
            // clear previous tabs
            Close();

            this.gameObject.SetActive(true);

            m_tabs.Clear();
            m_frames.Clear();

            // Add button listeners
            m_ScrollUpBtn.OnButtonPressed -= HandleScrollUp;
            m_ScrollDownBtn.OnButtonPressed -= HandleScrollDown;
            m_ScrollLeftBtn.OnButtonPressed -= HandleScrollLeft;
            m_ScrollRightBtn.OnButtonPressed -= HandleScrollRight;

            m_ScrollUpBtn.OnButtonPressed += HandleScrollUp;
            m_ScrollDownBtn.OnButtonPressed += HandleScrollDown;
            m_ScrollLeftBtn.OnButtonPressed += HandleScrollLeft;
            m_ScrollRightBtn.OnButtonPressed += HandleScrollRight;

            for (int topicIndex = 0; topicIndex < m_currLab.Topics.Count; topicIndex++) {
                m_frames.Add(new List<LabTaskFrame>());

                // create a LabTopicButton, HandleLabTopicTabPressed
                GameObject newTopicTabObj = Instantiate(m_topicTabPrefab, m_topicTabContainer.transform);
                newTopicTabObj.transform.localPosition -= new Vector3(0, m_ySpacing * topicIndex, 0);

                LabTopicTab newTopicTab = newTopicTabObj.GetComponent<LabTopicTab>();
                m_tabs.Add(newTopicTab);

                newTopicTab.Button.SetText("" + (topicIndex + 1));
                int topicTabIndex = topicIndex;
                newTopicTab.Button.OnButtonPressed += delegate { HandleLabTopicTabPressed(topicTabIndex, true); };

                float topicStartingOffset = 0.1f;
                RectTransform topicTabRect = newTopicTab.GetComponent<RectTransform>();
                float topicTabBuffer = m_ySpacing - topicTabRect.sizeDelta.y * topicTabRect.localScale.y;
                // m_topicTabContainer.sizeDelta = new Vector2((topicStartingOffset + topicTabRect.sizeDelta.x * topicTabRect.localScale.x + topicTabBuffer) * m_currLab.Topics.Count, m_topicTabContainer.sizeDelta.y);

                bool allTasksInTopicCorrect = true;

                // Create lab task tabs
                for (int taskIndex = 0; taskIndex < m_currLab.Topics[topicIndex].Tasks.Count; taskIndex++)
                {
                    float startingOffset = 0.1f;

                    GameObject newTabObj = Instantiate(m_tabPrefab, m_tabContainer.transform);
                    newTabObj.transform.localPosition += new Vector3(startingOffset + m_xSpacing * taskIndex, 0, 0);

                    int currTopicIndex = topicIndex;

                    LabTab newTab = newTabObj.GetComponent<LabTab>();
                    m_tabs[topicIndex].TaskTabs.Add(newTab);
                    m_tabs[topicIndex].TaskTabs[taskIndex].OnCompletionStateSubmitted += delegate { HandleLabTaskCompletionUpdated(currTopicIndex, false); };
                    m_tabs[topicIndex].TaskTabs[taskIndex].OnCompletionStateReset += delegate { HandleLabTaskCompletionUpdated(currTopicIndex, true); };

                    newTab.Button.SetText("Task " + (taskIndex + 1));
                    int currTabIndex = taskIndex;

                    m_tabs[topicIndex].TaskTabs[taskIndex].CompletedAndCorrect = LabMgr.Instance.Stats.LabMap[m_currLab.ID].CompletionState[topicIndex][taskIndex];
                    if (m_tabs[topicIndex].TaskTabs[taskIndex].CompletedAndCorrect)
                    {
                        newTab.ShowCompletionSprite();
                    }
                    else
                    {
                        newTab.HideCompletionSprite();
                        allTasksInTopicCorrect = false;
                    }

                    newTab.Button.OnButtonPressed += delegate { HandleLabTabPressed(currTopicIndex, currTabIndex); };

                    LabTaskFrame newFrame = PopulateLabTaskFrame(m_currLab.Topics[topicIndex].Tasks[currTabIndex], currTopicIndex, m_tabs[topicIndex].TaskTabs[taskIndex].CompletedAndCorrect);
                    // m_frames[topicIndex].Add(newFrame);

                    RectTransform tabRect = newTabObj.GetComponent<RectTransform>();
                    float tabBuffer = m_xSpacing - tabRect.sizeDelta.x * tabRect.localScale.x;
                    // m_tabContainer.sizeDelta = new Vector2((startingOffset + tabRect.sizeDelta.x * tabRect.localScale.x + tabBuffer) * m_currLab.Topics[topicIndex].Tasks.Count, m_tabContainer.sizeDelta.y);

                    if (newFrame != null)
                    {
                        newTab.RegisterFrame(newFrame);
                    }

                    // Remove "Next" button from final lab task
                    if (taskIndex == m_currLab.Topics[topicIndex].Tasks.Count - 1 && topicIndex == m_currLab.Topics.Count - 1)
                    {
                        newFrame.ConfigureAsLastTask();
                    }

                    // tabs start hidden
                    DeactivateTopicTab(topicIndex);
                }

                if (allTasksInTopicCorrect)
                {
                    newTopicTab.ShowCompletionSprite();
                }
                else
                {
                    newTopicTab.HideCompletionSprite();
                }
            }

            // Open first tab
            if (m_tabs.Count > 0) {
                HandleLabTopicTabPressed(0, false);
            }

            m_ScrollVerticalValidVisibleIndex = 0;
            m_ScrollHorizontalValidVisibleIndex = 0;

            Vector3 newVertPos = m_ScrollVerticalContainer.transform.localPosition;
            newVertPos.y = m_VerticalScrollOrigin;
            m_ScrollVerticalContainer.transform.localPosition = newVertPos;

            Vector3 newHorizPos = m_ScrollHorizontalContainer.transform.localPosition;
            newHorizPos.x = m_HorizontalScrollOrigin;
            m_ScrollHorizontalContainer.transform.localPosition = newHorizPos;
        }

        public override void Close() {
            this.gameObject.SetActive(false);

            m_wordBankRoutine.Stop();

            // Remove button listeners
            m_ScrollUpBtn.OnButtonPressed -= HandleScrollUp;
            m_ScrollDownBtn.OnButtonPressed -= HandleScrollDown;
            m_ScrollLeftBtn.OnButtonPressed -= HandleScrollLeft;
            m_ScrollRightBtn.OnButtonPressed -= HandleScrollRight;

            ResetWorldMods();

            for (int i = 0; i < m_tabs.Count; i++) {
                for (int j = 0; j < m_tabs[i].TaskTabs.Count; j++) {
                    m_tabs[i].TaskTabs[j].RemoveCompletionStateListeners();
                    Destroy(m_tabs[i].TaskTabs[j].gameObject);
                }
                Destroy(m_tabs[i].gameObject);
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

        private LabTaskFrame PopulateLabTaskFrame(TaskInfo taskInfo, int currTopicIndex, bool completed) {
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
                case TaskType.Text:
                    // populate text
                    newFrameObj = Instantiate(m_taskFramePrefabText, this.transform);
                    newFrame = newFrameObj.GetComponent<LabTaskFrame>();
                    evaluables = newFrame.GetEvaluables();
                    for (int i = 0; i < evaluables.Length; i++)
                    {
                        TextHub textHub = evaluables[i].GetComponent<TextHub>();

                        TextTaskDefinition newDef = new TextTaskDefinition();
                        newDef.MainText = taskInfo.TextOnly;

                        textHub.SetDefinition(newDef);
                    }

                    break;
                default:
                    framePopulated = false;
                    break;
            }

            newFrame.LoadCompleted(completed);

            if (framePopulated) {
                newFrameObj.SetActive(false);
                m_frames[currTopicIndex].Add(newFrame);
            }

            return newFrame;
        }

        private void ApplyWorldMods(TaskInfo mods) {
            World.Instance.ModMgr.SetActiveMods(mods);

            // Tools
            World.Instance.ModMgr.SetAllowedTools(mods.AllowedTools);

            // Sets
            // only override the previous state if this task explicitly sets a new state
            if (!mods.Sets.IsEmpty())
            {
                World.Instance.ModMgr.ApplySets();
            }

            // Limits
            World.Instance.ModMgr.SetLimits(mods.Limits);

            // Trail
            World.Instance.ModMgr.SetTrail(mods.TrailSettings);

            // Tutorial
            World.Instance.ModMgr.SetTutorial(mods.Tutorial);

            // prevent move ball functionality
            if (mods.GrabAllowed)
            {
                World.Instance?.ModMgr.EnableGraphBallInteractions();
            }
            else
            {
                World.Instance?.ModMgr.DisableGraphBallInteractions();
            }
        }

        private void ResetWorldMods() {
            // Tools
            World.Instance?.ModMgr.ResetToolRestrictions();

            // Sets
            // N/A

            // Limits
            World.Instance?.ModMgr.ResetLimits();

            // Tutorial
            World.Instance?.ModMgr.ResetTutorial();

            // Restore move ball functionality
            World.Instance?.ModMgr.EnableGraphBallInteractions();

            World.Instance.ModMgr.DeactivateMods();
        }

        #region Handlers

        private void HandlePreActivateLab(Tuple<LabInfo, int> info) {
            m_currLab = info.Item1;
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

        private void HandleLabTopicTabPressed(int newTopicIndex, bool fromPlayerAction) {
            
            if (m_activeTopicIndex != -1)
            {
                foreach (var tab in m_tabs[m_activeTopicIndex].TaskTabs)
                {
                    tab.EnableCollider();
                    tab.UIButton.interactable = true;
                }
            }

            if (m_activeTopicIndex == -1) {
                // no tab activated yet; activate new tab
                ActivateTopicTab(newTopicIndex, fromPlayerAction);
            }
            else if (m_activeTopicIndex == newTopicIndex) {
                // same as current tab; no change needed
            }
            else {
                // hide existing tab
                DeactivateTopicTab(m_activeTopicIndex);

                // activate new tab
                ActivateTopicTab(newTopicIndex, fromPlayerAction);
            }

            m_ScrollHorizontalValidVisibleIndex = 0;
            m_activeTabIndex = 0;
            m_activeTopicIndex = newTopicIndex;

            RefreshInteractableTopicTabs();
            RefreshInteractableTaskTabs();
        }

        private void HandleTaskResetPressed() {
            ApplyWorldMods(m_currLab.Topics[m_activeTopicIndex].Tasks[m_activeTabIndex]);
        }

        private void HandleWordBankOpened()
        {
            foreach (var tab in m_tabs[m_activeTopicIndex].TaskTabs)
            {
                tab.DisableCollider();
                tab.UIButton.interactable = false;
            }
        }

        private void HandleWordBankClosed(string str)
        {
            // wait for seconds
            m_wordBankRoutine.Replace(OnCloseWordBankRoutine());
        }

        private IEnumerator OnCloseWordBankRoutine()
        {
            // wait a beat
            yield return 0.75f;

            foreach (var tab in m_tabs[m_activeTopicIndex].TaskTabs)
            {
                tab.EnableCollider();
                tab.UIButton.interactable = true;
            }
        }


        private void HandleHomeButtonPressed(object sender, EventArgs args)
        {
            if (GameMgr.I.AudioEnabled)
            {
                Tablet.Instance.PlayUIAudio(m_homeAudioClip);
            }

            EventMgr.Events?.Dispatch(GameEvents.DeactivateLab);
            EventMgr.Events?.Dispatch(GameEvents.ClickLabHome);
        }

        private void HandleScrollUp(object sender, EventArgs args)
        {
            if (m_ScrollVerticalContainer.transform.localPosition.y <= m_VerticalScrollOrigin || Mathf.Approximately(m_ScrollVerticalContainer.transform.localPosition.x, m_VerticalScrollOrigin)) { return; }
            m_ScrollVerticalContainer.transform.localPosition = m_ScrollVerticalContainer.transform.localPosition - new Vector3(0, m_VerticalScrollSpacing, 0);
            m_ScrollVerticalValidVisibleIndex--;

            RefreshInteractableTopicTabs();
            RefreshInteractableTaskTabs();

            EventMgr.Events.Dispatch(GameEvents.ClickSectionScrollUp);
        }

        private void HandleScrollDown(object sender, EventArgs args)
        {
            float comparePos = ((m_currLab.Topics.Count - SCROLL_VERTICAL_NUM) * m_VerticalScrollSpacing) + m_VerticalScrollOrigin;
            if (m_ScrollVerticalContainer.transform.localPosition.y >= comparePos || Mathf.Approximately(m_ScrollVerticalContainer.transform.localPosition.x, comparePos)) { return; }
            m_ScrollVerticalContainer.transform.localPosition = m_ScrollVerticalContainer.transform.localPosition + new Vector3(0, m_VerticalScrollSpacing, 0);
            m_ScrollVerticalValidVisibleIndex++;

            RefreshInteractableTopicTabs();
            RefreshInteractableTaskTabs();

            EventMgr.Events.Dispatch(GameEvents.ClickSectionScrollDown);
        }

        private void HandleScrollLeft(object sender, EventArgs args)
        {
            if (m_ScrollHorizontalContainer.transform.localPosition.x >= m_HorizontalScrollOrigin || Mathf.Approximately(m_ScrollHorizontalContainer.transform.localPosition.x, m_HorizontalScrollOrigin)) { return; }
            m_ScrollHorizontalContainer.transform.localPosition = m_ScrollHorizontalContainer.transform.localPosition + new Vector3(m_HorizontalScrollSpacing, 0, 0);
            m_ScrollHorizontalValidVisibleIndex--;

            RefreshInteractableTaskTabs();

            EventMgr.Events.Dispatch(GameEvents.ClickTaskScrollLeft);
        }

        private void HandleScrollRight(object sender, EventArgs args)
        {
            float comparePos = (-(m_currLab.Topics[m_activeTopicIndex].Tasks.Count - SCROLL_HORIZONTAL_NUM) * m_HorizontalScrollSpacing) + m_HorizontalScrollOrigin;
            if (m_ScrollHorizontalContainer.transform.localPosition.x <= comparePos || Mathf.Approximately(m_ScrollHorizontalContainer.transform.localPosition.x, comparePos)) { return; }
            m_ScrollHorizontalContainer.transform.localPosition = m_ScrollHorizontalContainer.transform.localPosition - new Vector3(m_HorizontalScrollSpacing, 0, 0);
            m_ScrollHorizontalValidVisibleIndex++;

            RefreshInteractableTaskTabs();

            EventMgr.Events.Dispatch(GameEvents.ClickTaskScrollRight);
        }

        private void HandleLabTaskCompletionUpdated(int topicIndex, bool fromReset)
        {
            LabStats currStats;

            bool allCorrect = true;
            for (int i = 0; i < m_tabs[topicIndex].TaskTabs.Count; i++)
            {
                if (!m_tabs[topicIndex].TaskTabs[i].CompletedAndCorrect)
                {
                    allCorrect = false;
                }

                currStats = LabMgr.Instance.Stats.LabMap[m_currLab.ID];
                currStats.CompletionState[topicIndex][i] = m_tabs[topicIndex].TaskTabs[i].CompletedAndCorrect;
                LabMgr.Instance.Stats.LabMap[m_currLab.ID] = currStats;
            }

            if (allCorrect)
            {
                // Display green checkmark
                m_tabs[topicIndex].ShowCompletionSprite();
                EventMgr.Events.Dispatch(GameEvents.SectionCompleted);
            }
            else
            {
                // Remove green checkmark
                m_tabs[topicIndex].HideCompletionSprite();
            }

            if (fromReset)
            {
                EventMgr.Events.Dispatch(GameEvents.ClickResetQuiz);
            }

            // Refresh Player Progress
            currStats = LabMgr.Instance.Stats.LabMap[m_currLab.ID];
            currStats.RefreshProgress();
            LabMgr.Instance.Stats.LabMap[m_currLab.ID] = currStats;
            EventMgr.Events.Dispatch(GameEvents.LabProgressUpdated);

            if (!fromReset)
            {
                EventMgr.Events.Dispatch(GameEvents.ClickSubmitAnswer);
            }

            if (currStats.Progress == 1)
            {
                EventMgr.Events.Dispatch(GameEvents.LabCompleted);
            }
        }

        private void HandleTaskNextPressed()
        {
            if (m_activeTabIndex == m_tabs[m_activeTopicIndex].TaskTabs.Count - 1)
            {
                if (m_activeTopicIndex == m_tabs.Count - 1)
                {
                    // at end of lab; do nothing
                }
                else
                {
                    // at end of topic; move to next topic
                    HandleLabTopicTabPressed(m_activeTopicIndex + 1, false);

                    // shift left task bar down if there is room below
                    if (m_ScrollDownBtn.isActiveAndEnabled)
                    {
                        HandleScrollDown(this, EventArgs.Empty);
                    }
                }
            }
            else
            {
                HandleLabTabPressed(m_activeTopicIndex, m_activeTabIndex + 1);

                // shift top task bar right if there is room on the right
                if (m_ScrollRightBtn.isActiveAndEnabled)
                {
                    HandleScrollRight(this, EventArgs.Empty);
                }
            }
        }

        #endregion // Handlers

        private void ActivateTab(int topicIndex, int taskIndex)
        {
            // TODO: highlight active tab

            m_frames[topicIndex][taskIndex].gameObject.SetActive(true);

            m_tabs[topicIndex].TaskTabs[taskIndex].ButtonImage.sprite = GameDB.Instance.LabTaskTabActive;
            m_tabs[topicIndex].TaskTabs[taskIndex].ButtonRect.sizeDelta = new Vector2(m_tabs[topicIndex].TaskTabs[taskIndex].ButtonRect.sizeDelta.x, TAB_HEIGHT_ACTIVE);

            EventMgr.Events.Dispatch(GameEvents.TaskSwitched, taskIndex);
            EventMgr.Events.Dispatch(GameEvents.ClickSelectTask, m_currLab.Topics[topicIndex].Tasks[taskIndex]);

            //if (!m_tabs[topicIndex].TaskTabs[taskIndex].HasBeenEvaluated())
            //{
                ApplyWorldMods(m_currLab.Topics[topicIndex].Tasks[taskIndex]);
            //}

            EventMgr.Events.Dispatch(GameEvents.TaskAssigned, m_currLab.Topics[topicIndex].Tasks[taskIndex]);
        }

        private void DeactivateTab(int topicIndex, int taskIndex)
        {
            // TODO: un-highlight prev active tab

            m_frames[topicIndex][taskIndex].gameObject.SetActive(false);

            m_tabs[topicIndex].TaskTabs[taskIndex].ButtonImage.sprite = GameDB.Instance.LabTaskTabInactive;
            m_tabs[topicIndex].TaskTabs[taskIndex].ButtonRect.sizeDelta = new Vector2(m_tabs[topicIndex].TaskTabs[taskIndex].ButtonRect.sizeDelta.x, TAB_HEIGHT_INACTIVE);
        }

        private void ActivateTopicTab(int topicIndex, bool fromPlayerAction)
        {
            int newTaskIndex = 0;

            // TODO: highlight active topic tab
            m_topicHeader.SetText(m_currLab.Topics[topicIndex].TopicHeader);

            for (int i = 0; i < m_tabs[topicIndex].TaskTabs.Count; i++)
            {
                m_tabs[topicIndex].TaskTabs[i].gameObject.SetActive(true);
            }

            m_frames[topicIndex][newTaskIndex].gameObject.SetActive(true);

            m_tabs[topicIndex].ButtonImage.sprite = GameDB.Instance.LabTopicTabActive;
            m_tabs[topicIndex].ButtonImage.SetNativeSize();
            m_tabs[topicIndex].TaskTabs[newTaskIndex].ButtonImage.sprite = GameDB.Instance.LabTaskTabActive;
            m_tabs[topicIndex].TaskTabs[newTaskIndex].ButtonRect.sizeDelta = new Vector2(m_tabs[topicIndex].TaskTabs[newTaskIndex].ButtonRect.sizeDelta.x, TAB_HEIGHT_ACTIVE);

            m_ScrollHorizontalContainer.localPosition = new Vector3(m_HorizontalScrollOrigin, m_ScrollHorizontalContainer.localPosition.y, m_ScrollHorizontalContainer.localPosition.z);

            EventMgr.Events.Dispatch(GameEvents.TaskSwitched, newTaskIndex);
            EventMgr.Events.Dispatch(GameEvents.SectionSwitched, topicIndex);

            if (fromPlayerAction)
            {
                EventMgr.Events.Dispatch(GameEvents.ClickSelectSection, m_currLab.Topics[topicIndex]);
            }

            //if (!m_tabs[topicIndex].TaskTabs[0].HasBeenEvaluated())
            //{
            ApplyWorldMods(m_currLab.Topics[topicIndex].Tasks[newTaskIndex]);
            //}
        }

        private void DeactivateTopicTab(int topicIndex)
        {
            // TODO: un-highlight prev active topic tab
            m_topicHeader.SetText("");

            for (int i = 0; i < m_tabs[topicIndex].TaskTabs.Count; i++)
            {
                m_frames[topicIndex][i].gameObject.SetActive(false);
                m_tabs[topicIndex].TaskTabs[i].gameObject.SetActive(false);

                m_tabs[topicIndex].TaskTabs[i].ButtonImage.sprite = GameDB.Instance.LabTaskTabInactive;
                m_tabs[topicIndex].TaskTabs[i].ButtonRect.sizeDelta = new Vector2(m_tabs[topicIndex].TaskTabs[i].ButtonRect.sizeDelta.x, TAB_HEIGHT_INACTIVE);
            }

            m_tabs[topicIndex].ButtonImage.sprite = GameDB.Instance.LabTopicTabInactive;
            m_tabs[topicIndex].ButtonImage.SetNativeSize();
        }

        private void RefreshInteractableTopicTabs()
        {
            m_visibleSections.Clear();

            // sections
            for (int i = 0; i < m_tabs.Count; i++)
            {
                if (i >= m_ScrollVerticalValidVisibleIndex && i <= m_ScrollVerticalValidVisibleIndex + SCROLL_VERTICAL_NUM - 1)
                {
                    // not masked
                    m_tabs[i].EnableCollider();
                    m_visibleSections.Add(new IndexedTopicInfo(i, m_currLab.Topics[m_activeTopicIndex]));
                }
                else
                {
                    m_tabs[i].DisableCollider();
                }
            }

            EventMgr.Events.Dispatch(GameEvents.SectionListDisplayed, m_visibleSections);
        }

        private void RefreshInteractableTaskTabs()
        {
            m_visibleTasks.Clear();

            // tasks
            for (int i = 0; i < m_tabs[m_activeTopicIndex].TaskTabs.Count; i++)
            {
                if (i >= m_ScrollHorizontalValidVisibleIndex && i <= m_ScrollHorizontalValidVisibleIndex + SCROLL_HORIZONTAL_NUM - 1)
                {
                    // not masked
                    m_tabs[m_activeTopicIndex].TaskTabs[i].EnableCollider();
                    
                    m_visibleTasks.Add(new IndexedTaskInfo(i, m_currLab.Topics[m_activeTopicIndex].Tasks[i]));
                }
                else
                {
                    m_tabs[m_activeTopicIndex].TaskTabs[i].DisableCollider();
                }
            }

            // show/hide relevant buttons
            if (m_currLab.Topics[m_activeTopicIndex].Tasks.Count > SCROLL_HORIZONTAL_NUM) {
                if (m_ScrollHorizontalValidVisibleIndex > 0)
                {
                    // tasks to the left
                    m_ScrollLeftBtn.gameObject.SetActive(true);
                }
                else
                {
                    m_ScrollLeftBtn.gameObject.SetActive(false);
                }

                if (m_ScrollHorizontalValidVisibleIndex + SCROLL_HORIZONTAL_NUM < m_currLab.Topics[m_activeTopicIndex].Tasks.Count)
                {
                    m_ScrollRightBtn.gameObject.SetActive(true);
                }
                else
                {
                    m_ScrollRightBtn.gameObject.SetActive(false);
                }
            }
            else
            {
                m_ScrollLeftBtn.gameObject.SetActive(false);
                m_ScrollRightBtn.gameObject.SetActive(false);
            }

            if (m_currLab.Topics.Count > SCROLL_VERTICAL_NUM)
            {
                if (m_ScrollVerticalValidVisibleIndex > 0)
                {
                    // tasks above
                    m_ScrollUpBtn.gameObject.SetActive(true);
                }
                else
                {
                    m_ScrollUpBtn.gameObject.SetActive(false);
                }

                if (m_ScrollVerticalValidVisibleIndex + SCROLL_VERTICAL_NUM < m_currLab.Topics.Count)
                {
                    m_ScrollDownBtn.gameObject.SetActive(true);
                }
                else
                {
                    m_ScrollDownBtn.gameObject.SetActive(false);
                }
            }
            else
            {
                m_ScrollDownBtn.gameObject.SetActive(false);
                m_ScrollUpBtn.gameObject.SetActive(false);
            }

            EventMgr.Events.Dispatch(GameEvents.TaskListDisplayed, m_visibleTasks);
        }

    }

}