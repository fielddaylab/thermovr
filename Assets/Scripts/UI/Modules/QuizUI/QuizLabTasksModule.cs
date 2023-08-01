using System.Collections.Generic;
using ThermoVR.UI;
using UnityEngine;

namespace ThermoVR.Lab
{
    // Manages loading of task tabs
    public class QuizLabTasksModule : UIModule
    {
        [SerializeField] private GameObject m_tabContainer;
        [SerializeField] private GameObject m_tabPrefab;
        [SerializeField] private GameObject m_taskFramePrefabMC;
        [SerializeField] private GameObject m_taskFramePrefabWordBank;
        [SerializeField] private GameObject m_taskFramePrefabReachState;

        [SerializeField] private float m_xSpacing;

        private List<GameObject> m_tabs; // TODO: make these pools
        private List<LabTaskFrame> m_frames;
        private List<TaskInfo> m_tasks;

        private int m_activeTabIndex;

        #region IUIModule

        public override void Init() {
            base.Init();

            GameMgr.Events?.Register<Cartridge>(GameEvents.ActivateCartridge, HandleActivateCartridge);
            GameMgr.Events?.Register<Cartridge>(GameEvents.DeactivateCartridge, HandleDeactivateCartridge);

            m_tabs = new List<GameObject>();
            m_frames = new List<LabTaskFrame>();
            m_activeTabIndex = -1;
        }

        public override void Open() {
            base.Open();

            if (m_tasks == null) {
                return;
            }

            for (int i = 0; i < m_tasks.Count; i++) {
                GameObject newTab = Instantiate(m_tabPrefab, m_tabContainer.transform);
                newTab.transform.localPosition += new Vector3(m_xSpacing * i, 0, 0);
                m_tabs.Add(newTab);

                ThermoButton button = newTab.GetComponent<ThermoButton>();
                button.SetText("Task " + (i + 1));
                int tabIndex = i;
                button.OnButtonPressed += delegate { HandleLabTabPressed(tabIndex); };

                PopulateLabTaskFrame(m_tasks[tabIndex]);
            }
        }

        public override void Close() {
            base.Close();

            for (int i = 0; i < m_tabs.Count; i++) {
                Destroy(m_tabs[i]);
            }
            m_tabs.Clear();


            for (int i = 0; i < m_frames.Count; i++) {
                Destroy(m_frames[i].gameObject);
            }
            m_frames.Clear();
        }

        #endregion // IUIModule

        private void PopulateLabTaskFrame(TaskInfo taskInfo) {
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
                        newDef.CorrectID = taskInfo.CorrectID;

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
                        newDef.CorrectID = taskInfo.CorrectID;

                        wordBankHub.SetDefinition(newDef);
                    }

                    break;
                case TaskType.ReachState:
                    // populate reach state instructions
                    newFrameObj = Instantiate(m_taskFramePrefabReachState, this.transform);
                    newFrame = newFrameObj.GetComponent<LabTaskFrame>();
                    evaluables = newFrame.GetEvaluables();
                    for (int i = 0; i < evaluables.Length; i++) {
                        // Get Reach State Component
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
        }

        #region Handlers

        private void HandleActivateCartridge(Cartridge cartridge) {
            switch (cartridge.GetCartridgeType()) {
                case Cartridge.CartridgeType.Lab:
                    m_tasks = cartridge.GetInfo().Tasks;
                    break;
                case Cartridge.CartridgeType.Sandbox:
                    break;
                default:
                    break;
            }
        }

        private void HandleDeactivateCartridge(Cartridge cartridge) {
            m_tasks = null;
        }

        private void HandleLabTabPressed(int newTabIndex) {
            if (m_activeTabIndex == -1) {
                // no tab activated yet; activate new tab
                m_frames[newTabIndex].gameObject.SetActive(true);
            }
            else if (m_activeTabIndex == newTabIndex) {
                // same as current tab; no change needed
            }
            else {
                // hide existing tab
                m_frames[m_activeTabIndex].gameObject.SetActive(false);

                // activate new tab
                m_frames[newTabIndex].gameObject.SetActive(true);
            }

            m_activeTabIndex = newTabIndex;
        }

        #endregion // Handlers
    }

}