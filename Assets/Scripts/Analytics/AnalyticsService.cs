#if UNITY_EDITOR || DEVELOPMENT_BUILD
    #define DEVELOPMENT
#endif // UNITY_EDITOR || DEVELOPMENT_BUILD


using BeauUtil;
using System;
using System.Collections.Generic;
using UnityEngine;
using FieldDay;
using BeauUtil.Tags;
using BeauPools;
using ThermoVR.Lab;
using ThermoVR.Tools;
using Newtonsoft.Json;

namespace ThermoVR.Analytics
{

    public partial class AnalyticsService : MonoBehaviour
    {
        #region Inspector

        [SerializeField, Required] private string m_AppId = "THERMOVR";
        [SerializeField, Required] private string m_AppVersion = "1.0";
        [SerializeField] private FirebaseConsts m_Firebase = default(FirebaseConsts);

        #endregion // Inspector

        #region Logging Enums & Structs

        private enum GamePlatform
        {
            VR,
            WEB
        }

        public enum Hand {
            LEFT,
            RIGHT,
            WEB // desktop
        }

        private enum GraphElement
        {
            AXIS_NUMBERS,
            GRID_LINES,
            REGION_LABELS,
            AXIS_TRACKERS
        }

        [Serializable]
        public struct StateProperties
        {
            public string Region;
            public double P;
            public double V;
            public double T;
            public double u;
            public double s;
            public double h;
            public double x;
        }

        [Serializable]
        public struct HeadsetPos
        {
            public float pos_x;
            public float pos_y;
            public float pos_z;
            public float rot_x;
            public float rot_y;
            public float rot_z;
            public float rot_w;
        }

        [Serializable]
        public struct SliderSettings
        {
            public bool Enabled;
            public float SliderVal;
        }

        public struct SliderPanelLogData
        {
            public SliderSettings Insulation;
            public SliderSettings LowerStop;
            public SliderSettings UpperStop;
            public SliderSettings Weight;
            public SliderSettings NegativeWeight;
            public SliderSettings Heat;
            public SliderSettings Cooling;
            public SliderSettings ChamberPressure;
            public SliderSettings ChamberTemperature;
        }

        public struct LabLogData
        {
            public int Index;
            public string LabName;
            public string LabAuthor;
            public string PercentComplete;
            public bool IsActive;
            public List<SectionLogData> Sections;
        }

        public struct SectionLogData
        {
            public int Index;
            public string LabName; // Only set if not in sections array of a lab
            public string Description; // header
            public bool IsComplete;
            public bool IsActive;
            public List<TaskLogData> Tasks;
        }

        public struct AnswerSelectLogData
        {
            public uint SelectionIndex;
            public bool IsCorrectAnswer;

            public AnswerSelectLogData(uint index, bool correct)
            {
                SelectionIndex = index;
                IsCorrectAnswer = correct;
            }
        }


        #endregion // Logging Enums & Structs

        #region Logging Variables

        private OGDLog m_Log;

        private GamePlatform m_Platform;

        private LabInfo m_ActiveLabInfo;
        private int m_ActiveSectionIndex;
        private int m_ActiveTaskIndex;
        private Hand m_LastHandPress;

        private List<string> m_LastKnownWordBankStrs = new List<string>();

        [NonSerialized] private bool m_Debug;


        #endregion // Logging Variables

        #region Unity Callbacks

        private void Awake() {
            Initialize();
        }

        #endregion // Unity Callbacks

        #region IService

        protected void Initialize()
        {
            // General Events
            GameMgr.Events.Register<string>(GameEvents.NewNameGenerated, SetUserCode, this)
                .Register<bool>(GameEvents.HandStartPress, OnHandStartPress, this)
                .Register<LabInfo>(GameEvents.PreActivateLab, OnPreActivateLab, this)
                .Register<int>(GameEvents.SectionSwitched, OnSectionSwitched, this)
                .Register<int>(GameEvents.TaskSwitched, OnTaskSwitched, this)
                .Register<List<string>>(GameEvents.TaskChoiceSelected, OnTaskChoiceSelected);

            // Analytics Events
            GameMgr.Events.Register(GameEvents.StartGame, LogStartGame, this)
                .Register(GameEvents.SelectLab, LogSelectLab, this)
                .Register(GameEvents.ClickLabHome, LogClickLabHome, this)
                .Register<TaskInfo>(GameEvents.ClickSelectTask, LogClickSelectTask, this)
                .Register<TopicInfo>(GameEvents.ClickSelectSection, LogClickSelectSection, this)
                .Register(GameEvents.ClickLabScrollUp, LogClickLabScrollUp, this)
                .Register(GameEvents.ClickLabScrollDown, LogClickLabScrollDown, this)
                .Register(GameEvents.ClickSectionScrollUp, LogClickSectionScrollUp, this)
                .Register(GameEvents.ClickSectionScrollDown, LogClickSectionScrollDown, this)
                .Register(GameEvents.ClickTaskScrollLeft, LogClickTaskScrollLeft, this)
                .Register(GameEvents.ClickTaskScrollRight, LogClickTaskScrollRight, this)
                .Register(GameEvents.TargetStateReached, LogTargetStateAchieved, this)
                .Register<List<string>>(GameEvents.TargetStateLost, LogTargetStateLost, this)
                .Register<AnswerSelectLogData>(GameEvents.ClickSelectAnswer, LogClickSelectAnswer, this)
                .Register<AnswerSelectLogData>(GameEvents.ClickDeselectAnswer, LogClickDeselectAnswer, this)
                .Register(GameEvents.ClickSubmitAnswer, LogClickSubmitAnswer, this)
                .Register(GameEvents.ClickResetQuiz, LogClickResetQuiz, this)
                .Register(GameEvents.ClickOpenWordBank, LogClickOpenWordBank, this)
                .Register<List<string>>(GameEvents.WordBankDisplayed, LogWordBankDisplayed, this)
                .Register<string>(GameEvents.WordBankClosed, LogWordBankClosed, this)
                .Register<List<IndexedTaskInfo>>(GameEvents.TaskListDisplayed, LogTaskListDisplayed, this)
                .Register<List<IndexedTopicInfo>>(GameEvents.SectionListDisplayed, LogSectionListDisplayed, this)
                .Register<List<IndexedLabInfo>>(GameEvents.LabMenuDisplayed, LogLabMenuDisplayed, this);

            m_Log = new OGDLog(new OGDLogConsts() {
                AppId = m_AppId,
                AppVersion = m_AppVersion,
                ClientLogVersion = 1
            });
            m_Log.UseFirebase(m_Firebase);

            #if DEVELOPMENT
                m_Debug = true;
            #endif // DEVELOPMENT

            m_Log.SetDebug(m_Debug);

            m_Platform = GamePlatform.VR;
#if UNITY_WEBGL
            m_Platform = GamePlatform.WEB;
#elif UNITY_ANDROID
            m_Platform = GamePlatform.VR;
#endif
        }

        private void SetUserCode(string userCode)
        {
            Debug.Log("[Analytics] Setting user code: " + userCode);
            m_Log.Initialize(new OGDLogConsts() {
                AppId = m_AppId,
                AppVersion = m_AppVersion,
                ClientLogVersion = 1,
                AppBranch = BuildInfo.Branch()
            });
            m_Log.SetUserId(userCode);
        }

        protected void Shutdown()
        {
            // Game.Events?.DeregisterAll(this);
        }
        #endregion // IService

        #region GameState

        private void UpdateGameStateThermoProperties(StateProperties newProperties)
        {
            using (var gs = m_Log.OpenGameState())
            {
                // thermo attributes
                gs.Param("thermo_attributes", JsonConvert.SerializeObject(newProperties));
            }
        }

        private void UpdateGameStateHeadsetPos(HeadsetPos newPos)
        {
            // headset
            using (var gs = m_Log.OpenGameState())
            {
                // headset
                gs.Param("headset", JsonConvert.SerializeObject(newPos));
            }
        }

        private void UpdateGameStatePanelSettings(SliderPanelLogData updatedPanel)
        {
            // panel settings
            using (var gs = m_Log.OpenGameState())
            {
                gs.Param("slider_insulation", JsonConvert.SerializeObject(updatedPanel.Insulation));
                gs.Param("slider_lower_stop", JsonConvert.SerializeObject(updatedPanel.LowerStop));
                gs.Param("slider_upper_stop", JsonConvert.SerializeObject(updatedPanel.UpperStop));
                gs.Param("slider_weight", JsonConvert.SerializeObject(updatedPanel.Weight));
                gs.Param("slider_negative_weight", JsonConvert.SerializeObject(updatedPanel.NegativeWeight));
                gs.Param("slider_heat", JsonConvert.SerializeObject(updatedPanel.Heat));
                gs.Param("slider_cooling", JsonConvert.SerializeObject(updatedPanel.Cooling));
                gs.Param("slider_chamber_pressure", JsonConvert.SerializeObject(updatedPanel.ChamberPressure));
                gs.Param("slider_chamber_temperature", JsonConvert.SerializeObject(updatedPanel.ChamberTemperature));
            }
        }

        private void UpdateGameStateLab(LabLogData updatedLab, SectionLogData updatedSection, TaskLogData updatedTask)
        {
            // lab
            using (var gs = m_Log.OpenGameState())
            {
                gs.Param("current_lab", JsonConvert.SerializeObject(updatedLab));
                gs.Param("current_section", JsonConvert.SerializeObject(updatedSection));
                gs.Param("current_task", JsonConvert.SerializeObject(updatedTask));
            }
        }

        #endregion // GameState

        #region Log Events

        private void LogStartGame()
        {
            Debug.Log("[Analytics] event: start_game" + "\n" + "Platform: " + m_Platform);

            using (var e = m_Log.NewEvent("start_game"))
            {
                e.Param("mode", m_Platform.ToString());
            }
        }

        private void LogClickNewGame()
        {
            Debug.Log("[Analytics] event: click_new_game" + "\n" + "hand: " + m_LastHandPress);

            using (var e = m_Log.NewEvent("click_new_game"))
            {
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogClickResetSim(StateProperties resetToProperties)
        {
            Debug.Log("[Analytics] event: click_reset_sim" + "\n" + "hand: " + m_LastHandPress);

            using (var e = m_Log.NewEvent("click_new_game"))
            {
                e.Param("hand", m_LastHandPress.ToString());
                e.Param("default_state", JsonConvert.SerializeObject(resetToProperties));
            }
        }

        private void LogHeadsetOff() {
            Debug.Log("[Analytics] event: headset_off");

            m_Log.NewEvent("headset_off");
        }

        private void LogHeadsetOn()
        {
            Debug.Log("[Analytics] event: headset_on");

            m_Log.NewEvent("headset_on");
        }

        /*  TODO: 
            headset_data { array of ~30 frame samples, each has { pos_x, pos_y, pos_z, rot_x, rot_y, rot_z, rot_w } of headset position/rotation at each frame }
            left_hand_data { array of ~30 frame samples, each has { pos_x, pos_y, pos_z, rot_x, rot_y, rot_z, rot_w } of left hand position/rotation at each frame }
            right_hand_data { array of ~30 frame samples, each has { pos_x, pos_y, pos_z, rot_x, rot_y, rot_z, rot_w } of right hand position/rotation at each frame }
        */
        
        private void LogGrabTablet(Vector3 startPos, Quaternion startRot)
        {
            Debug.Log("[Analytics] event: grab_tablet");

            using (var e = m_Log.NewEvent("grab_tablet"))
            {
                e.Param("start_pos", JsonConvert.SerializeObject(startPos));
                e.Param("start_rot", JsonConvert.SerializeObject(startRot));
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogReleaseTablet(Vector3 endPos, Quaternion endRot)
        {
            Debug.Log("[Analytics] event: release_tablet");

            using (var e = m_Log.NewEvent("release_tablet"))
            {
                e.Param("end_pos", JsonConvert.SerializeObject(endPos));
                e.Param("end_rot", JsonConvert.SerializeObject(endRot));
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogGrabWorkstationHandle(Vector3 startPos, Quaternion startRot)
        {
            Debug.Log("[Analytics] event: grab_workstation_handle");

            using (var e = m_Log.NewEvent("grab_workstation_handle"))
            {
                e.Param("start_pos", JsonConvert.SerializeObject(startPos));
                e.Param("start_rot", JsonConvert.SerializeObject(startRot));
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogReleaseWorkstationHandle(Vector3 endPos, Quaternion endRot)
        {
            Debug.Log("[Analytics] event: release_workstation_handle");

            using (var e = m_Log.NewEvent("release_workstation_handle"))
            {
                e.Param("end_pos", JsonConvert.SerializeObject(endPos));
                e.Param("end_rot", JsonConvert.SerializeObject(endRot));
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogClickRotateGraphCW(int startDegrees, int endDegrees)
        {
            Debug.Log("[Analytics] event: click_rotate_graph_cw");

            using (var e = m_Log.NewEvent("click_rotate_graph_cw"))
            {
                e.Param("start_degrees", startDegrees);
                e.Param("end_degrees", endDegrees);
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogClickRotateGraphCCW(int startDegrees, int endDegrees)
        {
            Debug.Log("[Analytics] event: click_rotate_graph_ccw");

            using (var e = m_Log.NewEvent("click_rotate_graph_ccw"))
            {
                e.Param("start_degrees", startDegrees);
                e.Param("end_degrees", endDegrees);
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogGrabGraphBall()
        {
            Debug.Log("[Analytics] event: grab_graph_ball");

            using (var e = m_Log.NewEvent("grab_graph_ball"))
            {
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogReleaseGraphBall()
        {
            Debug.Log("[Analytics] event: release_graph_ball");

            using (var e = m_Log.NewEvent("release_graph_ball"))
            {
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogClickToolToggle(ToolType toolName, bool toolReset, bool toolEnabled)
        {
            Debug.Log("[Analytics] event: click_tool_toggle");

            using (var e = m_Log.NewEvent("click_tool_toggle"))
            {
                e.Param("tool_name", toolName.ToString());
                e.Param("tool_reset", toolReset);
                e.Param("hand", m_LastHandPress.ToString());
                e.Param("enabled", toolEnabled);
            }
        }

        private void LogClickToolIncrease(ToolType toolName, float endVal)
        {
            Debug.Log("[Analytics] event: click_tool_increase");

            using (var e = m_Log.NewEvent("click_tool_increase"))
            {
                e.Param("tool_name", toolName.ToString());
                e.Param("end_value", endVal);
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogClickToolDecrease(ToolType toolName, float endVal)
        {
            Debug.Log("[Analytics] event: click_tool_decrease");

            using (var e = m_Log.NewEvent("click_tool_decrease"))
            {
                e.Param("tool_name", toolName.ToString());
                e.Param("end_value", endVal);
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogGrabToolSlider(ToolType toolName, float startVal)
        {
            Debug.Log("[Analytics] event: grab_tool_slider");

            using (var e = m_Log.NewEvent("grab_tool_slider"))
            {
                e.Param("tool_name", toolName.ToString());
                e.Param("start_val", startVal);
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        // release_tool_slider { tool_name, end_value // in physical units, not 0-1, hand : enum(LEFT, RIGHT), auto_release : bool // true if slider was automatically released due to hand getting to far away, or sim was reset }
        private void LogReleaseToolSlider(ToolType toolName, float endVal, bool autoRelease)
        {
            Debug.Log("[Analytics] event: release_tool_slider");

            using (var e = m_Log.NewEvent("release_tool_slider"))
            {
                e.Param("tool_name", toolName.ToString());
                e.Param("end_value", endVal);
                e.Param("hand", m_LastHandPress.ToString());
                e.Param("auto_release", autoRelease);
            }
        }

        /*  TODO: 
            gaze_object_end { object : enum(TABLET, PISTON, GRAPH, CONTROLS), gaze_duration }
        */
        
        private void LogToggleSetting(GraphElement elementType, bool enabled)
        {
            Debug.Log("[Analytics] event: toggle_setting");

            using (var e = m_Log.NewEvent("toggle_setting"))
            {
                e.Param("setting", elementType.ToString());
                e.Param("enabled", enabled);
            }
        }

        private void LogClickSandboxMode()
        {
            Debug.Log("[Analytics] event: click_sandbox_mode");

            using (var e = m_Log.NewEvent("click_sandbox_mode"))
            {
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogToolLocked(ToolType tool)
        {
            Debug.Log("[Analytics] event: tool_locked");

            using (var e = m_Log.NewEvent("tool_locked"))
            {
                e.Param("tool", tool.ToString());
            }
        }


        private void LogToolUnlocked(ToolType tool)
        {
            Debug.Log("[Analytics] event: tool_unlocked");

            using (var e = m_Log.NewEvent("tool_unlocked"))
            {
                e.Param("tool", tool.ToString());
            }
        }

        private void LogClickLabMode(LabLogData lab)
        {
            Debug.Log("[Analytics] event: click_lab_mode");

            using (var e = m_Log.NewEvent("click_lab_mode"))
            {
                e.Param("initial_lab", JsonConvert.SerializeObject(lab));
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogClickLabScrollUp()
        {
            Debug.Log("[Analytics] event: click_lab_scroll_up");

            using (var e = m_Log.NewEvent("click_lab_scroll_up"))
            {
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogClickLabScrollDown()
        {
            Debug.Log("[Analytics] event: click_lab_scroll_down");

            using (var e = m_Log.NewEvent("click_lab_scroll_down"))
            {
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        // lab_menu_displayed { available_labs : array[Lab] // Lab objects won’t contain section array here }
        private void LogLabMenuDisplayed(List<IndexedLabInfo> visibleLabs)
        {
            List<LabLogData> visibleLabsData = new List<LabLogData>();

            for (int labIndex = 0; labIndex < visibleLabs.Count; labIndex++)
            {
                visibleLabsData.Add(LabInfoToLabLogData(visibleLabs[labIndex].Info, visibleLabs[labIndex].Index, false));
            }

            Debug.Log("[Analytics] event: lab_menu_displayed");

            using (var e = m_Log.NewEvent("lab_menu_displayed"))
            {
                e.Param("available_labs", JsonConvert.SerializeObject(visibleLabsData));
            }
        }

        private void LogSelectLab()
        {
            Debug.Log("[Analytics] event: select_lab" + "\n lab name: " + m_ActiveLabInfo.Name + " \n hand: " + m_LastHandPress);

            using (var e = m_Log.NewEvent("select_lab"))
            {
                e.Param("lab_name", m_ActiveLabInfo.Name);
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogClickLabHome()
        {
            Debug.Log("[Analytics] event: click_lab_home" + "\n hand: " + m_LastHandPress);

            using (var e = m_Log.NewEvent("click_lab_home"))
            {
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogClickSelectTask(TaskInfo task)
        {
            TaskLogData taskData = TaskInfoToTaskLogData(task, m_ActiveSectionIndex, m_ActiveTaskIndex, false);

            Debug.Log("[Analytics] event: click_select_task");

            using (var e = m_Log.NewEvent("click_select_task"))
            {
                e.Param("hand", m_LastHandPress.ToString());
                e.Param("task", JsonConvert.SerializeObject(taskData));
            }
        }

        private void LogClickTaskScrollLeft()
        {
            Debug.Log("[Analytics] event: click_task_scroll_left");

            using (var e = m_Log.NewEvent("click_task_scroll_left"))
            {
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogClickTaskScrollRight()
        {
            Debug.Log("[Analytics] event: click_task_scroll_right");

            using (var e = m_Log.NewEvent("click_task_scroll_right"))
            {
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogTaskListDisplayed(List<IndexedTaskInfo> tasks)
        {
            List<TaskLogData> tasksData = new List<TaskLogData>();
            foreach (IndexedTaskInfo task in tasks)
            {
                tasksData.Add(TaskInfoToTaskLogData(task.Info, m_ActiveSectionIndex, task.Index, false));
            }

            Debug.Log("[Analytics] event: task_list_displayed");

            using (var e = m_Log.NewEvent("task_list_displayed"))
            {
                e.Param("task_list", JsonConvert.SerializeObject(tasksData));
            }
        }

        private void LogClickSelectSection(TopicInfo section)
        {
            SectionLogData sectionData;

            List<TaskLogData> tasksData = new List<TaskLogData>();
            for (int taskIdx = 0; taskIdx < section.Tasks.Count; taskIdx++)
            {
                tasksData.Add(TaskInfoToTaskLogData(section.Tasks[taskIdx], m_ActiveSectionIndex, taskIdx, true));
            }

            sectionData = TopicInfoToSectionLogData(section, m_ActiveSectionIndex, tasksData, false);

            Debug.Log("[Analytics] event: click_select_section");

            using (var e = m_Log.NewEvent("click_select_section"))
            {
                e.Param("hand", m_LastHandPress.ToString());
                e.Param("section", JsonConvert.SerializeObject(sectionData));
            }
        }

        private void LogClickSectionScrollUp()
        {
            Debug.Log("[Analytics] event: click_section_scroll_up");

            using (var e = m_Log.NewEvent("click_section_scroll_up"))
            {
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogClickSectionScrollDown()
        {
            Debug.Log("[Analytics] event: click_section_scroll_down");

            using (var e = m_Log.NewEvent("click_section_scroll_down"))
            {
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        // section_list_displayed { available_sections : array[Section] // Section objects won’t contain task array here }
        private void LogSectionListDisplayed(List<IndexedTopicInfo> sections)
        {
            List<SectionLogData> sectionsData = new List<SectionLogData>();

            for (int sectionIdx = 0; sectionIdx < sections.Count; sectionIdx++)
            {
                List<TaskLogData> tasksData = null;
                // Section objects won’t contain task array here
                /*
                tasksData = new List<TaskLogData>();
                for (int taskIdx = 0; taskIdx < sections[sectionIdx].Info.Tasks.Count; taskIdx++)
                {
                    tasksData.Add(TaskInfoToTaskLogData(sections[sectionIdx].Info.Tasks[taskIdx], m_ActiveSectionIndex, taskIdx, true));
                }
                */
                sectionsData.Add(TopicInfoToSectionLogData(sections[sectionIdx].Info, sections[sectionIdx].Index, tasksData, false));
            }

            Debug.Log("[Analytics] event: section_list_displayed");

            using (var e = m_Log.NewEvent("section_list_displayed"))
            {
                e.Param("available_sections", JsonConvert.SerializeObject(sectionsData));
            }
        }

        private void LogTargetStateAchieved()
        {
            Dictionary<string, float> targetState = new Dictionary<string, float>();
            TaskInfo info = m_ActiveLabInfo.Topics[m_ActiveSectionIndex].Tasks[m_ActiveTaskIndex];
            foreach (var simTarget in info.Targets)
            {
                targetState.Add(simTarget.TargetID.ToString(), simTarget.TargetVal);
            }

            Debug.Log("[Analytics] event: target_state_achieved");

            using (var e = m_Log.NewEvent("target_state_achieved"))
            {
                e.Param("target_state", JsonConvert.SerializeObject(targetState));
            }
        }


        private void LogTargetStateLost(List<string> incorrectVars)
        {
            Dictionary<string, float> targetState = new Dictionary<string, float>();
            TaskInfo info = m_ActiveLabInfo.Topics[m_ActiveSectionIndex].Tasks[m_ActiveTaskIndex];
            foreach (var simTarget in info.Targets)
            {
                targetState.Add(simTarget.TargetID.ToString(), simTarget.TargetVal);
            }

            Debug.Log("[Analytics] event: target_state_lost");

            using (var e = m_Log.NewEvent("target_state_lost"))
            {
                e.Param("target_state", JsonConvert.SerializeObject(targetState));
                e.Param("incorrect_variables", JsonConvert.SerializeObject(incorrectVars));
            }
        }

        /* TODO: 
            constant_variable_achieved { TODO }
            constant_variable_lost { TODO }
        */

        private void LogClickSelectAnswer(AnswerSelectLogData answerSelectData)
        {
            TaskLogData taskData = TaskInfoToTaskLogData(m_ActiveLabInfo.Topics[m_ActiveSectionIndex].Tasks[m_ActiveTaskIndex], m_ActiveSectionIndex, m_ActiveTaskIndex, false);

            Debug.Log("[Analytics] event: click_select_answer");

            using (var e = m_Log.NewEvent("click_select_answer"))
            {
                e.Param("quiz_task", JsonConvert.SerializeObject(taskData));
                e.Param("selection_index", answerSelectData.SelectionIndex);
                e.Param("is_correct_answer", answerSelectData.IsCorrectAnswer);
            }
        }

        private void LogClickDeselectAnswer(AnswerSelectLogData answerSelectData)
        {
            TaskLogData taskData = TaskInfoToTaskLogData(m_ActiveLabInfo.Topics[m_ActiveSectionIndex].Tasks[m_ActiveTaskIndex], m_ActiveSectionIndex, m_ActiveTaskIndex, false);

            Debug.Log("[Analytics] event: click_deselect_answer");

            using (var e = m_Log.NewEvent("click_deselect_answer"))
            {
                e.Param("quiz_task", JsonConvert.SerializeObject(taskData));
                e.Param("selection_index", answerSelectData.SelectionIndex);
                e.Param("is_correct_answer", answerSelectData.IsCorrectAnswer);
            }
        }

        private void LogClickSubmitAnswer()
        {
            TaskLogData taskData = TaskInfoToTaskLogData(m_ActiveLabInfo.Topics[m_ActiveSectionIndex].Tasks[m_ActiveTaskIndex], m_ActiveSectionIndex, m_ActiveTaskIndex, false);

            bool isCorrect = LabMgr.Instance.Stats.LabMap[m_ActiveLabInfo.ID].CompletionState[m_ActiveSectionIndex][m_ActiveTaskIndex];

            Debug.Log("[Analytics] event: click_submit_answer");

            using (var e = m_Log.NewEvent("click_submit_answer"))
            {
                e.Param("quiz_task", JsonConvert.SerializeObject(taskData));
                e.Param("is_correct_answer", isCorrect);
            }
        }

        private void LogClickResetQuiz()
        {
            TaskLogData taskData = TaskInfoToTaskLogData(m_ActiveLabInfo.Topics[m_ActiveSectionIndex].Tasks[m_ActiveTaskIndex], m_ActiveSectionIndex, m_ActiveTaskIndex, false);

            bool wasCorrect = LabMgr.Instance.Stats.LabMap[m_ActiveLabInfo.ID].CompletionState[m_ActiveSectionIndex][m_ActiveTaskIndex];

            Debug.Log("[Analytics] event: click_reset_quiz");

            using (var e = m_Log.NewEvent("click_reset_quiz"))
            {
                e.Param("quiz_task", JsonConvert.SerializeObject(taskData));
                e.Param("was_correct_answer", wasCorrect);
            }
        }

        private void LogClickOpenWordBank()
        {
            TaskLogData taskData = TaskInfoToTaskLogData(m_ActiveLabInfo.Topics[m_ActiveSectionIndex].Tasks[m_ActiveTaskIndex], m_ActiveSectionIndex, m_ActiveTaskIndex, false);

            Debug.Log("[Analytics] event: click_open_word_bank");

            using (var e = m_Log.NewEvent("click_open_word_bank"))
            {
                e.Param("quiz_task", JsonConvert.SerializeObject(taskData));
            }
        }

        private void LogWordBankDisplayed(List<string> words)
        {
            m_LastKnownWordBankStrs = words;

            Debug.Log("[Analytics] event: word_bank_displayed");

            using (var e = m_Log.NewEvent("word_bank_displayed"))
            {
                e.Param("words", JsonConvert.SerializeObject(words));
            }
        }

        private void LogWordBankClosed(string selectedWord)
        {
            Debug.Log("[Analytics] event: word_bank_closed");

            using (var e = m_Log.NewEvent("word_bank_closed"))
            {
                e.Param("words", JsonConvert.SerializeObject(m_LastKnownWordBankStrs));
                e.Param("selected_word", JsonConvert.SerializeObject(selectedWord));
            }
        }

        #endregion // Log Events

        #region Other Events

        private void OnPreActivateLab(LabInfo lab)
        {
            m_ActiveLabInfo = lab;
        }

        private void OnHandStartPress(bool leftHand)
        {
            m_LastHandPress = leftHand ? Hand.LEFT : Hand.RIGHT;
        }

        private void OnSectionSwitched(int newSectionIndex)
        {
            m_ActiveSectionIndex = newSectionIndex;
        }

        private void OnTaskSwitched(int newTaskIndex)
        {
            m_ActiveTaskIndex = newTaskIndex;
        }

        private void OnTaskChoiceSelected(List<string> selections)
        {
            LabMgr.Instance.Stats.LabMap[m_ActiveLabInfo.ID].SelectedOptionsState[m_ActiveSectionIndex][m_ActiveTaskIndex] = selections;
        }

        #endregion // Other Events

        #region Helpers

        private TaskLogData TaskInfoToTaskLogData(TaskInfo info, int topicIndex, int taskIndex, bool partOfSection)
        {
            TaskCategory taskCategory = (TaskCategory)info.TaskType;
            if (taskCategory == TaskCategory.TARGET_STATE || taskCategory == TaskCategory.CONSTANT_VARIABLE)
            {
                // target task
                TargetTaskLogData taskData = new TargetTaskLogData();
                taskData.Category = taskCategory;
                taskData.LabName = m_ActiveLabInfo.Name;
                if (!partOfSection) { taskData.SectionIndex = topicIndex; }
                taskData.Index = taskIndex;
                taskData.IsActive = true;
                taskData.IsComplete = LabMgr.Instance.Stats.LabMap[m_ActiveLabInfo.ID].CompletionState[topicIndex][taskIndex];
                taskData.AvailableTools = m_ActiveLabInfo.Topics[topicIndex].Tasks[taskIndex].AllowedTools;
                taskData.Prompts = new List<string>() { m_ActiveLabInfo.Topics[topicIndex].Tasks[taskIndex].InitialConditions };
                for (int i = 0; i < m_ActiveLabInfo.Topics[topicIndex].Tasks[taskIndex].TaskQuestions.Count; i++)
                {
                    taskData.Prompts.Add(m_ActiveLabInfo.Topics[topicIndex].Tasks[taskIndex].TaskQuestions[i]);
                }

                if (taskCategory == TaskCategory.TARGET_STATE)
                {
                    taskData.TargetStateTarget = new Dictionary<string, float>();
                    foreach (var simTarget in info.Targets)
                    {
                        taskData.TargetStateTarget.Add(simTarget.TargetID.ToString(), simTarget.TargetVal);
                    }
                }
                else
                {
                    // TODO: Implement if we ever add constant variable tasks
                    // taskData.ConstantVariableTarget = array[str], the vars/tools that must be made constant
                }

                return taskData;
            }
            else
            {
                // quiz task
                QuizTaskLogData taskData = new QuizTaskLogData();
                taskData.Category = taskCategory;
                taskData.LabName = m_ActiveLabInfo.Name;
                if (!partOfSection) { taskData.SectionIndex = topicIndex; }
                taskData.Index = taskIndex;
                taskData.IsActive = true;
                taskData.IsComplete = LabMgr.Instance.Stats.LabMap[m_ActiveLabInfo.ID].CompletionState[topicIndex][taskIndex];
                taskData.AvailableTools = m_ActiveLabInfo.Topics[topicIndex].Tasks[taskIndex].AllowedTools;
                taskData.Prompts = new List<string>() { m_ActiveLabInfo.Topics[topicIndex].Tasks[taskIndex].InitialConditions };
                for (int i = 0; i < m_ActiveLabInfo.Topics[topicIndex].Tasks[taskIndex].TaskQuestions.Count; i++)
                {
                    taskData.Prompts.Add(m_ActiveLabInfo.Topics[topicIndex].Tasks[taskIndex].TaskQuestions[i]);
                }

                taskData.Options = info.SecondaryTexts;
                taskData.SelectedOptions = LabMgr.Instance.Stats.LabMap[m_ActiveLabInfo.ID].SelectedOptionsState[m_ActiveSectionIndex][m_ActiveTaskIndex];
                taskData.Answer = new List<string>();
                foreach(var id in info.CorrectIDs)
                {
                    taskData.Answer.Add(id.ToString());
                }

                return taskData;
            }
        }

        private SectionLogData TopicInfoToSectionLogData(TopicInfo info, int topicIndex, List<TaskLogData> taskData, bool partOfLab)
        {
            SectionLogData sectionData = new SectionLogData();
            sectionData.Index = topicIndex;
            if (!partOfLab) { sectionData.LabName = m_ActiveLabInfo.Name; }
            sectionData.Description = info.TopicHeader;
            
            bool allCorrect = true;
            for (int i = 0; i < LabMgr.Instance.Stats.LabMap[m_ActiveLabInfo.ID].CompletionState[topicIndex].Length; i++)
            {
                if (!LabMgr.Instance.Stats.LabMap[m_ActiveLabInfo.ID].CompletionState[topicIndex][i])
                {
                    allCorrect = false;
                }
            }
            sectionData.IsComplete = allCorrect;
            sectionData.IsActive = true;
            if (taskData != null) { sectionData.Tasks = taskData; }

            return sectionData;
        }

        private LabLogData LabInfoToLabLogData(LabInfo info, int labIndex, bool includeSections)
        {
            LabLogData labData = new LabLogData();

            labData.Index = labIndex;
            labData.LabName = info.Name;
            labData.LabAuthor = info.Author;
            labData.PercentComplete = LabMgr.Instance.Stats.LabMap[info.ID].ToString();
            labData.IsActive = info.ID == m_ActiveLabInfo.ID;
            if (includeSections) { labData.Sections = null; }  // TODO: sections

            return labData;
        }

        #endregion // Helpers
    }
}
