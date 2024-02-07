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

        private enum Hand {
            LEFT,
            RIGHT
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

        #endregion // Logging Enums & Structs

        #region Logging Variables

        private OGDLog m_Log;

        private GamePlatform m_Platform;

        [NonSerialized] private bool m_Debug;


        #endregion // Logging Variables

        private void Awake() {
            Initialize();
        }

        #region IService

        protected void Initialize()
        {
            // General Events
            //Game.Events.Register(GameEvents.StoryEvalBegin, OnFeedbackBegin, this)
            //    .Register<string>(GameEvents.ProfileStarting, SetUserCode, this)

            
            // Analytics Events
            // text click
            //Game.Events.Register(GameEvents.TextClicked, LogTextClick, this)
            // display text dialog
                
            // SceneHelper.OnSceneLoaded += LogSceneChanged;

            // CrashHandler.OnCrash += OnCrash;

            // NetworkStats.OnError.Register(OnNetworkError);

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

            LogStartGame();
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

        private void LogClickNewGame(Hand inHand)
        {
            Debug.Log("[Analytics] event: click_new_game" + "\n" + "hand: " + inHand);

            using (var e = m_Log.NewEvent("click_new_game"))
            {
                e.Param("hand", inHand.ToString());
            }
        }

        private void LogClickResetSim(Hand inHand, StateProperties resetToProperties)
        {
            Debug.Log("[Analytics] event: click_reset_sim" + "\n" + "hand: " + inHand);

            using (var e = m_Log.NewEvent("click_new_game"))
            {
                e.Param("hand", inHand.ToString());
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
        
        private void LogGrabTablet(Vector3 startPos, Quaternion startRot, Hand inHand)
        {
            Debug.Log("[Analytics] event: grab_tablet");

            using (var e = m_Log.NewEvent("grab_tablet"))
            {
                e.Param("start_pos", JsonConvert.SerializeObject(startPos));
                e.Param("start_rot", JsonConvert.SerializeObject(startRot));
                e.Param("hand", inHand.ToString());
            }
        }

        private void LogReleaseTablet(Vector3 endPos, Quaternion endRot, Hand inHand)
        {
            Debug.Log("[Analytics] event: release_tablet");

            using (var e = m_Log.NewEvent("release_tablet"))
            {
                e.Param("end_pos", JsonConvert.SerializeObject(endPos));
                e.Param("end_rot", JsonConvert.SerializeObject(endRot));
                e.Param("hand", inHand.ToString());
            }
        }

        private void LogGrabWorkstationHandle(Vector3 startPos, Quaternion startRot, Hand inHand)
        {
            Debug.Log("[Analytics] event: grab_workstation_handle");

            using (var e = m_Log.NewEvent("grab_workstation_handle"))
            {
                e.Param("start_pos", JsonConvert.SerializeObject(startPos));
                e.Param("start_rot", JsonConvert.SerializeObject(startRot));
                e.Param("hand", inHand.ToString());
            }
        }

        private void LogReleaseWorkstationHandle(Vector3 endPos, Quaternion endRot, Hand inHand)
        {
            Debug.Log("[Analytics] event: release_workstation_handle");

            using (var e = m_Log.NewEvent("release_workstation_handle"))
            {
                e.Param("end_pos", JsonConvert.SerializeObject(endPos));
                e.Param("end_rot", JsonConvert.SerializeObject(endRot));
                e.Param("hand", inHand.ToString());
            }
        }

        private void LogClickRotateGraphCW(int startDegrees, int endDegrees, Hand inHand)
        {
            Debug.Log("[Analytics] event: click_rotate_graph_cw");

            using (var e = m_Log.NewEvent("click_rotate_graph_cw"))
            {
                e.Param("start_degrees", startDegrees);
                e.Param("end_degrees", endDegrees);
                e.Param("hand", inHand.ToString());
            }
        }

        private void LogClickRotateGraphCCW(int startDegrees, int endDegrees, Hand inHand)
        {
            Debug.Log("[Analytics] event: click_rotate_graph_ccw");

            using (var e = m_Log.NewEvent("click_rotate_graph_ccw"))
            {
                e.Param("start_degrees", startDegrees);
                e.Param("end_degrees", endDegrees);
                e.Param("hand", inHand.ToString());
            }
        }

        private void LogGrabGraphBall(Hand inHand)
        {
            Debug.Log("[Analytics] event: grab_graph_ball");

            using (var e = m_Log.NewEvent("grab_graph_ball"))
            {
                e.Param("hand", inHand.ToString());
            }
        }

        private void LogReleaseGraphBall(Hand inHand)
        {
            Debug.Log("[Analytics] event: release_graph_ball");

            using (var e = m_Log.NewEvent("release_graph_ball"))
            {
                e.Param("hand", inHand.ToString());
            }
        }

        private void LogClickToolToggle(ToolType toolName, bool toolReset, Hand inHand, bool toolEnabled)
        {
            Debug.Log("[Analytics] event: click_tool_toggle");

            using (var e = m_Log.NewEvent("click_tool_toggle"))
            {
                e.Param("tool_name", toolName.ToString());
                e.Param("tool_reset", toolReset);
                e.Param("hand", inHand.ToString());
                e.Param("enabled", toolEnabled);
            }
        }

        private void LogClickToolIncrease(ToolType toolName, float endVal, Hand inHand)
        {
            Debug.Log("[Analytics] event: click_tool_increase");

            using (var e = m_Log.NewEvent("click_tool_increase"))
            {
                e.Param("tool_name", toolName.ToString());
                e.Param("end_value", endVal);
                e.Param("hand", inHand.ToString());
            }
        }

        private void LogClickToolDecrease(ToolType toolName, float endVal, Hand inHand)
        {
            Debug.Log("[Analytics] event: click_tool_decrease");

            using (var e = m_Log.NewEvent("click_tool_decrease"))
            {
                e.Param("tool_name", toolName.ToString());
                e.Param("end_value", endVal);
                e.Param("hand", inHand.ToString());
            }
        }

        private void LogGrabToolSlider(ToolType toolName, float startVal, Hand inHand)
        {
            Debug.Log("[Analytics] event: grab_tool_slider");

            using (var e = m_Log.NewEvent("grab_tool_slider"))
            {
                e.Param("tool_name", toolName.ToString());
                e.Param("start_val", startVal);
                e.Param("hand", inHand.ToString());
            }
        }

        // release_tool_slider { tool_name, end_value // in physical units, not 0-1, hand : enum(LEFT, RIGHT), auto_release : bool // true if slider was automatically released due to hand getting to far away, or sim was reset }
        private void LogReleaseToolSlider(ToolType toolName, float endVal, Hand inHand, bool autoRelease)
        {
            Debug.Log("[Analytics] event: release_tool_slider");

            using (var e = m_Log.NewEvent("release_tool_slider"))
            {
                e.Param("tool_name", toolName.ToString());
                e.Param("end_value", endVal);
                e.Param("hand", inHand.ToString());
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

        private void LogClickSandboxMode(Hand inHand)
        {
            Debug.Log("[Analytics] event: click_sandbox_mode");

            using (var e = m_Log.NewEvent("click_sandbox_mode"))
            {
                e.Param("hand", inHand.ToString());
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

        private void LogClickLabMode(LabLogData lab, Hand inHand)
        {
            Debug.Log("[Analytics] event: click_lab_mode");

            using (var e = m_Log.NewEvent("click_lab_mode"))
            {
                e.Param("initial_lab", JsonConvert.SerializeObject(lab));
                e.Param("hand", inHand.ToString());
            }
        }

        private void LogClickLabScrollUp(Hand inHand)
        {
            Debug.Log("[Analytics] event: click_lab_scroll_up");

            using (var e = m_Log.NewEvent("click_lab_scroll_up"))
            {
                e.Param("hand", inHand.ToString());
            }
        }

        private void LogClickLabScrollDown(Hand inHand)
        {
            Debug.Log("[Analytics] event: click_lab_scroll_down");

            using (var e = m_Log.NewEvent("click_lab_scroll_down"))
            {
                e.Param("hand", inHand.ToString());
            }
        }

        // lab_menu_displayed { available_labs : array[Lab] // Lab objects won’t contain section array here }
        private void LogLabMenuDisplayed(List<LabLogData> availableLabs)
        {
            Debug.Log("[Analytics] event: lab_menu_displayed");

            using (var e = m_Log.NewEvent("lab_menu_displayed"))
            {
                e.Param("available_labs", JsonConvert.SerializeObject(availableLabs));
            }
        }

        private void LogSelectLab(string labName, Hand inHand)
        {
            Debug.Log("[Analytics] event: select_lab");

            using (var e = m_Log.NewEvent("select_lab"))
            {
                e.Param("lab_name", labName);
                e.Param("hand", inHand.ToString());
            }
        }

        private void LogClickLabHome(Hand inHand)
        {
            Debug.Log("[Analytics] event: click_lab_home");

            using (var e = m_Log.NewEvent("click_lab_home"))
            {
                e.Param("hand", inHand.ToString());
            }
        }

        private void LogClickSelectTask(Hand inHand, TaskLogData task)
        {
            Debug.Log("[Analytics] event: click_select_task");

            using (var e = m_Log.NewEvent("click_select_task"))
            {
                e.Param("hand", inHand.ToString());
                e.Param("task", JsonConvert.SerializeObject(task));
            }
        }

        private void LogClickTaskScrollLeft(Hand inHand)
        {
            Debug.Log("[Analytics] event: click_task_scroll_left");

            using (var e = m_Log.NewEvent("click_task_scroll_left"))
            {
                e.Param("hand", inHand.ToString());
            }
        }

        private void LogClickTaskScrollRight(Hand inHand)
        {
            Debug.Log("[Analytics] event: click_task_scroll_right");

            using (var e = m_Log.NewEvent("click_task_scroll_right"))
            {
                e.Param("hand", inHand.ToString());
            }
        }

        private void LogTaskListDisplayed(List<TaskLogData> tasks)
        {
            Debug.Log("[Analytics] event: task_list_displayed");

            using (var e = m_Log.NewEvent("task_list_displayed"))
            {
                e.Param("task_list", JsonConvert.SerializeObject(tasks));
            }
        }

        private void LogClickSelectSection(Hand inHand, SectionLogData section)
        {
            Debug.Log("[Analytics] event: click_select_section");

            using (var e = m_Log.NewEvent("click_select_section"))
            {
                e.Param("hand", inHand.ToString());
                e.Param("section", JsonConvert.SerializeObject(section));
            }
        }

        private void LogClickTaskScrollUp(Hand inHand)
        {
            Debug.Log("[Analytics] event: click_section_scroll_up");

            using (var e = m_Log.NewEvent("click_section_scroll_up"))
            {
                e.Param("hand", inHand.ToString());
            }
        }

        private void LogClickTaskScrollDown(Hand inHand)
        {
            Debug.Log("[Analytics] event: click_section_scroll_down");

            using (var e = m_Log.NewEvent("click_section_scroll_down"))
            {
                e.Param("hand", inHand.ToString());
            }
        }

        // section_list_displayed { available_sections : array[Section] // Section objects won’t contain task array here }
        private void LogSectionListDisplayed(List<SectionLogData> sections)
        {
            Debug.Log("[Analytics] event: section_list_displayed");

            using (var e = m_Log.NewEvent("section_list_displayed"))
            {
                e.Param("available_selections", JsonConvert.SerializeObject(sections));
            }
        }

        private void LogTargetStateAchieved(Dictionary<string, float> targetState)
        {
            Debug.Log("[Analytics] event: target_state_achieved");

            using (var e = m_Log.NewEvent("target_state_achieved"))
            {
                e.Param("target_state", JsonConvert.SerializeObject(targetState));
            }
        }


        private void LogTargetStateLost(Dictionary<string, float> targetState, List<string> incorrectVars)
        {
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

        private void LogClickSelectAnswer(QuizTaskLogData quizTask, int selectionIndex, bool isCorrect)
        {
            Debug.Log("[Analytics] event: click_select_answer");

            using (var e = m_Log.NewEvent("click_select_answer"))
            {
                e.Param("quiz_task", JsonConvert.SerializeObject(quizTask));
                e.Param("selection_index", selectionIndex);
                e.Param("is_correct_answer", isCorrect);
            }
        }

        private void LogClickDeselectAnswer(QuizTaskLogData quizTask, int selectionIndex, bool isCorrect)
        {
            Debug.Log("[Analytics] event: click_deselect_answer");

            using (var e = m_Log.NewEvent("click_deselect_answer"))
            {
                e.Param("quiz_task", JsonConvert.SerializeObject(quizTask));
                e.Param("selection_index", selectionIndex);
                e.Param("is_correct_answer", isCorrect);
            }
        }

        private void LogClickSubmitAnswer(QuizTaskLogData quizTask, bool isCorrect)
        {
            Debug.Log("[Analytics] event: click_submit_answer");

            using (var e = m_Log.NewEvent("click_submit_answer"))
            {
                e.Param("quiz_task", JsonConvert.SerializeObject(quizTask));
                e.Param("is_correct_answer", isCorrect);
            }
        }

        private void LogClickResetQuiz(QuizTaskLogData quizTask, bool wasCorrect)
        {
            Debug.Log("[Analytics] event: click_reset_quiz");

            using (var e = m_Log.NewEvent("click_reset_quiz"))
            {
                e.Param("quiz_task", JsonConvert.SerializeObject(quizTask));
                e.Param("was_correct_answer", wasCorrect);
            }
        }

        private void LogClickOpenWordBank(QuizTaskLogData quizTask)
        {
            Debug.Log("[Analytics] event: click_open_word_bank");

            using (var e = m_Log.NewEvent("click_open_word_bank"))
            {
                e.Param("quiz_task", JsonConvert.SerializeObject(quizTask));
            }
        }

        private void LogWordBankDisplayed(List<string> words)
        {
            Debug.Log("[Analytics] event: word_bank_displayed");

            using (var e = m_Log.NewEvent("word_bank_displayed"))
            {
                e.Param("words", JsonConvert.SerializeObject(words));
            }
        }

        private void LogWordBankClosed(List<string> words, string selectedWord)
        {
            Debug.Log("[Analytics] event: word_bank_closed");

            using (var e = m_Log.NewEvent("word_bank_closed"))
            {
                e.Param("words", JsonConvert.SerializeObject(words));
                e.Param("selected_word", JsonConvert.SerializeObject(selectedWord));
            }
        }

        #endregion // Log Events

        #region Other Events



        #endregion // Other Events

        #region Helpers



        #endregion // Helpers
    }
}
