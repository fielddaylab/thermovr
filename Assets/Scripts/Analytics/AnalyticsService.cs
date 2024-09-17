#if UNITY_EDITOR || DEVELOPMENT_BUILD
    #define DEVELOPMENT
#endif // UNITY_EDITOR || DEVELOPMENT_BUILD


using BeauUtil;
using System;
using System.Collections.Generic;
using UnityEngine;
using BeauUtil.Tags;
using Newtonsoft.Json;
using ThermoVR.Lab;
using ThermoVR.Tools;
using ThermoVR.UI.GraphElements;
using ThermoVR.Controls;
using ThermoVR.State;
using OGD;
using ThermoVR.UI;
using FieldDay;

namespace ThermoVR.Analytics
{
    [DefaultExecutionOrder(-10)]
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
            DESKTOP
        }

        private enum Location
        {
            TITLE_SCREEN,
            TABLET
        }

        private enum GraphElement
        {
            AXIS_NUMBERS,
            GRID_LINES,
            REGION_LABELS,
            AXIS_TRACKERS
        }

        private enum LogToolType
        {
            INSULATION,
            LOWER_STOP,
            UPPER_STOP,
            INCREASE_WEIGHT,
            DECREASE_WEIGHT,
            HEAT,
            COOLING,
            CHAMBER_TERMPERATURE,
            CHAMBER_PRESSURE,
            UNKOWN
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

        [Serializable]
        public struct LabLogData
        {
            public int Index;
            public string LabName;
            public string LabAuthor;
            public float PercentComplete;
            public bool IsActive;
            public List<SectionLogData> Sections;
        }

        [Serializable]
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

        private LabInfo m_ActiveLabInfo;
        private int m_ActiveLabIndex;
        private LabLogData m_ActiveLabLogData;
        private int m_ActiveSectionIndex;
        private int m_ActiveTaskIndex;
        private Hand m_LastHandPress;
        private LogToolType m_LastInputProxyType;
        private LogToolType m_LastKnownSliderToolType;
        private Hand m_LastKnownSliderHand;

        private bool m_IsGameMode;
        private STuple<float, float, float> m_LastKnownGameModeTarget;
        private int m_LastKnownGameModeScore;

        private List<string> m_LastKnownWordBankStrs = new List<string>();

        [NonSerialized] private bool m_Debug;


        #endregion // Logging Variables

        #region GameStateVars

        private float m_GSElapsedTime;
        private StateProperties m_GSProperties;
        private PositionDataFrame m_GSHeadset;
        private SliderPanelLogData m_GSPanel;
        private LabLogData m_GSLab;
        private SectionLogData m_GSSection;
        private TaskLogData m_GSTask;
        private int m_GSPlayScore;
        private string m_GSPlatform;
        private string m_GSTabletMode;

        private JsonBuilder m_JsonBuilder = new JsonBuilder(1024);

        #endregion // GameStateVars

        #region Unity Callbacks

        private void Awake() {
            Initialize();
        }

        #endregion // Unity Callbacks

        #region IService

        protected void Initialize()
        {
            m_LastInputProxyType = LogToolType.UNKOWN;

            // General Events
            EventMgr.Events.Register<string>(GameEvents.NewNameGenerated, SetUserCode, this)
                .Register<Hand>(GameEvents.HandStartPress, OnHandStartPress, this)
                .Register<STuple<LabInfo, int>>(GameEvents.PreActivateLab, OnPreActivateLab, this)
                .Register<int>(GameEvents.SectionSwitched, OnSectionSwitched, this)
                .Register<int>(GameEvents.TaskSwitched, OnTaskSwitched, this)
                .Register<List<string>>(GameEvents.TaskChoiceSelected, OnTaskChoiceSelected)
                .Register<StateProperties>(GameEvents.StatePropertiesUpdated, OnStatePropertiesUpdated)
                .Register<PositionDataFrame>(GameEvents.HeadsetPosUpdated, OnHeadsetPosUpdated)
                .Register<SliderPanelLogData>(GameEvents.SliderPanelUpdated, OnSliderPanelUpdated)
                .Register(GameEvents.LabProgressUpdated, OnLabProgressUpdated)
                .Register<float>(GameEvents.ElapsedTimeUpdated, OnElapsedTimeUpdated)
                .Register<STuple<float, float, float>>(GameEvents.GameModeCompleteGenerateTarget, OnGameModeCompleteGenerateTarget)
                .Register<int>(GameEvents.GameModeScoreUpdated, OnGameModeScoreUpdated)
                .Register(GameEvents.GameModeStarted, OnGameModeStarted)
                .Register(GameEvents.GameModeExited, OnGameModeExited)
                .Register<UIID>(GameEvents.TabletModeSwitched, OnTabletModeSwitched)
            ;

            // Analytics Events
            EventMgr.Events.Register(GameEvents.StartGame, LogStartGame, this)
                .Register(GameEvents.StartSession, LogStartSession, this)
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
                //.Register(GameEvents.TargetStateTaskBegan, LogTargetStateTaskBegan, this)
                //.Register(GameEvents.TargetStateTaskEnded, LogTargetStateTaskEnded, this)
                .Register(GameEvents.TargetStateEntered, LogTargetStateEntered, this)
                .Register(GameEvents.TargetStateCompleted, LogTargetStateCompleted, this)
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
                .Register<List<IndexedLabInfo>>(GameEvents.LabMenuDisplayed, LogLabMenuDisplayed, this)
                .Register<StateProperties>(GameEvents.ResetSimClicked, LogClickResetSim, this)
                .Register<STuple<PositionDataFrame, Hand>>(GameEvents.TabletGrabbed, LogGrabTablet, this)
                .Register<STuple<PositionDataFrame, Hand>>(GameEvents.TabletReleased, LogReleaseTablet, this)
                .Register<STuple<PositionDataFrame, Hand>>(GameEvents.WorkspaceHandleGrabbed, LogGrabWorkstationHandle, this)
                .Register<STuple<PositionDataFrame, Hand>>(GameEvents.WorkspaceHandleReleased, LogReleaseWorkstationHandle, this)
                .Register<STuple<float, float>>(GameEvents.RotateGraphClickedCW, LogClickRotateGraphCW, this)
                .Register<STuple<float, float>>(GameEvents.RotateGraphClickedCCW, LogClickRotateGraphCCW, this)
                .Register<Hand>(GameEvents.GraphBallGrabbed, LogGrabGraphBall, this)
                .Register<Hand>(GameEvents.GraphBallReleased, LogReleaseGraphBall, this)
                .Register(GameEvents.SandboxModeClicked, LogClickSandboxMode, this)
                .Register(GameEvents.LabModeClicked, LogClickLabMode, this)
                .Register(GameEvents.TaskCompleted, LogCompleteTask, this)
                .Register(GameEvents.SectionCompleted, LogCompleteSection, this)
                .Register(GameEvents.LabCompleted, LogCompleteLab, this)
                .Register(GameEvents.HeadsetOn, LogHeadsetOn, this)
                .Register(GameEvents.HeadsetOff, LogHeadsetOff, this)
                .Register<STuple<ToolType, bool, bool, int>>(GameEvents.ToolTogglePressed, LogClickToolToggle, this)
                .Register<STuple<ToolType, float, int>>(GameEvents.ClickToolIncrease, LogClickToolIncrease, this)
                .Register<STuple<ToolType, float, int>>(GameEvents.ClickToolDecrease, LogClickToolDecrease, this)
                .Register<STuple<ToolType, int, List<double>>> (GameEvents.MoveToolSlider, LogMoveToolSlider, this)
                .Register<STuple<ToolType, float, Hand, bool, int>>(GameEvents.ReleaseToolSlider, LogReleaseToolSlider, this)
                .Register<STuple<ToolType, float, Hand, int>>(GameEvents.GrabToolSlider, LogGrabToolSlider, this)
                .Register<Tool>(GameEvents.AllowTool, LogToolUnlocked, this)
                .Register<Tool>(GameEvents.DisallowTool, LogToolLocked, this)
                .Register(GameEvents.SettingsViewClicked, LogClickViewSettings, this)
                .Register<GraphSettingUpdate>(GameEvents.UpdateGraphSetting, LogClickToggleSetting)
                .Register<STuple<GazeTargetType, float>>(GameEvents.GazeEnd, LogGazeObjectEnd)
                .Register<PositionDataFrame[]>(GameEvents.ViewportData, LogViewportData)
                .Register<PositionDataFrame[]>(GameEvents.LeftHandData, LogLeftHandData)
                .Register<PositionDataFrame[]>(GameEvents.RightHandData, LogRightHandData)
                .Register<SimStateDataFrame[]>(GameEvents.SimStateData, LogSimStateData)
                .Register<ToolType>(GameEvents.EditToolValStarted, LogClickEditToolVal, this)
                .Register<float>(GameEvents.ProxyInputSubmitted, LogSetToolVal, this)
                .Register<string>(GameEvents.SetInvalidToolVal, LogSetInvalidToolVal, this)
                .Register(GameEvents.CancelEditToolVal, LogCancelEditToolVal, this)
                .Register(GameEvents.ClickGameMode, LogClickGameMode, this)
                .Register(GameEvents.ClickGameStart, LogClickGameStart, this)
                .Register(GameEvents.ClickGameStop, LogClickGameStop, this)
                .Register(GameEvents.ClickGameScoreReset, LogClickGameScoreReset, this)
                .Register(GameEvents.NewGameTargetAssigned, LogNewGameTargetAssigned, this)
                .Register(GameEvents.EnterNudgeMode, LogEnterNudgeMode, this)
                .Register(GameEvents.ExitNudgeMode, LogExitNudgeMode, this)
                .Register(GameEvents.TitleScreenDisplayed, LogTitleScreenDisplayed, this)
                .Register(GameEvents.TitleScreenClosed, LogTitleScreenClosed, this)
                .Register(GameEvents.ClickCloseTitleScreen, LogClickCloseTitleScreen, this)
                .Register<bool>(GameEvents.ClickDisplayCredits, LogClickDisplayCredits, this)
                .Register<bool>(GameEvents.ClickCloseCredits, LogClickCloseCredits, this)
                .Register<GraphSettingUpdate>(GameEvents.ClickToggleTitleSetting, LogClickToggleTitleSetting, this)
                .Register(GameEvents.ClickConfigTab, LogClickConfigTab, this)
                .Register(GameEvents.ClickControlsTab, LogClickControlsTab, this)
                .Register(GameEvents.NudgeHintDisplayed, LogNudgeHintDisplayed, this)
                .Register(GameEvents.NudgeHintHidden, LogNudgeHintHidden, this)
                .Register<TaskInfo>(GameEvents.TaskAssigned, LogTaskAssigned, this)
                .Register(GameEvents.TaskNextPressed, LogClickNextTask, this)
                .Register(GameEvents.ClickClearHeatMeter, LogClickClearHeatMeter, this)
                .Register(GameEvents.ClickClearWorkMeter, LogClickClearWorkMeter, this)
            ;

            m_Log = new OGDLog(
                new OGDLogConsts()
                {
                    AppId = m_AppId,
                    AppVersion = m_AppVersion,
                    ClientLogVersion = 4
                },
                new OGDLog.MemoryConfig
                (
                    OGDLog.MemoryConfig.Default.EventParameterBufferSize * 4,
                    OGDLog.MemoryConfig.Default.GameStateBufferSize * 3,
                    OGDLog.MemoryConfig.Default.PlayerDataBufferSize * 2
                )
            );
            m_Log.UseFirebase(m_Firebase);

            #if DEVELOPMENT
                m_Debug = true;
            #endif // DEVELOPMENT

            m_Log.SetDebug(m_Debug);

            var schedulingConfig = OGDLog.SchedulingConfig.Default;
            schedulingConfig.FlushDelay = 2;
            m_Log.ConfigureScheduling(schedulingConfig);

            m_GSPlatform = GamePlatform.VR.ToString();
#if UNITY_WEBGL
            m_GSPlatform = GamePlatform.DESKTOP.ToString();
            m_LastHandPress = Hand.MOUSE;
#elif UNITY_ANDROID
            m_GSPlatform = GamePlatform.VR.ToString();
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

        private void LateUpdate()
        {
            UpdateAllGameState();
        }

        protected void Shutdown()
        {
            // Game.Events?.DeregisterAll(this);
        }
        #endregion // IService

        #region GameState

        private void UpdateGameStateThermoProperties(StateProperties newProperties)
        {
            m_GSProperties = newProperties;
        }

        private void UpdateGameStateHeadsetPos(PositionDataFrame newPos)
        {
            m_GSHeadset = newPos;
        }

        private void UpdateGameStatePanelSettings(SliderPanelLogData updatedPanel)
        {
            m_GSPanel = updatedPanel;
        }

        private void UpdateGameStateLab(LabLogData updatedLab, SectionLogData updatedSection, TaskLogData updatedTask)
        {
            m_GSLab = updatedLab;
            m_GSSection = updatedSection;
            m_GSTask = updatedTask;
        }

        private void UpdateGameStateSection(SectionLogData updatedSection, TaskLogData updatedTask)
        {
            m_GSSection = updatedSection;
            m_GSTask = updatedTask;
        }

        private void UpdateGameStateTask(TaskLogData updatedTask)
        {
            m_GSTask = updatedTask;
        }

        private void UpdateAllGameState()
        {
            try
            {
                using (var gs = m_Log.OpenGameState(m_JsonBuilder))
                {
                    gs.Field("seconds_from_launch", m_GSElapsedTime);

                    WriteStateProperties(gs, "thermo_attributes", m_GSProperties);
                    WritePositionDataFrame(gs, "headset", m_GSHeadset);

                    WriteSliderSettings(gs, "slider_insulation", m_GSPanel.Insulation);
                    WriteSliderSettings(gs, "slider_lower_stop", m_GSPanel.LowerStop);
                    WriteSliderSettings(gs, "slider_upper_stop", m_GSPanel.UpperStop);
                    WriteSliderSettings(gs, "slider_weight", m_GSPanel.Weight);
                    WriteSliderSettings(gs, "slider_negative_weight", m_GSPanel.NegativeWeight);
                    WriteSliderSettings(gs, "slider_heat", m_GSPanel.Heat);
                    WriteSliderSettings(gs, "slider_cooling", m_GSPanel.Cooling);
                    WriteSliderSettings(gs, "slider_chamber_pressure", m_GSPanel.ChamberPressure);
                    WriteSliderSettings(gs, "slider_chamber_temperature", m_GSPanel.ChamberTemperature);

                    WriteLabLogData(gs, "current_lab", m_GSLab);
                    WriteSectionLogData(gs, "current_section", m_GSSection);
                    WriteTaskLogData(gs, "current_task", m_GSTask);
                    
                    gs.Field("play_score", m_LastKnownGameModeScore);
                    gs.Field("tablet_mode", m_GSTabletMode);
                    gs.Field("platform", m_GSPlatform);
                }
            }
            catch
            {
                Debug.LogWarning("[Analytics] Unable to update Game State! This may occur when the buffer sizes are too small for the amount of data being transmitted.");
            }
        }

        #endregion // GameState

        #region Log Events

        private void LogStartGame()
        {
            Debug.Log("[Analytics] event: game_start");

            using (var e = m_Log.NewEvent("game_start"))
            {

            }
        }

        private void LogStartSession()
        {
            Debug.Log("[Analytics] event: session_start");

            using (var e = m_Log.NewEvent("session_start"))
            {

            }
        }

        /// <summary>
        /// We don't currently have a "new game" button in the interface. But if we add one,
        /// this is what would log it.
        /// </summary>
        private void LogClickNewGame()
        {
            Debug.Log("[Analytics] event: click_new_game" + "\n" + "hand: " + m_LastHandPress);

            using (var e = m_Log.NewEvent("click_new_game"))
            {
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogClickResetSim(StateProperties resetTo)
        {
            Debug.Log("[Analytics] event: click_reset_sim" + "\n" + "hand: " + m_LastHandPress);

            using (var e = m_Log.NewEvent("click_reset_sim", m_JsonBuilder))
            {
                e.Field("hand", m_LastHandPress.ToString());
                WriteStateProperties(e, "default_state", resetTo);
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

        // headset_data { array of ~30 frame samples, each has { pos_x, pos_y, pos_z, rot_x, rot_y, rot_z, rot_w } of headset position/rotation at each frame }
        private void LogViewportData(PositionDataFrame[] data)
        {
            Debug.Log("[Analytics] event: viewport_data");

            using (var e = m_Log.NewEvent("viewport_data", m_JsonBuilder))
            {
                e.BeginArray("data");
                foreach (var f in data) {
                    WritePositionDataFrame(e, f);
                }
                e.EndArray();
            }
        }

        // left_hand_data { array of ~30 frame samples, each has { pos_x, pos_y, pos_z, rot_x, rot_y, rot_z, rot_w } of left hand position/rotation at each frame }
        private void LogLeftHandData(PositionDataFrame[] data)
        {
            Debug.Log("[Analytics] event: left_hand_data");

            using (var e = m_Log.NewEvent("left_hand_data", m_JsonBuilder))
            {
                e.BeginArray("data");
                foreach (var f in data) {
                    WritePositionDataFrame(e, f);
                }
                e.EndArray();
            }
        }

        // right_hand_data { array of ~30 frame samples, each has { pos_x, pos_y, pos_z, rot_x, rot_y, rot_z, rot_w } of right hand position/rotation at each frame }
        private void LogRightHandData(PositionDataFrame[] data)
        {
            Debug.Log("[Analytics] event: right_hand_data");

            using (var e = m_Log.NewEvent("right_hand_data", m_JsonBuilder))
            {
                e.BeginArray("data");
                foreach (var f in data) {
                    WritePositionDataFrame(e, f);
                }
                e.EndArray();
            }
        }

        // simulation_data { array of ~30 frame samples, each has { P, V, T, u, s, h, x } of sim vars at each frame }
        private unsafe void LogSimStateData(SimStateDataFrame[] data)
        {
            // TODO: fix so that no more Out of Memory errors are thrown
            Debug.Log("[Analytics] event: simulation_data");

            using (var e = m_Log.NewEvent("simulation_data", m_JsonBuilder))
            {
                e.BeginArray("data");
                foreach (var f in data) {
                    WriteSimStateDataFrame(e, f);
                }
                e.EndArray();
            }
        }

        private unsafe void LogGrabTablet(STuple<PositionDataFrame, Hand> args)
        {
            Debug.Log("[Analytics] event: grab_tablet");

            using (var e = m_Log.NewEvent("grab_tablet", m_JsonBuilder))
            {
                DoPositionDataFrame(e, args.Item1, "start_pos", "start_rot");
                e.Field("hand", args.Item2.ToString());
            }
        }

        private unsafe void LogReleaseTablet(STuple<PositionDataFrame, Hand> args)
        {
            Debug.Log("[Analytics] event: release_tablet");

            using (var e = m_Log.NewEvent("release_tablet", m_JsonBuilder))
            {
                DoPositionDataFrame(e, args.Item1, "end_pos", "end_rot");
                e.Field("hand", args.Item2.ToString());
            }
        }

        private unsafe void LogGrabWorkstationHandle(STuple<PositionDataFrame, Hand> args)
        {
            Debug.Log("[Analytics] event: grab_workstation_handle");

            using (var e = m_Log.NewEvent("grab_workstation_handle", m_JsonBuilder))
            {
                DoPositionDataFrame(e, args.Item1, "start_pos", "start_rot");
                e.Field("hand", args.Item2.ToString());
            }
        }

        private unsafe void LogReleaseWorkstationHandle(STuple<PositionDataFrame, Hand> args)
        {
            Debug.Log("[Analytics] event: release_workstation_handle");

            using (var e = m_Log.NewEvent("release_workstation_handle", m_JsonBuilder))
            {
                DoPositionDataFrame(e, args.Item1, "end_pos", "end_rot");
                e.Field("hand", args.Item2.ToString());
            }
        }

        private void LogClickRotateGraphCW(STuple<float, float> degrees)
        {
            Debug.Log("[Analytics] event: click_rotate_graph_cw");

            using (var e = m_Log.NewEvent("click_rotate_graph_cw"))
            {
                e.Param("start_degrees", degrees.Item1);
                e.Param("end_degrees", degrees.Item2);
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogClickRotateGraphCCW(STuple<float, float> degrees)
        {
            Debug.Log("[Analytics] event: click_rotate_graph_ccw");

            using (var e = m_Log.NewEvent("click_rotate_graph_ccw"))
            {
                e.Param("start_degrees", degrees.Item1);
                e.Param("end_degrees", degrees.Item2);
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogGrabGraphBall(Hand handType)
        {
            Debug.Log("[Analytics] event: grab_graph_ball");

            using (var e = m_Log.NewEvent("grab_graph_ball"))
            {
                e.Param("hand", handType.ToString());
            }
        }

        private void LogReleaseGraphBall(Hand handType)
        {
            Debug.Log("[Analytics] event: release_graph_ball");

            using (var e = m_Log.NewEvent("release_graph_ball"))
            {
                e.Param("hand", handType.ToString());
            }
        }

        private void LogClickToolToggle(STuple<ToolType, bool, bool, int> args)
        {
            LogToolType type = ToolTypeToLogToolType(args.Item1, args.Item4);

            Debug.Log("[Analytics] event: click_tool_toggle");

            using (var e = m_Log.NewEvent("click_tool_toggle"))
            {
                e.Param("tool_name", type.ToString());
                e.Param("tool_enabled", args.Item2);
                e.Param("tool_reset", args.Item3);
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogClickToolIncrease(STuple<ToolType, float, int> args)
        {
            LogToolType type = ToolTypeToLogToolType(args.Item1, args.Item3);

            Debug.Log("[Analytics] event: click_tool_increase");

            using (var e = m_Log.NewEvent("click_tool_increase"))
            {
                e.Param("tool_name", type.ToString());
                e.Param("end_value", args.Item2);
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogClickToolDecrease(STuple<ToolType, float, int> args)
        {
            LogToolType type = ToolTypeToLogToolType(args.Item1, args.Item3);

            Debug.Log("[Analytics] event: click_tool_decrease");

            using (var e = m_Log.NewEvent("click_tool_decrease"))
            {
                e.Param("tool_name", type.ToString());
                e.Param("end_value", args.Item2);
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogGrabToolSlider(STuple<ToolType, float, Hand, int> args)
        {
            LogToolType type = ToolTypeToLogToolType(args.Item1, args.Item4);

            m_LastKnownSliderToolType = type;
            m_LastKnownSliderHand = args.Item3;

            Debug.Log("[Analytics] event: grab_tool_slider");

            using (var e = m_Log.NewEvent("grab_tool_slider"))
            {
                e.Param("tool_name", type.ToString());
                e.Param("start_val", args.Item2);
                e.Param("hand", args.Item3.ToString());
            }
        }

        private void LogMoveToolSlider(STuple<ToolType, int, List<double>> args) {
            LogToolType type = ToolTypeToLogToolType(args.Item1, args.Item2);

            Debug.Log("[Analytics] event: move_tool_slider");

            using (var e = m_Log.NewEvent("move_tool_slider", m_JsonBuilder)) {
                e.Field("tool_name", type.ToString());
                e.Field("hand", m_LastKnownSliderHand.ToString());
                e.BeginArray("movement_data");
                foreach (var item in args.Item3) {
                    e.Item(item);
                }
                e.EndArray();
            }

            //using (var e = m_Log.NewEvent("move_tool_slider")) {
            //    e.Param("tool_name", type.ToString());
            //    e.Param("hand", m_LastKnownSliderHand.ToString());
            //    e.Json("movement_data", JsonConvert.SerializeObject(args.Item3));
            //}
        }

        // release_tool_slider { tool_name, end_value // in physical units, not 0-1, hand : enum(LEFT, RIGHT), auto_release : bool // true if slider was automatically released due to hand getting to far away, or sim was reset }
        private void LogReleaseToolSlider(STuple<ToolType, float, Hand, bool, int> args)
        {
            LogToolType type = ToolTypeToLogToolType(args.Item1, args.Item5);

            m_LastKnownSliderToolType = LogToolType.UNKOWN;

            Debug.Log("[Analytics] event: release_tool_slider");

            using (var e = m_Log.NewEvent("release_tool_slider"))
            {
                e.Param("tool_name", type.ToString());
                e.Param("end_value", args.Item2);
                e.Param("hand", args.Item3.ToString());
                e.Param("auto_release", args.Item4);
            }
        }

        // gaze_object_end { object : enum(TABLET, PISTON, GRAPH, CONTROLS), gaze_duration }
        private void LogGazeObjectEnd(STuple<GazeTargetType, float> args)
        {
            Debug.Log("[Analytics] event: gaze_object_end");

            using (var e = m_Log.NewEvent("gaze_object_end"))
            {
                e.Param("object", args.Item1.ToString());
                e.Param("gaze_duration", args.Item2);
            }
        }

        private void LogClickViewSettings()
        {
            Debug.Log("[Analytics] event: click_view_settings");

            using (var e = m_Log.NewEvent("click_view_settings"))
            {
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogClickToggleSetting(GraphSettingUpdate settingUpdate)
        {
            Debug.Log("[Analytics] event: click_toggle_setting");

            using (var e = m_Log.NewEvent("click_toggle_setting"))
            {
                e.Param("setting", settingUpdate.GraphElementID.ToString());
                e.Param("enabled", settingUpdate.ToggleVal);
            }
        }

        private void LogClickToggleTitleSetting(GraphSettingUpdate settingUpdate)
        {
            Debug.Log("[Analytics] event: click_toggle_title_setting");

            using (var e = m_Log.NewEvent("click_toggle_title_setting"))
            {
                e.Param("setting", settingUpdate.GraphElementID.ToString());
                e.Param("enabled", settingUpdate.ToggleVal);
            }
        }

        private void LogToolLocked(Tool tool)
        {
            LogToolType type = ToolTypeToLogToolType(tool.tool_type, tool.unique_id);

            Debug.Log("[Analytics] event: tool_locked");

            using (var e = m_Log.NewEvent("tool_locked"))
            {
                e.Param("tool", type.ToString());
            }
        }


        private void LogToolUnlocked(Tool tool)
        {
            LogToolType type = ToolTypeToLogToolType(tool.tool_type, tool.unique_id);

            Debug.Log("[Analytics] event: tool_unlocked");

            using (var e = m_Log.NewEvent("tool_unlocked"))
            {
                e.Param("tool", type.ToString());
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

        private void LogClickLabMode()
        {
            Debug.Log("[Analytics] event: click_lab_mode");

            if (m_ActiveLabInfo.ID.Equals(string.Empty))
            {
                using (var e = m_Log.NewEvent("click_lab_mode"))
                {
                    e.Param("initial_lab", "null");
                    e.Param("hand", m_LastHandPress.ToString());
                }
            }
            else
            {
                LabLogData newLab = LabInfoToLabLogData(m_ActiveLabInfo, m_ActiveLabIndex, true);

                using (var e = m_Log.NewEvent("click_lab_mode"))
                {
                    e.Json("initial_lab", JsonConvert.SerializeObject(newLab));
                    e.Param("hand", m_LastHandPress.ToString());
                }
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
                e.Json("available_labs", JsonConvert.SerializeObject(visibleLabsData));
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

        private void LogClickSelectSection(TopicInfo section)
        {
            SectionLogData sectionData;

            sectionData = TopicInfoToSectionLogData(section, m_ActiveSectionIndex, false);

            Debug.Log("[Analytics] event: click_select_section");

            using (var e = m_Log.NewEvent("click_select_section"))
            {
                e.Param("hand", m_LastHandPress.ToString());
                e.Json("section", JsonConvert.SerializeObject(sectionData));
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
                sectionsData.Add(TopicInfoToSectionLogData(sections[sectionIdx].Info, sections[sectionIdx].Index, false, true));
            }

            Debug.Log("[Analytics] event: section_list_displayed");

            using (var e = m_Log.NewEvent("section_list_displayed"))
            {
                e.Json("available_sections", JsonConvert.SerializeObject(sectionsData));
            }
        }

        private void LogClickSelectTask(TaskInfo task)
        {
            TaskLogData taskData = TaskInfoToTaskLogData(task, m_ActiveSectionIndex, m_ActiveTaskIndex, false);

            Debug.Log("[Analytics] event: click_select_task");

            using (var e = m_Log.NewEvent("click_select_task"))
            {
                e.Param("hand", m_LastHandPress.ToString());
                e.Json("task", JsonConvert.SerializeObject(taskData));
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
                e.Json("task_list", JsonConvert.SerializeObject(tasksData));
            }
        }

        /*
        private void LogTargetStateTaskBegan()
        {
            Debug.Log("[Analytics] event: target_state_task_began");

            using (var e = m_Log.NewEvent("target_state_task_began"))
            {

            }
        }

        private void LogTargetStateTaskEnded()
        {
            Debug.Log("[Analytics] event: target_state_task_ended");

            using (var e = m_Log.NewEvent("target_state_task_ended"))
            {

            }
        }
        */

        private void LogTargetStateEntered()
        {
            TaskInfo info = m_ActiveLabInfo.Topics[m_ActiveSectionIndex].Tasks[m_ActiveTaskIndex];

            Debug.Log("[Analytics] event: target_state_entered");

            using (var e = m_Log.NewEvent("target_state_entered", m_JsonBuilder))
            {
                e.BeginObject("target_state");
                foreach (var simTarget in info.Targets) {
                    e.Field(simTarget.TargetID.ToString(), simTarget.TargetVal);
                }
                e.EndObject();
            }
        }

        private void LogTargetStateCompleted()
        {
            Dictionary<string, float> targetState = new Dictionary<string, float>();
            Dictionary<string, float> targetTolerances = new Dictionary<string, float>();

            if (!m_IsGameMode)
            {
                TaskInfo info = m_ActiveLabInfo.Topics[m_ActiveSectionIndex].Tasks[m_ActiveTaskIndex];
                foreach (var simTarget in info.Targets)
                {
                    targetState.Add(simTarget.TargetID.ToString(), simTarget.TargetVal);
                    targetTolerances.Add(simTarget.TargetID.ToString(), simTarget.TargetRange);
                }
            }


            Debug.Log("[Analytics] event: target_state_completed");

            using (var e = m_Log.NewEvent("target_state_completed"))
            {
                e.Json("target_state", JsonConvert.SerializeObject(targetState));
                e.Json("tolerances", JsonConvert.SerializeObject(targetTolerances));
                if (m_IsGameMode) {
                    e.Param("score_value", m_LastKnownGameModeScore);
                } else {
                    e.Param("score_value", "null");
                }
            }
        }


        private void LogTargetStateLost(List<string> incorrectVars)
        {
            TaskInfo info = m_ActiveLabInfo.Topics[m_ActiveSectionIndex].Tasks[m_ActiveTaskIndex];

            Debug.Log("[Analytics] event: target_state_lost");

            using (var e = m_Log.NewEvent("target_state_lost", m_JsonBuilder))
            {
                e.BeginObject("target_state");
                foreach (var simTarget in info.Targets) {
                    e.Field(simTarget.TargetID.ToString(), simTarget.TargetVal);
                }
                e.EndObject();

                e.BeginArray("incorrect_variables");
                foreach(var v in incorrectVars) {
                    e.Item(v);
                }
                e.EndArray();
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
                e.Json("quiz_task", JsonConvert.SerializeObject(taskData));
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
                e.Json("quiz_task", JsonConvert.SerializeObject(taskData));
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
                e.Json("quiz_task", JsonConvert.SerializeObject(taskData));
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
                e.Json("quiz_task", JsonConvert.SerializeObject(taskData));
                e.Param("was_correct_answer", wasCorrect);
            }
        }

        private void LogClickOpenWordBank()
        {
            TaskLogData taskData = TaskInfoToTaskLogData(m_ActiveLabInfo.Topics[m_ActiveSectionIndex].Tasks[m_ActiveTaskIndex], m_ActiveSectionIndex, m_ActiveTaskIndex, false);

            Debug.Log("[Analytics] event: click_open_word_bank");

            using (var e = m_Log.NewEvent("click_open_word_bank"))
            {
                e.Json("quiz_task", JsonConvert.SerializeObject(taskData));
            }
        }

        private void LogWordBankDisplayed(List<string> words)
        {
            m_LastKnownWordBankStrs = words;

            Debug.Log("[Analytics] event: word_bank_displayed");

            using (var e = m_Log.NewEvent("word_bank_displayed"))
            {
                e.Json("words", JsonConvert.SerializeObject(words));
            }
        }

        private void LogWordBankClosed(string selectedWord)
        {
            Debug.Log("[Analytics] event: word_bank_closed");

            using (var e = m_Log.NewEvent("word_bank_closed"))
            {
                e.Json("words", JsonConvert.SerializeObject(m_LastKnownWordBankStrs));
                e.Param("selected_word", selectedWord);
            }
        }

        private void LogCompleteLab()
        {
            LabLogData completedLab = LabInfoToLabLogData(m_ActiveLabInfo, m_ActiveLabIndex, false);
            Debug.Log("[Analytics] event: complete_lab");

            using (var e = m_Log.NewEvent("complete_lab"))
            {
                e.Json("lab", JsonConvert.SerializeObject(completedLab));
            }
        }

        private void LogCompleteSection()
        {
            TopicInfo section = m_ActiveLabInfo.Topics[m_ActiveSectionIndex];

            SectionLogData sectionData = TopicInfoToSectionLogData(section, m_ActiveSectionIndex, false);

            Debug.Log("[Analytics] event: complete_section");

            using (var e = m_Log.NewEvent("complete_section"))
            {
                e.Json("section", JsonConvert.SerializeObject(sectionData));
            }
        }

        private void LogCompleteTask()
        {
            var data = TaskInfoToTaskLogData(m_ActiveLabInfo.Topics[m_ActiveSectionIndex].Tasks[m_ActiveTaskIndex], m_ActiveSectionIndex, m_ActiveTaskIndex, false);

            Debug.Log("[Analytics] event: complete_task");

            using (var e = m_Log.NewEvent("complete_task"))
            {
                e.Json("task", JsonConvert.SerializeObject(data));
            }
        }

        private void LogClickEditToolVal(ToolType toolType)
        {
            Debug.Log("[Analytics] event: click_edit_tool_val");

            m_LastInputProxyType = ToolTypeToLogToolType(toolType, 0);

            using (var e = m_Log.NewEvent("click_edit_tool_val"))
            {
                e.Param("tool_name", m_LastInputProxyType.ToString());
                e.Param("hand", Hand.MOUSE.ToString());
            }
        }

        private void LogSetToolVal(float endVal)
        {
            Debug.Log("[Analytics] event: set_tool_val");

            using (var e = m_Log.NewEvent("set_tool_val"))
            {
                e.Param("tool_name", m_LastInputProxyType.ToString());
                e.Param("new_value", endVal);
                e.Param("auto_release", false);
                e.Param("hand", Hand.MOUSE.ToString());
            }
        }

        private void LogSetInvalidToolVal(string badVal)
        {
            Debug.Log("[Analytics] event: set_invalid_tool_val");

            using (var e = m_Log.NewEvent("set_invalid_tool_val"))
            {
                e.Param("tool_name", m_LastInputProxyType.ToString());
                e.Param("bad_value", badVal);
                e.Param("hand", Hand.MOUSE.ToString());
            }
        }

        private void LogCancelEditToolVal()
        {
            Debug.Log("[Analytics] event: cancel_edit_tool_val");

            using (var e = m_Log.NewEvent("cancel_edit_tool_val"))
            {
                e.Param("tool_name", m_LastInputProxyType.ToString());
                e.Param("hand", Hand.MOUSE.ToString());
            }
        }

        private void LogClickGameMode()
        {
            Debug.Log("[Analytics] event: click_game_mode");

            using (var e = m_Log.NewEvent("click_game_mode"))
            {
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogClickGameStart()
        {
            Debug.Log("[Analytics] event: click_game_start");

            using (var e = m_Log.NewEvent("click_game_start"))
            {
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogNewGameTargetAssigned()
        {
            Debug.Log("[Analytics] event: click_game_target_assigned");

            using (var e = m_Log.NewEvent("click_game_target_assigned"))
            {
                e.Param("p", m_LastKnownGameModeTarget.Item1);
                e.Param("v", m_LastKnownGameModeTarget.Item2);
                e.Param("t", m_LastKnownGameModeTarget.Item3);
            }
        }

        private void LogClickGameStop()
        {
            Debug.Log("[Analytics] event: click_game_stop");

            using (var e = m_Log.NewEvent("click_game_stop"))
            {
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogClickGameScoreReset()
        {
            Debug.Log("[Analytics] event: click_game_score_reset");

            using (var e = m_Log.NewEvent("click_game_score_reset"))
            {
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogEnterNudgeMode()
        {
            Debug.Log("[Analytics] event: enter_nudge_mode");

            using (var e = m_Log.NewEvent("enter_nudge_mode"))
            {
                e.Param("tool", m_LastKnownSliderToolType.ToString());
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogExitNudgeMode()
        {
            Debug.Log("[Analytics] event: exit_nudge_mode");

            using (var e = m_Log.NewEvent("exit_nudge_mode"))
            {
                e.Param("tool", m_LastKnownSliderToolType.ToString());
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogTitleScreenDisplayed()
        {
            Debug.Log("[Analytics] event: title_screen_displayed");

            using (var e = m_Log.NewEvent("title_screen_displayed"))
            {

            }
        }

        private void LogTitleScreenClosed()
        {
            Debug.Log("[Analytics] event: title_screen_closed");

            using (var e = m_Log.NewEvent("title_screen_closed"))
            {

            }
        }

        private void LogClickCloseTitleScreen()
        {
            Debug.Log("[Analytics] event: click_close_title_screen");

            using (var e = m_Log.NewEvent("click_close_title_screen"))
            {

            }
        }

        private void LogClickDisplayCredits(bool isTitle)
        {
            Debug.Log("[Analytics] event: click_display_credits");

            using (var e = m_Log.NewEvent("click_display_credits"))
            {
                e.Param("location", isTitle ? Location.TITLE_SCREEN.ToString() : Location.TABLET.ToString());
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogClickCloseCredits(bool isTitle)
        {
            Debug.Log("[Analytics] event: click_close_credits");

            using (var e = m_Log.NewEvent("click_close_credits"))
            {
                e.Param("location", isTitle ? Location.TITLE_SCREEN.ToString() : Location.TABLET.ToString());
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogClickConfigTab()
        {
            Debug.Log("[Analytics] event: click_config_tab");

            using (var e = m_Log.NewEvent("click_config_tab"))
            {

            }
        }

        private void LogClickControlsTab()
        {
            Debug.Log("[Analytics] event: click_controls_tab");

            using (var e = m_Log.NewEvent("click_controls_tab"))
            {

            }
        }

        private void LogNudgeHintDisplayed()
        {
            Debug.Log("[Analytics] event: nudge_hint_displayed");

            using (var e = m_Log.NewEvent("nudge_hint_displayed"))
            {

            }
        }

        private void LogNudgeHintHidden()
        {
            Debug.Log("[Analytics] event: nudge_hint_hidden");

            using (var e = m_Log.NewEvent("nudge_hint_hidden"))
            {

            }
        }

        private void LogTaskAssigned(TaskInfo task)
        {
            TaskLogData taskData = TaskInfoToTaskLogData(task, m_ActiveSectionIndex, m_ActiveTaskIndex, false);

            Debug.Log("[Analytics] event: click_task_assigned");

            using (var e = m_Log.NewEvent("click_task_assigned"))
            {
                e.Param("hand", m_LastHandPress.ToString());
                e.Json("task", JsonConvert.SerializeObject(taskData));
            }
        }

        private void LogClickNextTask()
        {
            Debug.Log("[Analytics] event: click_next_task");

            using (var e = m_Log.NewEvent("click_next_task"))
            {
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogClickClearHeatMeter()
        {
            Debug.Log("[Analytics] event: click_clear_heat_meter");

            using (var e = m_Log.NewEvent("click_clear_heat_meter"))
            {
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        private void LogClickClearWorkMeter()
        {
            Debug.Log("[Analytics] event: click_clear_work_meter");

            using (var e = m_Log.NewEvent("click_clear_work_meter"))
            {
                e.Param("hand", m_LastHandPress.ToString());
            }
        }

        #endregion // Log Events

        #region Other Events

        private void OnPreActivateLab(STuple<LabInfo, int> labInfo)
        {
            m_ActiveLabInfo = labInfo.Item1;
            m_ActiveLabIndex = labInfo.Item2;

            m_ActiveSectionIndex = 0;
            m_ActiveTaskIndex = 0;

            TopicInfo section = m_ActiveLabInfo.Topics[m_ActiveSectionIndex];
            LabLogData newLabData = LabInfoToLabLogData(m_ActiveLabInfo, m_ActiveLabIndex, false);
            SectionLogData newSectionData = TopicInfoToSectionLogData(section, m_ActiveSectionIndex, false, true);
            TaskLogData newTaskData = TaskInfoToTaskLogData(m_ActiveLabInfo.Topics[m_ActiveSectionIndex].Tasks[m_ActiveTaskIndex], m_ActiveSectionIndex, m_ActiveTaskIndex, false);
            UpdateGameStateLab(newLabData, newSectionData, newTaskData);
        }

        private void OnHandStartPress(Hand handType)
        {
            m_LastHandPress = handType;
        }

        private void OnSectionSwitched(int newSectionIndex)
        {
            m_ActiveSectionIndex = newSectionIndex;

            TopicInfo section = m_ActiveLabInfo.Topics[m_ActiveSectionIndex];
            SectionLogData newSectionData = TopicInfoToSectionLogData(section, m_ActiveSectionIndex, false, true);
            TaskLogData newTaskData = TaskInfoToTaskLogData(m_ActiveLabInfo.Topics[m_ActiveSectionIndex].Tasks[m_ActiveTaskIndex], m_ActiveSectionIndex, m_ActiveTaskIndex, false);
            UpdateGameStateSection(newSectionData, newTaskData);
        }

        private void OnTaskSwitched(int newTaskIndex)
        {
            m_ActiveTaskIndex = newTaskIndex;

            TaskLogData newData = TaskInfoToTaskLogData(m_ActiveLabInfo.Topics[m_ActiveSectionIndex].Tasks[m_ActiveTaskIndex], m_ActiveSectionIndex, m_ActiveTaskIndex, false);
            UpdateGameStateTask(newData);
        }

        private void OnTaskChoiceSelected(List<string> selections)
        {
            LabMgr.Instance.Stats.LabMap[m_ActiveLabInfo.ID].SelectedOptionsState[m_ActiveSectionIndex][m_ActiveTaskIndex] = selections;
        }

        private void OnStatePropertiesUpdated(StateProperties newProperties)
        {
            UpdateGameStateThermoProperties(newProperties);
        }

        private void OnHeadsetPosUpdated(PositionDataFrame newPos)
        {
            UpdateGameStateHeadsetPos(newPos);
        }

        private void OnSliderPanelUpdated(SliderPanelLogData newData)
        {
            UpdateGameStatePanelSettings(newData);
        }

        private void OnLabProgressUpdated()
        {
            m_GSLab.PercentComplete = LabMgr.Instance.Stats.LabMap[m_ActiveLabInfo.ID].Progress;
        }

        private void OnElapsedTimeUpdated(float newTime)
        {
            m_GSElapsedTime = newTime;
        }

        private void OnGameModeCompleteGenerateTarget(STuple<float, float, float> pvt)
        {
            m_LastKnownGameModeTarget = pvt;
        }

        private void OnGameModeScoreUpdated(int newScore)
        {
            m_LastKnownGameModeScore = newScore;
        }

        private void OnGameModeStarted()
        {
            m_IsGameMode = true;
        }

        private void OnGameModeExited()
        {
            m_IsGameMode = false;
        }

        private void OnTabletModeSwitched(UIID newMode)
        {
            switch (newMode) {
                case UIID.Sandbox:
                    m_GSTabletMode = TabletMode.SANDBOX.ToString();
                    break;
                case UIID.Lab:
                    m_GSTabletMode = TabletMode.LAB.ToString();
                    break;
                case UIID.Game:
                    m_GSTabletMode = TabletMode.GAME.ToString();
                    break;
                case UIID.Graph:
                    m_GSTabletMode = TabletMode.SETTINGS.ToString();
                    break;
                default:
                    break;
            }

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
                taskData.AvailableTools = AllowedToolsToLogAllowedTools(m_ActiveLabInfo.Topics[topicIndex].Tasks[taskIndex].AllowedTools);
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
                taskData.AvailableTools = AllowedToolsToLogAllowedTools(m_ActiveLabInfo.Topics[topicIndex].Tasks[taskIndex].AllowedTools);
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

        private List<string> AllowedToolsToLogAllowedTools(List<ToolType> allowedTools)
        {
            List<string> logList = new List<string>();
            if (allowedTools == null) { return logList; }

            foreach (ToolType tool in allowedTools)
            {
                if (tool == ToolType.Stops)
                {
                    logList.Add(ToolTypeToLogToolType(tool, 1).ToString());
                    logList.Add(ToolTypeToLogToolType(tool, 2).ToString());
                }
                else
                {
                    logList.Add(ToolTypeToLogToolType(tool, 0).ToString());
                }
            }

            return logList;
        }

        private SectionLogData TopicInfoToSectionLogData(TopicInfo info, int topicIndex, bool partOfLab, bool overrideTaskSections = false)
        {
            List<TaskLogData> taskData = new List<TaskLogData>();
            if (!overrideTaskSections)
            {
                for (int taskIdx = 0; taskIdx < info.Tasks.Count; taskIdx++)
                {
                    taskData.Add(TaskInfoToTaskLogData(info.Tasks[taskIdx], topicIndex, taskIdx, true));
                }
            }

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
            labData.PercentComplete = LabMgr.Instance.Stats.LabMap.ContainsKey(info.ID) ? LabMgr.Instance.Stats.LabMap[info.ID].Progress : 0;
            labData.IsActive = info.ID == m_ActiveLabInfo.ID;
            if (includeSections) {
                List<SectionLogData> allSections = new List<SectionLogData>();
                for (int s = 0; s < info.Topics.Count; s++)
                {
                    TopicInfo section = info.Topics[s];
                    SectionLogData newSectionData = TopicInfoToSectionLogData(section, s, false);
                    allSections.Add(newSectionData);
                }
                labData.Sections = allSections;
            }
            else { labData.Sections = null; }
            return labData;
        }

        private LogToolType ToolTypeToLogToolType(ToolType inType, int uniqueStopID)
        {
            switch(inType)
            {
                case ToolType.Burner:
                    return LogToolType.HEAT;
                case ToolType.Coil:
                    return LogToolType.COOLING;
                case ToolType.Insulator:
                    return LogToolType.INSULATION;
                case ToolType.NegativeWeight:
                    return LogToolType.DECREASE_WEIGHT;
                case ToolType.Weight:
                    return LogToolType.INCREASE_WEIGHT;
                case ToolType.Stops:
                    if (uniqueStopID == 1) { return LogToolType.UPPER_STOP; }
                    else if (uniqueStopID == 2) { return LogToolType.LOWER_STOP; }
                    else { return LogToolType.UNKOWN; }
                case ToolType.SurroundingPressure:
                    return LogToolType.CHAMBER_PRESSURE;
                case ToolType.SurroundingTemperature:
                    return LogToolType.CHAMBER_TERMPERATURE;
                default:
                    return LogToolType.UNKOWN;
            }
        }

        #endregion // Helpers

        #region Types

        static private unsafe void WritePositionDataFrame(OGDExtensions.JsonScope scope, in PositionDataFrame frame) {
            scope.BeginObject();
            DoPositionDataFrame(scope, frame, "pos", "rot");
            scope.EndObject();
        }

        static private unsafe void WritePositionDataFrame(OGDExtensions.JsonScope scope, string fieldName, in PositionDataFrame frame) {
            scope.BeginObject(fieldName);
            DoPositionDataFrame(scope, frame, "pos", "rot");
            scope.EndObject();
        }

        static private unsafe void DoPositionDataFrame(OGDExtensions.JsonScope scope, in PositionDataFrame frame, string posFieldName, string rotFieldName) {
            scope.BeginArray(posFieldName)
                .Item(frame.pos[0])
                .Item(frame.pos[1])
                .Item(frame.pos[2])
                .EndArray();
            scope.BeginArray(rotFieldName)
                .Item(frame.rot[0])
                .Item(frame.rot[1])
                .Item(frame.rot[2])
                .Item(frame.rot[3])
                .EndArray();
        }

        static private unsafe void WriteSimStateDataFrame(OGDExtensions.JsonScope scope, in SimStateDataFrame frame) {
            scope.BeginObject()
                .Field("P", frame.P)
                .Field("V", frame.V)
                .Field("T", frame.T)
                .Field("u", frame.u)
                .Field("s", frame.s)
                .Field("h", frame.h)
                .Field("x", frame.x)
                .EndObject();
        }

        static private unsafe void WriteStateProperties(OGDExtensions.JsonScope scope, string fieldName, in StateProperties properties) {
            scope.BeginObject(fieldName)
                .Field("Region", properties.Region)
                .Field("P", properties.P)
                .Field("V", properties.V)
                .Field("T", properties.T)
                .Field("u", properties.u)
                .Field("s", properties.s)
                .Field("h", properties.h)
                .Field("x", properties.x)
                .EndObject();
        }

        static private unsafe void WriteSliderSettings(OGDExtensions.JsonScope scope, string fieldName, in SliderSettings settings) {
            scope.BeginObject(fieldName)
                .Field("Enabled", settings.Enabled)
                .Field("SliderVal", settings.SliderVal)
                .EndObject();
        }

        static private unsafe void WriteLabLogData(OGDExtensions.JsonScope scope, string fieldName, in LabLogData data) {
            scope.BeginObject(fieldName);
            DoLabLogData(scope, data);
            scope.EndObject();
        }

        static private unsafe void DoLabLogData(OGDExtensions.JsonScope scope, in LabLogData data) {
            scope.Field("Index", data.Index)
                .Field("LabName", data.LabName)
                .Field("LabAuthor", data.LabAuthor)
                .Field("PercentComplete", data.PercentComplete, 3)
                .Field("IsActive", data.IsActive);

            scope.BeginArray("Sections");
            if (data.Sections != null) {
                foreach (var section in data.Sections) {
                    WriteSectionLogData(scope, section);
                }
            }
            scope.EndArray();
        }

        #region SectionLogData

        static private unsafe void WriteSectionLogData(OGDExtensions.JsonScope scope, in SectionLogData data) {
            scope.BeginObject();
            DoSectionLogData(scope, data);
            scope.EndObject();
        }

        static private unsafe void WriteSectionLogData(OGDExtensions.JsonScope scope, string fieldName, in SectionLogData data) {
            scope.BeginObject(fieldName);
            DoSectionLogData(scope, data);
            scope.EndObject();
        }

        static private unsafe void DoSectionLogData(OGDExtensions.JsonScope scope, in SectionLogData data) {
            scope.Field("Index", data.Index)
                .Field("LabName", data.LabName)
                .Field("Description", data.Description)
                .Field("IsComplete", data.IsComplete)
                .Field("IsActive", data.IsActive);

            scope.BeginArray("Tasks");
            if (data.Tasks != null) {
                foreach (var task in data.Tasks) {
                    WriteTaskLogData(scope, task);
                }
            }
            scope.EndArray();
        }

        #endregion // SectionLogData

        #region TaskLogData

        static private unsafe void WriteTaskLogData(OGDExtensions.JsonScope scope, string fieldName, in TaskLogData data) {
            if (data == null) {
                scope.Field(fieldName, (string) null);
            } else {
                scope.BeginObject(fieldName);
                DoTaskLogData(scope, data);
                scope.EndObject();
            }
        }

        static private unsafe void WriteTaskLogData(OGDExtensions.JsonScope scope, in TaskLogData data) {
            scope.BeginObject();
            DoTaskLogData(scope, data);
            scope.EndObject();
        }
        
        static private unsafe void DoTaskLogData(OGDExtensions.JsonScope scope, in TaskLogData data) {
            scope.Field("Category", data.Category.ToString())
                .Field("LabName", data.LabName)
                .Field("Index", data.Index)
                .Field("IsActive", data.IsActive)
                .Field("IsComplete", data.IsComplete);

            scope.BeginArray("AvailableTools");
            if (data.AvailableTools != null) {
                foreach (var tool in data.AvailableTools) {
                    scope.Item(tool);
                }
            }
            scope.EndArray();

            scope.BeginArray("Prompts");
            if (data.Prompts != null) {
                foreach (var prompt in data.Prompts) {
                    scope.Item(prompt);
                }
            }
            scope.EndArray();
        }

        #endregion // TaskLogData

        #endregion // Types
    }
}
