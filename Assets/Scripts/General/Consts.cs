using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using BeauUtil;

static public class GameEvents
{
    static public readonly StringHash32 InitialLoadComplete = "world:initial-load-compelete";
    static public readonly StringHash32 GatherPressables = "world:gather-pressables";
    static public readonly StringHash32 RegisterPressable = "world:register-pressable"; // Pressable
    static public readonly StringHash32 RegisterMovable = "world:register-movable"; // Touchable
    static public readonly StringHash32 CheckForPress = "world:check-for-press"; // bool
    static public readonly StringHash32 LabLoaded = "world:lab-loaded"; // LabInfo
    static public readonly StringHash32 ColliderReleased = "world:collider-released"; // Collider
    static public readonly StringHash32 ColliderGrabbed = "world:collider-grabbed"; // Collider
    static public readonly StringHash32 ObjectReleased = "world:object-released"; // GameObject
    static public readonly StringHash32 ObjectGrabbed = "world:object-grabbed"; // GameObject
    static public readonly StringHash32 PreActivateLab = "world:preactivate-lab"; // LabInfo
    static public readonly StringHash32 ActivateLab = "world:activate-lab"; // LabInfo
    static public readonly StringHash32 DeactivateLab = "world:deactivate-lab"; // LabInfo
    static public readonly StringHash32 TaskResetPressed = "world:task-reset-pressed";
    static public readonly StringHash32 BeginLab = "world:begin-lab";
    static public readonly StringHash32 PressedToolToggle = "world:pressed-tool-toggle";
    static public readonly StringHash32 WarpPVT = "sim:warp-pvt";
    static public readonly StringHash32 ActivateTool = "sim:activate-tool"; // Tool
    static public readonly StringHash32 DeactivateTool = "sim:deactivate-tool"; // Tool
    static public readonly StringHash32 AllowTool = "sim:allow-tool"; // Tool
    static public readonly StringHash32 DisallowTool = "sim:disallow-tool"; // Tool
    static public readonly StringHash32 UpdateAllowedTools = "sim:update-allowed-tools"; // List<ToolType>
    static public readonly StringHash32 ResetPressed = "sim:reset-pressed";
    static public readonly StringHash32 ResetSimClicked = "sim:reset-sim-clicked";
    static public readonly StringHash32 ResetToolRestrictions = "sim:reset-tool-restrictions";
    static public readonly StringHash32 UpdateVaporFlow = "sim:vapor-update-flow"; // double
    static public readonly StringHash32 UpdateVarText = "sim:update-var-text"; // VarUpdate
    static public readonly StringHash32 UpdateGraphSetting = "sim:update-graph-setting"; // GraphSettingUpdate
    static public readonly StringHash32 UISwitched = "sim:ui-switched"; // GraphSettingUpdate

    static public readonly StringHash32 TryNewName = "game:new-name";
    static public readonly StringHash32 NewNameGenerated = "game:new-name-generated";
    static public readonly StringHash32 StartSession = "game:start-session";
    static public readonly StringHash32 StartGame = "game:start-game";

    static public readonly StringHash32 HandStartPress = "world:hand-start-press";
    
    static public readonly StringHash32 SelectLab = "lab:select-lab";
    static public readonly StringHash32 ClickLabHome = "lab:click-lab-home";
    static public readonly StringHash32 ClickSelectTask = "lab:click-select-task";
    static public readonly StringHash32 ClickSelectSection = "lab:click-select-section";
    static public readonly StringHash32 SectionSwitched = "lab:section-switched"; // int
    static public readonly StringHash32 TaskSwitched = "lab:task-switched"; // int
    static public readonly StringHash32 ClickTaskScrollLeft = "lab:task-scroll-left";
    static public readonly StringHash32 ClickTaskScrollRight = "lab:task-scroll-right";
    static public readonly StringHash32 ClickSectionScrollUp = "lab:section-scroll-up";
    static public readonly StringHash32 ClickSectionScrollDown = "lab:section-scroll-down";
    static public readonly StringHash32 ClickLabScrollUp = "lab:lab-scroll-up";
    static public readonly StringHash32 ClickLabScrollDown = "lab:lab-scroll-down";
    static public readonly StringHash32 TaskListDisplayed = "lab:task-list-displayed"; // List<IndexedTaskInfo>
    static public readonly StringHash32 SectionListDisplayed = "lab:section-list-displayed"; // List<IndexedTopicInfo>
    static public readonly StringHash32 LabMenuDisplayed = "lab:lab-menu-displayed"; // List<IndexedLabInfo>
    static public readonly StringHash32 TaskChoiceSelected = "lab:task-choice-selected";
    static public readonly StringHash32 TargetStateReached = "lab:target-state-reached";
    static public readonly StringHash32 TargetStateLost = "lab:target-state-lost";
    static public readonly StringHash32 ClickSelectAnswer = "lab:click-select-answer";
    static public readonly StringHash32 ClickDeselectAnswer = "lab:click-deselect-answer";
    static public readonly StringHash32 ClickSubmitAnswer = "lab:click-submit-answer";
    static public readonly StringHash32 ClickResetQuiz = "lab:click-reset-quiz";
    static public readonly StringHash32 ClickOpenWordBank = "lab:click-open-word-bank";
    static public readonly StringHash32 WordBankDisplayed = "lab:word-bank-displayed";
    static public readonly StringHash32 WordBankClosed = "lab:word-bank-closed";
    static public readonly StringHash32 SandboxModeClicked = "lab:sandbox-mode-clicked";
    static public readonly StringHash32 LabModeClicked = "lab:lab-mode-clicked";
    static public readonly StringHash32 SettingsViewClicked = "lab:settings-view-clicked";
    static public readonly StringHash32 TaskCompleted = "lab:task-completed";
    static public readonly StringHash32 SectionCompleted = "lab:section-completed";
    static public readonly StringHash32 LabCompleted = "lab:lab-completed";

    static public readonly StringHash32 TabletGrabbed = "world:tablet-grabbed";
    static public readonly StringHash32 TabletReleased = "world:tablet-released";
    static public readonly StringHash32 WorkspaceHandleGrabbed = "world:workspace-handle-grabbed";
    static public readonly StringHash32 WorkspaceHandleReleased = "world:workspace-handle-released";
    static public readonly StringHash32 RotateGraphClickedCW = "world:rotate-graph-clicked-cw";
    static public readonly StringHash32 RotateGraphClickedCCW = "world:rotate-graph-clicked-ccw";
    static public readonly StringHash32 GraphBallGrabbed = "world:graph-ball-grabbed";
    static public readonly StringHash32 GraphBallReleased = "world:graph-ball-released";
    static public readonly StringHash32 HeadsetOn = "world:headset-on";
    static public readonly StringHash32 HeadsetOff = "world:headset-off";
    static public readonly StringHash32 GazeEnd = "world:gaze-end";

    static public readonly StringHash32 ToolTogglePressed = "sim:tool-toggle-pressed";
    static public readonly StringHash32 ClickToolIncrease = "sim:click-tool-increase";
    static public readonly StringHash32 ClickToolDecrease = "sim:click-tool-decrease";
    static public readonly StringHash32 GrabToolSlider = "sim:grab-tool-slider";
    static public readonly StringHash32 ReleaseToolSlider = "sim:release-tool-slider";
}

static public class ObjectIDs
{
    static public readonly string CenterEyeAnchor = "CenterEyeAnchor";
}

static public class Units
{
    static public readonly string Weight = "kg";
    static public readonly string Heat = "kJ/s";
    static public readonly string TemperatureK = "K";
    static public readonly string TemperatureC = "°C";
    static public readonly string Percent = "%";

    static public readonly string Pressure = "kPa";
    static public readonly string AmbientPressure = "kPa";
    static public readonly string Volume = "m³/kg";
    static public readonly string InternalEnergy = "kJ/kg";
    static public readonly string Entropy = "kJ/kgK";
    static public readonly string Enthalpy = "kJ/kg";
    static public readonly string Quality = "%";
}

// Used by graph elements, dial readouts, and sensor readouts
static public class DigitFormat
{
    static public readonly string Weight = "{0:#.00e+0}";
    static public readonly string Heat = "{0:0.00}";
    static public readonly string TemperatureK = "{0:#.00e+0}";
    static public readonly string Percent = "{0:0.00}";

    static public readonly string Pressure = "{0:#.00e+0}";
    static public readonly string AmbientPressure = "{0:#.00e+0}";
    static public readonly string Volume = "{0:#.00e+0}";
    static public readonly string InternalEnergy = "{0:#.##e+0}";
    static public readonly string Entropy = "{0:#.##e+0}";
    static public readonly string Enthalpy = "{0:#.##e+0}";
    static public readonly string Quality = "{0:0.00}";
}