using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using BeauUtil;

static public class GameEvents
{
    static public readonly StringHash32 GatherPressables = "world:gather-pressables";
    static public readonly StringHash32 RegisterPressable = "world:register-pressable"; // Pressable
    static public readonly StringHash32 RegisterMovable = "world:register-movable"; // Touchable
    static public readonly StringHash32 CheckForPress = "world:check-for-press"; // bool
    static public readonly StringHash32 LabLoaded = "world:lab-loaded"; // LabInfo
    static public readonly StringHash32 ColliderReleased = "world:collider-released"; // Collider
    static public readonly StringHash32 ColliderGrabbed = "world:collider-grabbed"; // Collider
    static public readonly StringHash32 ActivateCartridge = "world:activate-cartridge"; // Cartridge
    static public readonly StringHash32 DeactivateCartridge = "world:deactivate-cartridge"; // Cartridge
    static public readonly StringHash32 TaskResetPressed = "world:task-reset-pressed";
    static public readonly StringHash32 BeginLab = "world:begin-lab";
    static public readonly StringHash32 PressedToolToggle = "world:pressed-tool-toggle";
    static public readonly StringHash32 WarpPVT = "sim:warp-pvt";
    static public readonly StringHash32 ActivateTool = "sim:activate-tool"; // Tool
    static public readonly StringHash32 DeactivateTool = "sim:deactivate-tool"; // Tool
    static public readonly StringHash32 UpdateAllowedTools = "sim:update-allowed-tools"; // List<ToolType>
    static public readonly StringHash32 ResetPressed = "sim:reset-pressed";
    static public readonly StringHash32 ResetToolRestrictions = "sim:reset-tool-restrictions";
    static public readonly StringHash32 UpdateVaporFlow = "sim:vapor-update-flow"; // double
    static public readonly StringHash32 UpdateVarText = "sim:update-var-text"; // VarUpdate
    static public readonly StringHash32 UpdateGraphSetting = "sim:update-graph-setting"; // GraphSettingUpdate

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