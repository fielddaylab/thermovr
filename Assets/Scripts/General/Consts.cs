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
    static public readonly StringHash32 ColliderReleased = "world:cartridge-released"; // Collider
    static public readonly StringHash32 ActivateCartridge = "world:activate-cartridge"; // Cartridge
    static public readonly StringHash32 ActivateTool = "sim:activate-tool"; // Tool
    static public readonly StringHash32 DetachTool = "sim:detach-tool"; // Tool
    static public readonly StringHash32 StoreTool = "sim:store-tool"; // Tool
    static public readonly StringHash32 UpdateToolText = "tool:update-text"; // Tool
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
    static public readonly string AmbientPressure = "psi";
    static public readonly string Volume = "m³/kg";
    static public readonly string InternalEnergy = "kJ/kg";
    static public readonly string Entropy = "kJ/kgK";
    static public readonly string Enthalpy = "kJ/kg";
    static public readonly string Quality = "%";
}