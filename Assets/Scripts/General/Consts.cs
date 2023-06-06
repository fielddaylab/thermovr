using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using BeauUtil;

static public class GameEvents
{
    static public readonly StringHash32 ActivateTool = "sim:activate-tool"; // Tool
    static public readonly StringHash32 DetachTool = "sim:detach-tool"; // Tool
    static public readonly StringHash32 StoreTool = "sim:store-tool"; // Tool
    static public readonly StringHash32 UpdateToolText = "tool:update-text"; // Tool
    static public readonly StringHash32 UpdateVaporFlow = "sim:vapor-update-flow"; // double delta_pressure
}
