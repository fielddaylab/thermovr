using System.Collections;
using System.Collections.Generic;
using ThermoVR.Tools;
using UnityEngine;

namespace ThermoVR.Analytics
{
    public enum TaskCategory
    {
        MULTIPLE_CHOICE,
        MULTIPLE_SELECT,
        WORD_BANK,
        REACH_STATE,
        CONSTANT_VARIABLE
    }

    [System.Serializable]
    public class TaskLogData
    {
        public TaskCategory Category;
        public string LabName; // only set if not in the tasks array of section
        public int SectionIndex;
        public bool IsActive;
        public List<ToolType> AvailableTools;
        public List<string> Prompts;
    }
}