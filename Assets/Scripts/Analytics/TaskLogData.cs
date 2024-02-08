using System.Collections;
using System.Collections.Generic;
using ThermoVR.Lab;
using ThermoVR.Tools;
using UnityEngine;

namespace ThermoVR.Analytics
{
    public enum TaskCategory
    {
        MULTIPLE_CHOICE,
        WORD_BANK,
        TARGET_STATE,
        MULTIPLE_SELECT,
        CONSTANT_VARIABLE
    }

    public struct IndexedTaskInfo
    {
        public int Index;
        public TaskInfo Info;

        public IndexedTaskInfo(int index, TaskInfo info)
        {
            Index = index;
            Info = info;
        }
    }

    public struct IndexedTopicInfo
    {
        public int Index;
        public TopicInfo Info;

        public IndexedTopicInfo(int index, TopicInfo info)
        {
            Index = index;
            Info = info;
        }
    }

    public struct IndexedLabInfo
    {
        public int Index;
        public LabInfo Info;

        public IndexedLabInfo(int index, LabInfo info)
        {
            Index = index;
            Info = info;
        }
    }


    [System.Serializable]
    public class TaskLogData
    {
        public TaskCategory Category;
        public string LabName; // only set if not in the tasks array of section
        public int SectionIndex;
        public int Index; // task index
        public bool IsActive;
        public bool IsComplete;
        public List<ToolType> AvailableTools;
        public List<string> Prompts;
    }
}