using BeauUtil;
using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using ThermoVR.State;
using ThermoVR.Tools;
using UnityEngine;

namespace ThermoVR.Lab
{
    [Serializable]
    public struct LabInfo {
        public string Name;
        public List<TaskInfo> Tasks;
    }

    [Serializable]
    public struct TaskInfo {
        // Common
        public TaskType TaskType;
        public string TaskQuestion;
        public List<ToolType> AllowedTools;

        // Multiple choice, word bank, multi-select
        public List<string> SecondaryTexts;
        public List<uint> CorrectIDs;

        // Reach State
        public List<SimStateTarget> Targets;
    }

    public enum TaskType {
        MultipleChoice, // single-select multiple choice
        WordBank,
        ReachState, // reach a given state in the sim
        MultipleChoiceMulti // multi-select multiple choice
    }

    public class LabLoader : MonoBehaviour
    {
        #region consts

        private static uint TYPE_INDEX = 1;
        private static uint TOOL_INDEX = 2;
        private static uint SET_INDEX = 3;
        private static uint LIMIT_INDEX = 4;
        private static uint TARGET_INDEX = 5;
        private static uint QUIZ_INDEX = 6;
        private static uint EFFICIENCY_INDEX = 7;

        private static uint NUM_TASK_SECTIONS = 8; // 7 + 1 leading delim

        private static string GROUP_DELIM = "||";
        private static string TASK_INFO_DELIM = "|";
        private static string QUIZ_ANSWER_DELIM = ",";
        private static string STATE_REQ_GROUP_DELIM = ",";
        private static string STATE_REQ_CHUNK_DELIM = ":";
        private static string CORRECT_MARKER = "(c)";


        private static string MC_KEY = "multiple-choice";
        private static string WORD_BANK_KEY = "word-bank";
        private static string REACH_STATE_KEY = "reach-state";
        private static string MC_MULTI_KEY = "multi-select-choice";

        #endregion // Consts

        #region Inspector

        [SerializeField] private TextAsset[] m_initialLabs;
        [SerializeField] private bool m_verboseDebug = false;

        #endregion // Inspector

        #region Unity Callbacks

        private void Start() {
            if (GameMgr.I.IsAlphaRelease) {
                // enable lab mode
                for (int i = 0; i < m_initialLabs.Length; i++) {
                    LoadLab(m_initialLabs[i]);
                }
            }
        }

        #endregion // Unity Callbacks

        #region Editor

        [ContextMenu("Convert Labs To JSON")]
        /// <summary>
        /// Takes Text Assets and converts them to JSON
        /// </summary>
        private void ConvertLabsToJSON() {
            for (int i = 0; i < m_initialLabs.Length; i++) {
                if (TryParseLab(m_initialLabs[i], out LabInfo info)) {
                    string labJSON = JsonUtility.ToJson(info);
                    TextIO.WriteString("Assets/Content/Labs/" + info.Name + ".json", labJSON);
                    Debug.Log("[LabLoader] Converted lab " + info.Name);
                }
                else {
                    Debug.Log("[LabLoad] Failed to convert asset " + m_initialLabs[i].name);
                }
            }
        }

        #endregion // Editor

        #region Lab Parsing

        private void LoadLab(TextAsset labInfoAsset) {
            if (TryParseLab(labInfoAsset, out LabInfo newLabInfo)) {
                Debug.Log("[LabLoad] Asset " + labInfoAsset.name + " loaded successfully.");

                GameMgr.Events.Dispatch(GameEvents.LabLoaded, newLabInfo);
            }
            else {
                Debug.Log("[LabLoad] Failed to load asset " + labInfoAsset.name);
            }
        }

        /// <summary>
        /// Parses a lab
        /// </summary>
        /// <param name="labInfoAsset"></param>
        private bool TryParseLab(TextAsset labInfoAsset, out LabInfo newLabInfo) {
            Debug.Log("[LabLoad] Attempting to parse asset " + labInfoAsset.name + "...");

            bool succeeded = true;
            newLabInfo = new LabInfo();
            if (!TryParseLabFromJSON(labInfoAsset.text, out newLabInfo)) {
                succeeded = false;
            }

            if (!succeeded) {
                succeeded = true;
                newLabInfo.Tasks = new List<TaskInfo>();

                try {
                    List<string> groups = TextIO.TextAssetToList(labInfoAsset, GROUP_DELIM);

                    for (int i = 0; i < groups.Count; i++) {
                        string currGroup = groups[i];

                        if (currGroup.Contains("LAB-NAME")) {
                            ParseLabName(currGroup, ref newLabInfo);
                        }
                        else if (currGroup.Contains("TASK")) {
                            ParseTaskInfo(currGroup, ref newLabInfo);
                        }
                        else if (currGroup.Contains("END")) {
                            // end of file
                            break;
                        }
                        else {
                            // unknown section
                            continue;
                        }
                    }
                }
                catch (InvalidDataException e) {
                    succeeded = false;
                }
            }

            return succeeded;
        }

        /// <summary>
        /// Parses a Lab from JSON
        /// </summary>
        /// <param name="labInfoJSON"></param>
        private bool TryParseLabFromJSON(string labInfoJSON, out LabInfo info) {
            try {
                info = JsonUtility.FromJson<LabInfo>(labInfoJSON);
                return true;
            }
            catch (Exception) {
                // Invalid JSON
                info = new LabInfo();
                return false;
            }
        }

        private void ParseLabName(string group, ref LabInfo labInfo) {
            string labName = group.Substring(group.IndexOf(":") + 1).Trim();
            Debug.Log("[LabLoad] Lab Name: " + labName);

            labInfo.Name = labName;
        }

        #endregion // Lab Parsing

        #region Task Parsing

        private void ParseTaskInfo(string group, ref LabInfo labInfo) {
            TaskInfo newTaskInfo = new TaskInfo();
            newTaskInfo.SecondaryTexts = new List<string>();
            newTaskInfo.CorrectIDs = new List<uint>();
            newTaskInfo.Targets = new List<SimStateTarget>();

            string[] sections = group.Split(TASK_INFO_DELIM);

            if (sections.Length != NUM_TASK_SECTIONS) {
                Debug.Log("[LabLoad] Task was in invalid format! Expecting " + NUM_TASK_SECTIONS + " fields, found " + sections.Length);
                throw new InvalidDataException();
            }

            // Type
            string typeInfo = sections[TYPE_INDEX].Trim();
            if (typeInfo.Contains(MC_KEY)) {
                newTaskInfo.TaskType = TaskType.MultipleChoice;
            }
            else if (typeInfo.Contains(MC_MULTI_KEY)) {
                newTaskInfo.TaskType = TaskType.MultipleChoiceMulti;
            }
            else if (typeInfo.Contains(WORD_BANK_KEY)) {
                newTaskInfo.TaskType = TaskType.WordBank;
            }
            else if (typeInfo.Contains(REACH_STATE_KEY)) {
                newTaskInfo.TaskType = TaskType.ReachState;
            }
            if (m_verboseDebug) { Debug.Log("[LabLoad] Type Info: " + newTaskInfo.TaskType); }

            // Tools
            string toolInfo = sections[TOOL_INDEX].Trim();
            ParseTaskTools(ref toolInfo, ref newTaskInfo);
            if (m_verboseDebug) { Debug.Log("[LabLoad] Tool Info: " + toolInfo); }

            // Set
            string setInfo = sections[SET_INDEX].Trim();
            ParseTaskSets(ref setInfo, ref newTaskInfo);
            if (m_verboseDebug) { Debug.Log("[LabLoad] Set Info: " + setInfo); }

            // Limit
            string limitInfo = sections[LIMIT_INDEX].Trim();
            ParseTaskLimits(ref limitInfo, ref newTaskInfo);
            if (m_verboseDebug) { Debug.Log("[LabLoad] Limit Info: " + limitInfo); }

            // TODO: move quiz Requirements section to the target into block, and rearrange
            // Target
            string targetInfo = sections[TARGET_INDEX].Trim();
            if (m_verboseDebug) { Debug.Log("[LabLoad] Target Info: " + targetInfo); }

            // Quiz
            string quizInfo = sections[QUIZ_INDEX].Trim();
            // get quiz question
            string iterateQuizInfo = "";
            if (quizInfo.Contains("Question:")) {
                ParseTaskQuestion(ref quizInfo, ref iterateQuizInfo, ref newTaskInfo);
            }
            // get quiz options
            if (quizInfo.Contains("Answers:")) {
                ParseTaskAnswers(ref quizInfo, ref iterateQuizInfo, ref newTaskInfo);
            }
            // get reach state requirements
            else if (quizInfo.Contains("Requirements:")) {
                ParseTaskRequirements(ref quizInfo, ref iterateQuizInfo, ref newTaskInfo);
            }
            if (m_verboseDebug) { Debug.Log("[LabLoad] Quiz Info: " + quizInfo); }

            // Efficiency
            string efficiencyInfo = sections[EFFICIENCY_INDEX].Trim();
            if (m_verboseDebug) { Debug.Log("[LabLoad] Efficiency Info: " + efficiencyInfo); }

            // TODO: labInfo.etc = etcInfo

            labInfo.Tasks.Add(newTaskInfo);
        }

        private void ParseTaskTools(ref string toolInfo, ref TaskInfo newTaskInfo) {
            int startIndex = toolInfo.IndexOf('[') + 1;
            if (startIndex == 0) {
                Debug.Log("[LabLoad] Allowed tools definition is invalid");
                return;
            }
            string iterateToolInfo = toolInfo.Substring(startIndex);
            int endIndex = iterateToolInfo.IndexOf(']');
            if (endIndex == -1) {
                Debug.Log("[LabLoad] Allowed tools definition is invalid");
                return;
            }
            int length = endIndex;
            iterateToolInfo = iterateToolInfo.Substring(0, length);
            string[] tools = iterateToolInfo.Split(',');
            bool allowAll = iterateToolInfo.ToLower().Contains("all");
            List<ToolType> allowedTools = new List<ToolType>();

            if (allowAll) {
                allowedTools.Clear();
                foreach (ToolType tool in Enum.GetValues(typeof(ToolType))) {
                    allowedTools.Add(tool);
                }
            }
            else {
                foreach (var tool in tools) {
                    if (Enum.TryParse(tool, true, out ToolType newTool)) {
                        allowedTools.Add(newTool);
                    }
                    else {
                        Debug.Log("[LabLoad] Unrecognized tool type: " + tool);
                    }
                }
            }

            newTaskInfo.AllowedTools = allowedTools;
            if (m_verboseDebug) { Debug.Log("[LabLoad] Task Tools: " + newTaskInfo.AllowedTools); }
        }

        private void ParseTaskSets(ref string setInfo, ref TaskInfo newTaskInfo) {
            // if (m_verboseDebug) { Debug.Log("[LabLoad] Task Sets: " + newTaskInfo.Tools); }
        }

        private void ParseTaskLimits(ref string limitInfo, ref TaskInfo newTaskInfo) {
            // if (m_verboseDebug) { Debug.Log("[LabLoad] Task Limits: " + newTaskInfo.Tools); }
        }

        private void ParseTaskQuestion(ref string quizInfo, ref string iterateQuizInfo, ref TaskInfo newTaskInfo) {
            int preIndex = quizInfo.IndexOf("Question:");
            iterateQuizInfo = quizInfo.Substring(preIndex);
            int startIndex = iterateQuizInfo.IndexOf('"') + 1;
            iterateQuizInfo = iterateQuizInfo.Substring(startIndex);
            int endIndex = iterateQuizInfo.IndexOf('"');
            int length = endIndex;
            iterateQuizInfo = iterateQuizInfo.Substring(0, length);
            newTaskInfo.TaskQuestion = iterateQuizInfo;
            if (m_verboseDebug) { Debug.Log("[LabLoad] Task Question: " + newTaskInfo.TaskQuestion); }
        }

        private void ParseTaskAnswers(ref string quizInfo, ref string iterateQuizInfo, ref TaskInfo newTaskInfo) {
            int preIndex = quizInfo.IndexOf("Answers:");
            iterateQuizInfo = quizInfo.Substring(preIndex);
            int startIndex = iterateQuizInfo.IndexOf('[') + 1;
            int endIndex = iterateQuizInfo.IndexOf(']');
            int length = endIndex - startIndex;
            iterateQuizInfo = iterateQuizInfo.Substring(startIndex, length);
            string[] rawAnswers = iterateQuizInfo.Split(QUIZ_ANSWER_DELIM);

            for (int i = 0; i < rawAnswers.Length; i++) {
                if (m_verboseDebug) { Debug.Log("[LabLoad] Parsing answer: " + rawAnswers[i]); }

                if (rawAnswers[i].ToLower().Contains(CORRECT_MARKER)) {
                    int correctIndex = i;
                    newTaskInfo.CorrectIDs.Add((uint)correctIndex);
                    if (m_verboseDebug) { Debug.Log("[LabLoad] correct index: " + correctIndex); }
                }

                int startAnswerIndex = rawAnswers[i].IndexOf('"') + 1;
                rawAnswers[i] = rawAnswers[i].Substring(startAnswerIndex);
                int endAnswerIndex = rawAnswers[i].IndexOf('"');
                length = endAnswerIndex;
                string parsedAnswer = rawAnswers[i].Substring(0, length);
                newTaskInfo.SecondaryTexts.Add(parsedAnswer);
            }
            if (m_verboseDebug) { Debug.Log("[LabLoad] Question Answers: " + newTaskInfo.SecondaryTexts); }
        }

        private void ParseTaskRequirements(ref string quizInfo, ref string iterateQuizInfo, ref TaskInfo newTaskInfo) {
            int preIndex = quizInfo.IndexOf("Requirements:");
            iterateQuizInfo = quizInfo.Substring(preIndex);
            int startIndex = iterateQuizInfo.IndexOf('[') + 1;
            int endIndex = iterateQuizInfo.IndexOf(']');
            int length = endIndex - startIndex;
            iterateQuizInfo = iterateQuizInfo.Substring(startIndex, length);
            string[] rawReqs = iterateQuizInfo.Split(STATE_REQ_GROUP_DELIM);

            for (int i = 0; i < rawReqs.Length; i++) {
                if (m_verboseDebug) { Debug.Log("[LabLoad] Parsing requirements: " + rawReqs[i]); }

                int startAnswerIndex = rawReqs[i].IndexOf('(') + 1;
                rawReqs[i] = rawReqs[i].Substring(startAnswerIndex);
                int endAnswerIndex = rawReqs[i].IndexOf(')');
                length = endAnswerIndex;
                string reqGroup = rawReqs[i].Substring(0, length);

                string[] reqChunks = reqGroup.Split(STATE_REQ_CHUNK_DELIM);
                VarID reqID = (VarID)System.Enum.Parse(typeof(VarID), reqChunks[0]);
                float reqVal = float.Parse(reqChunks[1]);
                float reqRange = float.Parse(reqChunks[2]);

                SimStateTarget newTarget = new SimStateTarget(reqID, reqVal, reqRange);
                newTaskInfo.Targets.Add(newTarget);
            }
            if (m_verboseDebug) { Debug.Log("[LabLoad] Question Answers: " + newTaskInfo.SecondaryTexts); }
        }

        #endregion // Task Parsing
    }
}