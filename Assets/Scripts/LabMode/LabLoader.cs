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
    public struct LabInfo
    {
        public string Name;
        public string Author;
        public List<TopicInfo> Topics;
    }

    /// <summary>
    /// Defines a limit for a given 
    /// </summary>
    [Serializable]
    public struct LimitPair
    {
        public double Limit;
    }

    [Serializable]
    public struct SetGroup
    {
        public double P;
        public double V;
        public double T;

        public SetGroup(double p, double v, double t) {
            P = p;
            V = v;
            T = t;
        }

        public bool IsEmpty()
        {
            return P == -1 && V == -1 && T == -1;
        }
    }

    [Serializable]
    public struct LimitsGroup
    {
        // One for each limitable sim variable
        public LimitBounds Pressure;
        public LimitBounds Temperature;
        public LimitBounds Volume;
        public LimitBounds InternalEnergy;
        public LimitBounds Entropy;
        public LimitBounds Enthalpy;
        public LimitBounds Quality;

        public LimitsGroup(bool emptyInit) {
            Pressure = new LimitBounds(emptyInit);
            Temperature = new LimitBounds(emptyInit);
            Volume = new LimitBounds(emptyInit);
            InternalEnergy = new LimitBounds(emptyInit);
            Entropy = new LimitBounds(emptyInit);
            Enthalpy = new LimitBounds(emptyInit);
            Quality = new LimitBounds(emptyInit);
        }
    }

    [Serializable]
    public struct LimitBounds
    {
        public double Ceiling;
        public double Floor;

        public LimitBounds(bool emptyInit) {
            Ceiling = -1;
            Floor = -1;
        }
    }

    [Serializable]
    public struct TopicInfo
    {
        public string TopicHeader;
        public List<TaskInfo> Tasks;
    }

    [Serializable]
    public struct TaskInfo
    {
        // Common
        public TaskType TaskType;
        public string InitialConditions;
        public List<string> TaskQuestions;
        public List<ToolType> AllowedTools;
        public bool GrabAllowed;
        public SetGroup Sets; // p, v, and t values to set
        public LimitsGroup Limits; // stores limits for simulation variables

        // Multiple choice, word bank, multi-select
        public List<string> SecondaryTexts;
        public List<uint> CorrectIDs;

        // Reach State
        public List<SimStateTarget> Targets;
    }

    public enum TaskType
    {
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

        private static uint TOPIC_HEADER_INDEX = 1;

        private static uint NUM_TOPIC_SECTIONS = 2; // 1 + 1 leading delim

        private static string GROUP_DELIM = "||";
        private static string TASK_INFO_DELIM = "|";
        private static string TOPIC_INFO_DELIM = "|";
        private static string QUIZ_ANSWER_DELIM = ",";
        private static string STATE_REQ_GROUP_DELIM = ",";
        private static string STATE_REQ_CHUNK_DELIM = ":";
        private static string CORRECT_MARKER = "(c)";
        private static string TOOL_DELIM = ",";
        private static string SET_GROUP_DELIM = ",";
        private static string SET_CHUNK_DELIM = ":";
        private static string LIMIT_GROUP_DELIM = ",";
        private static string LIMIT_CHUNK_DELIM = ":";


        private static string MC_KEY = "multiple-choice";
        private static string WORD_BANK_KEY = "word-bank";
        private static string REACH_STATE_KEY = "reach-state";
        private static string MC_MULTI_KEY = "multi-select-choice";

        private static string QUESTIONS_START_KEY = "--Questions:";
        private static string QUESTIONS_END_KEY = "Questions-->";

        private static string QUESTION_START_KEY = "Question:";

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
                newLabInfo.Topics = new List<TopicInfo>();

                try
                {
                    List<string> groups = TextIO.TextAssetToList(labInfoAsset, GROUP_DELIM);

                    for (int i = 0; i < groups.Count; i++) {
                        string currGroup = groups[i];

                        if (currGroup.Contains("LAB-NAME")) {
                            ParseLabName(currGroup, ref newLabInfo);
                        }
                        else if (currGroup.Contains("LAB-AUTHOR"))
                        {
                            ParseLabAuthor(currGroup, ref newLabInfo);
                        }
                        else if (currGroup.Contains("TOPIC"))
                        {
                            ParseTopicInfo(currGroup, ref newLabInfo);
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
            string labName = group.Substring(group.IndexOf("LAB-NAME:") + "LAB-NAME:".Length).Trim();
            Debug.Log("[LabLoad] Lab Name: " + labName);

            labInfo.Name = labName;
        }

        private void ParseLabAuthor(string group, ref LabInfo labInfo)
        {
            string labAuthor = group.Substring(group.IndexOf("LAB-AUTHOR:") + "LAB-AUTHOR:".Length).Trim();
            Debug.Log("[LabLoad] Lab Author: " + labAuthor);

            labInfo.Author = labAuthor;
        }

        #endregion // Lab Parsing

        #region Task Parsing

        private void ParseTopicInfo(string group, ref LabInfo labInfo) {
            TopicInfo newTopicInfo = new TopicInfo();
            newTopicInfo.Tasks = new List<TaskInfo>();

            string[] sections = group.Split(TOPIC_INFO_DELIM);

            if (sections.Length != NUM_TOPIC_SECTIONS)
            {
                Debug.Log("[LabLoad] TOPIC was in invalid format! Expecting " + NUM_TOPIC_SECTIONS + " fields, found " + sections.Length);
                throw new InvalidDataException();
            }

            // HEADER
            string headerInfo = sections[TOPIC_HEADER_INDEX].Trim();
            string iterateTopicInfo = "";
            if (headerInfo.Contains("Header:")) { 
                ParseTopicHeader(ref headerInfo, ref iterateTopicInfo, ref newTopicInfo);
            }
            if (m_verboseDebug) { Debug.Log("[LabLoad] Topic Header Info: " + newTopicInfo); }

            labInfo.Topics.Add(newTopicInfo);
        }

        private void ParseTopicHeader(ref string headerInfo, ref string iterateTopicInfo, ref TopicInfo newTopicInfo)
        {
            int preIndex = headerInfo.IndexOf("Header:");
            iterateTopicInfo = headerInfo.Substring(preIndex);
            int startIndex = iterateTopicInfo.IndexOf('"') + 1;
            iterateTopicInfo = iterateTopicInfo.Substring(startIndex);
            int endIndex = iterateTopicInfo.IndexOf('"');
            int length = endIndex;
            iterateTopicInfo = iterateTopicInfo.Substring(0, length);
            newTopicInfo.TopicHeader = iterateTopicInfo;
            if (m_verboseDebug) { Debug.Log("[LabLoad] Topic Header: " + newTopicInfo.TopicHeader); }
        }

        private void ParseTaskInfo(string group, ref LabInfo labInfo) {
            TaskInfo newTaskInfo = new TaskInfo();
            newTaskInfo.TaskQuestions = new List<string>();
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
            if (quizInfo.Contains("Initial Conditions:"))
            {
                ParseTaskInitConditions(ref quizInfo, ref iterateQuizInfo, ref newTaskInfo);
            }
            if (quizInfo.Contains(QUESTIONS_START_KEY)) {
                ParseTaskQuestions(ref quizInfo, ref iterateQuizInfo, ref newTaskInfo);
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

            labInfo.Topics[labInfo.Topics.Count - 1].Tasks.Add(newTaskInfo);
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
            string[] tools = iterateToolInfo.Split(TOOL_DELIM);
            bool allowAll = iterateToolInfo.ToLower().Contains("all");
            List<ToolType> allowedTools = new List<ToolType>();

            newTaskInfo.GrabAllowed = false;

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

            foreach (var tool in tools)
            {
                if (tool.Trim().Equals("Grab"))
                {
                    // enable grab
                    newTaskInfo.GrabAllowed = true;
                }
            }

            newTaskInfo.AllowedTools = allowedTools;
            if (m_verboseDebug) { Debug.Log("[LabLoad] Task Tools: " + newTaskInfo.AllowedTools); }
        }

        private void ParseTaskSets(ref string setInfo, ref TaskInfo newTaskInfo) {
            // Example format is [(p:5), (v:5), (t:5)]
            SetGroup newSets = new SetGroup(-1, -1, -1);
            newTaskInfo.Sets = newSets;

            int startIndex = setInfo.IndexOf('[') + 1;
            if (startIndex == 0) {
                Debug.Log("[LabLoad] SetInfo definition is invalid");
                return;
            }
            string iterateSetInfo = setInfo.Substring(startIndex);
            int endIndex = iterateSetInfo.IndexOf(']');
            if (endIndex == -1) {
                Debug.Log("[LabLoad] SetInfo definition is invalid");
                return;
            }
            int length = endIndex;
            iterateSetInfo = iterateSetInfo.Substring(0, length);
            string[] setGroups = iterateSetInfo.Split(SET_GROUP_DELIM);
            double p, v, t;
            p = v = t = -1;

            foreach (string group in setGroups) {
                if (group.Contains('p')) {
                    if (TryExtractSetVal(group, out double val)) {
                        p = val * 1000; // defined in kPa
                    }
                    else {
                        Debug.Log("[LabLoad] SetInfo for p is invalid");
                        continue;
                    }
                }
                else if (group.Contains("v")) {
                    if (TryExtractSetVal(group, out double val)) {
                        v = val;
                    }
                    else {
                        Debug.Log("[LabLoad] SetInfo for v is invalid");
                        continue;
                    }
                }
                else if (group.Contains("t")) {
                    if (TryExtractSetVal(group, out double val)) {
                        t = val;
                    }
                    else {
                        Debug.Log("[LabLoad] SetInfo for t is invalid");
                        continue;
                    }
                }
            }

            newSets.P = p;
            newSets.V = v;
            newSets.T = t;
            newTaskInfo.Sets = newSets;
            if (m_verboseDebug) { Debug.Log("[LabLoad] Task Sets: " + newTaskInfo.Sets); }
        }

        private bool TryExtractSetVal(string group, out double val) {
            int startValIndex, endValIndex;
            startValIndex = group.IndexOf(SET_CHUNK_DELIM) + 1;
            endValIndex = group.IndexOf(')');
            if (startValIndex == 0 || endValIndex == -1) {
                val = -1;
                return false;
            }

            string valStr = group.Substring(startValIndex, endValIndex - startValIndex);
            double.TryParse(valStr, out val);
            return true;
        }

        private void ParseTaskLimits(ref string limitInfo, ref TaskInfo newTaskInfo) {
            // Example format is [(pressure:10:ceiling), (pressure:1:floor), (enthalpy:10:floor)]
            LimitsGroup newLimits = new LimitsGroup(true);
            newTaskInfo.Limits = newLimits;

            int startIndex = limitInfo.IndexOf('[') + 1;
            if (startIndex == 0) {
                Debug.Log("[LabLoad] LimitInfo definition is invalid");
                return;
            }
            string iterateLimitInfo = limitInfo.Substring(startIndex);
            int endIndex = iterateLimitInfo.IndexOf(']');
            if (endIndex == -1) {
                Debug.Log("[LabLoad] LimitInfo definition is invalid");
                return;
            }
            int length = endIndex;
            iterateLimitInfo = iterateLimitInfo.Substring(0, length);
            string[] limitGroups = iterateLimitInfo.Split(LIMIT_GROUP_DELIM);

            for (int i = 0; i < limitGroups.Length; i++) {
                string group = limitGroups[i];

                // remove parentheses from front and end
                group = group.Trim();
                group = group.Substring(1, group.Length - 2);

                // parse the three components: varID, value, and floor/ceiling
                string[] limitChunks = group.Split(LIMIT_CHUNK_DELIM);
                VarID simVar;
                double val = -1;
                int boundaryType = 0; // floor = -1, ceiling = 1, default = 0
                bool validInputs = true;

                validInputs = Enum.TryParse(limitChunks[0], true, out simVar) ? validInputs : false;
                validInputs = double.TryParse(limitChunks[1], out val) ? validInputs : false;
                if (limitChunks[2].Contains('c')) {
                    boundaryType = 1;
                }
                else if (limitChunks[2].Contains('f')) {
                    boundaryType = -1;
                }
                validInputs = boundaryType != 0 ? validInputs : false;


                if (validInputs) {
                    RecordLimit(simVar, val, boundaryType, ref newLimits);
                }
            }

            newTaskInfo.Limits = newLimits;
            if (m_verboseDebug) { Debug.Log("[LabLoad] Task Limits: " + newTaskInfo.Limits); }
        }

        private void RecordLimit(VarID simVar, double val, int boundaryType, ref LimitsGroup limitGroup) {
            switch (simVar) {
                case VarID.Pressure:
                    if (boundaryType == 1) { limitGroup.Pressure.Ceiling = val; }
                    // else boundary type = -1
                    else { limitGroup.Pressure.Floor = val; }
                    break;
                case VarID.Temperature:
                    if (boundaryType == 1) { limitGroup.Temperature.Ceiling = val; }
                    // else boundary type = -1
                    else { limitGroup.Temperature.Floor = val; }
                    break;
                case VarID.Volume:
                    if (boundaryType == 1) { limitGroup.Volume.Ceiling = val; }
                    // else boundary type = -1
                    else { limitGroup.Volume.Floor = val; }
                    break;
                case VarID.InternalEnergy:
                    if (boundaryType == 1) { limitGroup.InternalEnergy.Ceiling = val; }
                    // else boundary type = -1
                    else { limitGroup.InternalEnergy.Floor = val; }
                    break;
                case VarID.Entropy:
                    if (boundaryType == 1) { limitGroup.Entropy.Ceiling = val; }
                    // else boundary type = -1
                    else { limitGroup.Entropy.Floor = val; }
                    break;
                case VarID.Enthalpy:
                    if (boundaryType == 1) { limitGroup.Enthalpy.Ceiling = val; }
                    // else boundary type = -1
                    else { limitGroup.Enthalpy.Floor = val; }
                    break;
                case VarID.Quality:
                    if (boundaryType == 1) { limitGroup.Quality.Ceiling = val; }
                    // else boundary type = -1
                    else { limitGroup.Quality.Floor = val; }
                    break;
                default:
                    break;
            }
        }

        private void ParseTaskInitConditions(ref string quizInfo, ref string iterateQuizInfo, ref TaskInfo newTaskInfo)
        {
            int preIndex = quizInfo.IndexOf("Initial Conditions:");
            iterateQuizInfo = quizInfo.Substring(preIndex);
            int startIndex = iterateQuizInfo.IndexOf('"') + 1;
            iterateQuizInfo = iterateQuizInfo.Substring(startIndex);
            int endIndex = iterateQuizInfo.IndexOf('"');
            int length = endIndex;
            iterateQuizInfo = iterateQuizInfo.Substring(0, length);
            newTaskInfo.InitialConditions = iterateQuizInfo;

            if (m_verboseDebug) { Debug.Log("[LabLoad] Initial Conditions: " + newTaskInfo.InitialConditions); }
        }

        private void ParseTaskQuestions(ref string quizInfo, ref string iterateQuizInfo, ref TaskInfo newTaskInfo) {
            // isolate questions chunk
            int questionsChunkPreIndex = quizInfo.IndexOf(QUESTIONS_START_KEY);
            int questionsChunkEndIndex = quizInfo.IndexOf(QUESTIONS_END_KEY);

            if (questionsChunkEndIndex == -1 || questionsChunkPreIndex == -1)
            {
                return;
            }

            iterateQuizInfo = quizInfo.Substring(questionsChunkPreIndex, questionsChunkEndIndex - questionsChunkPreIndex);

            int preIndex;
            int startIndex;
            int endIndex;
            int length;

            string qStr;

            int iters = 0; // prevent infinite loop if incorrect format

            while (iterateQuizInfo.Contains(QUESTION_START_KEY) || iters > 10)
            {
                preIndex = iterateQuizInfo.IndexOf(QUESTION_START_KEY);
                qStr = iterateQuizInfo.Substring(preIndex);
                startIndex = qStr.IndexOf('"') + 1;
                qStr = qStr.Substring(startIndex);
                endIndex = qStr.IndexOf('"');
                length = endIndex;
                qStr = qStr.Substring(0, length);
                newTaskInfo.TaskQuestions.Add(qStr);
                if (m_verboseDebug) { Debug.Log("[LabLoad] Task Question: " + qStr); }

                iterateQuizInfo = iterateQuizInfo.Substring(preIndex + startIndex + length);

                iters++;
            }
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