using BeauUtil;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using UnityEngine;

namespace ThermoVR.Lab
{
    public struct LabInfo {
        public string Name;
        public List<TaskInfo> Tasks;
    }

    public struct TaskInfo {
        public TaskType TaskType;
        public string TaskQuestion;
        public List<string> SecondaryTexts;
        public uint CorrectID;
    }

    public enum TaskType {
        MultipleChoice,
        WordBank,
        ReachState // reach a given state in the sim
    }

    public class LabLoader : MonoBehaviour
    {
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
        private static string CORRECT_MARKER = "(c)";


        private static string MC_KEY = "multiple-choice";
        private static string WORD_BANK_KEY = "word-bank";
        private static string REACH_STATE_KEY = "reach-state";

        [SerializeField] private TextAsset[] m_initialLabs;
        [SerializeField] private bool m_verboseDebug = false;

        private void Start() {
            for (int i = 0; i < m_initialLabs.Length; i++) {
                LoadLab(m_initialLabs[i]);
            }
        }

        private void LoadLab(TextAsset labInfoAsset) {
            Debug.Log("[LabLoad] Attempting to load asset " + labInfoAsset.name + "...");

            bool succeeded = true;
            LabInfo newLabInfo = new LabInfo();
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
            catch (InvalidDataException e){
                succeeded = false;
            }

            // Terminate

            if (succeeded) {
                Debug.Log("[LabLoad] Asset " + labInfoAsset.name + " loaded successfully.");

                GameMgr.Events.Dispatch(GameEvents.LabLoaded, newLabInfo);
            }
            else {
                Debug.Log("[LabLoad] Failed to load asset " + labInfoAsset.name);
            }
        }

        private void ParseLabName(string group, ref LabInfo labInfo) {
            string labName = group.Substring(group.IndexOf(":") + 1).Trim();
            Debug.Log("[LabLoad] Lab Name: " + labName);

            labInfo.Name = labName;
        }

        private void ParseTaskInfo(string group, ref LabInfo labInfo) {
            TaskInfo newTaskInfo = new TaskInfo();
            newTaskInfo.SecondaryTexts = new List<string>();

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
            else if (typeInfo.Contains(WORD_BANK_KEY)) {
                newTaskInfo.TaskType = TaskType.WordBank;
            }
            else if (typeInfo.Contains(REACH_STATE_KEY)) {
                newTaskInfo.TaskType = TaskType.ReachState;
            }
            if (m_verboseDebug) { Debug.Log("[LabLoad] Type Info: " + newTaskInfo.TaskType); }

            // Tools
            string toolInfo = sections[TOOL_INDEX].Trim();
            if (m_verboseDebug) { Debug.Log("[LabLoad] Tool Info: " + toolInfo); }

            // Set
            string setInfo = sections[SET_INDEX].Trim();
            if (m_verboseDebug) { Debug.Log("[LabLoad] Set Info: " + setInfo); }

            // Limit
            string limitInfo = sections[LIMIT_INDEX].Trim();
            if (m_verboseDebug) { Debug.Log("[LabLoad] Limit Info: " + limitInfo); }

            // Target
            string targetInfo = sections[TARGET_INDEX].Trim();
            if (m_verboseDebug) { Debug.Log("[LabLoad] Target Info: " + targetInfo); }

            // Quiz
            string quizInfo = sections[QUIZ_INDEX].Trim();
            // get quiz question
            string iterateQuizInfo;
            if (quizInfo.Contains("Question:")) {
                int preIndex = quizInfo.IndexOf("Question:");
                iterateQuizInfo = quizInfo.Substring(preIndex);
                int startIndex = iterateQuizInfo.IndexOf('"');
                iterateQuizInfo = iterateQuizInfo.Substring(startIndex + 1);
                int endIndex = iterateQuizInfo.IndexOf('"');
                iterateQuizInfo = iterateQuizInfo.Substring(0, endIndex);
                newTaskInfo.TaskQuestion = iterateQuizInfo;
                if (m_verboseDebug) { Debug.Log("[LabLoad] Task Question: " + newTaskInfo.TaskQuestion); }
            }
            // get quiz options
            if (quizInfo.Contains("Answers:")) {
                int preIndex = quizInfo.IndexOf("Answers:");
                iterateQuizInfo = quizInfo.Substring(preIndex);
                int startIndex = iterateQuizInfo.IndexOf('[');
                int endIndex = iterateQuizInfo.IndexOf(']');
                Debug.Log("[LabLoad] pre: " + preIndex + " || start: " + startIndex + " || end: " + endIndex + " || string: " + iterateQuizInfo + " length + " + iterateQuizInfo.Length);
                iterateQuizInfo = iterateQuizInfo.Substring(startIndex + 1, endIndex);
                string[] rawAnswers = iterateQuizInfo.Split(QUIZ_ANSWER_DELIM);

                for (int i = 0; i < rawAnswers.Length; i++) {
                    if (m_verboseDebug) { Debug.Log("[LabLoad] Parsing answer: " + rawAnswers[i]); }

                    if (rawAnswers[i].ToLower().Contains(CORRECT_MARKER)) {
                        int correctIndex = i;
                        newTaskInfo.CorrectID = (uint)correctIndex;
                        if (m_verboseDebug) { Debug.Log("[LabLoad] correct index: " + correctIndex); }

                    }

                    int startAnswerIndex = rawAnswers[i].IndexOf('"');
                    rawAnswers[i] = rawAnswers[i].Substring(startAnswerIndex + 1);
                    int endAnswerIndex = rawAnswers[i].IndexOf('"');
                    string parsedAnswer = rawAnswers[i].Substring(0, endAnswerIndex);
                    newTaskInfo.SecondaryTexts.Add(parsedAnswer);
                }
                if (m_verboseDebug) { Debug.Log("[LabLoad] Question Answers: " + newTaskInfo.SecondaryTexts); }
            }
            if (m_verboseDebug) { Debug.Log("[LabLoad] Quiz Info: " + quizInfo); }

            // Efficiency
            string efficiencyInfo = sections[EFFICIENCY_INDEX].Trim();
            if (m_verboseDebug) { Debug.Log("[LabLoad] Efficiency Info: " + efficiencyInfo); }

            // TODO: labInfo.etc = etcInfo

            labInfo.Tasks.Add(newTaskInfo);
        }
    }
}