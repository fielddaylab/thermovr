using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class LabLoader : MonoBehaviour
{
    private static uint TOOL_INDEX = 1;
    private static uint SET_INDEX = 2;
    private static uint LIMIT_INDEX = 3;
    private static uint TARGET_INDEX = 4;
    private static uint QUIZ_INDEX = 5;
    private static uint EFFICIENCY_INDEX = 6;

    private static uint NUM_TASK_SECTIONS = 7; // 6 + 1 leading delim

    private static string GROUP_DELIM = "||";
    private static string TASK_INFO_DELIM = "|";

    [SerializeField] private TextAsset m_initialLab;
    [SerializeField] private bool m_verboseDebug = false;

    private void Start() {
        LoadLab(m_initialLab);
    }

    private void LoadLab(TextAsset labInfo) {
        Debug.Log("[LabLoad] Attempting to load asset " + labInfo.name + "...");

        List<string> groups = TextIO.TextAssetToList(labInfo, GROUP_DELIM);

        for (int i = 0; i < groups.Count; i++) {
            string currGroup = groups[i];

            if (currGroup.Contains("LAB-NAME")) {
                ParseLabName(currGroup);
            }
            else if (currGroup.Contains("TASK")) {
                ParseTaskInfo(currGroup);
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

        // Terminate
        Debug.Log("[LabLoad] Asset " + labInfo.name + " loaded successfully.");
    }

    private void ParseLabName(string group) {
        string labName = group.Substring(group.IndexOf(":") + 1).Trim();
        Debug.Log("[LabLoad] Lab Name: " + labName);
    }

    private void ParseTaskInfo(string group) {
        string[] sections = group.Split(TASK_INFO_DELIM);

        if (sections.Length != NUM_TASK_SECTIONS) {
            Debug.Log("[LabLoad] Task was in invalid format! Expecting " + NUM_TASK_SECTIONS + " fields, found " + sections.Length);
            return;
        }

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
        if (m_verboseDebug) { Debug.Log("[LabLoad] Quiz Info: " + quizInfo); }

        // Efficiency
        string efficiencyInfo = sections[EFFICIENCY_INDEX].Trim();
        if (m_verboseDebug) { Debug.Log("[LabLoad] Efficiency Info: " + efficiencyInfo); }
    }
}
