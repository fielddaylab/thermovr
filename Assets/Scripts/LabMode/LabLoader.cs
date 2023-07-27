using System.Collections;
using System.Collections.Generic;
using System.IO;
using UnityEngine;

namespace ThermoVR.Lab
{
    public struct LabInfo {
        public string Name;
    }


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

        [SerializeField] private TextAsset[] m_initialLabs;
        [SerializeField] private bool m_verboseDebug = false;

        private void Start() {
            for (int i = 0; i < m_initialLabs.Length; i++) {
                LoadLab(m_initialLabs[i]);
            }

            GameMgr.Events?.Register<Cartridge>(GameEvents.ActivateCartridge, HandleActivateCartridge);

        }

        private void LoadLab(TextAsset labInfoAsset) {
            Debug.Log("[LabLoad] Attempting to load asset " + labInfoAsset.name + "...");

            bool succeeded = true;
            LabInfo newInfo = new LabInfo();

            try {
                List<string> groups = TextIO.TextAssetToList(labInfoAsset, GROUP_DELIM);

                for (int i = 0; i < groups.Count; i++) {
                    string currGroup = groups[i];

                    if (currGroup.Contains("LAB-NAME")) {
                        ParseLabName(currGroup, ref newInfo);
                    }
                    else if (currGroup.Contains("TASK")) {
                        ParseTaskInfo(currGroup, ref newInfo);
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

                GameMgr.Events.Dispatch(GameEvents.LabLoaded, newInfo);
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
            string[] sections = group.Split(TASK_INFO_DELIM);

            if (sections.Length != NUM_TASK_SECTIONS) {
                Debug.Log("[LabLoad] Task was in invalid format! Expecting " + NUM_TASK_SECTIONS + " fields, found " + sections.Length);
                throw new InvalidDataException();
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

            // TODO: labInfo.etc = etcInfo
        }


        #region Handlers

        private void HandleActivateCartridge(Cartridge cartridge) {
            switch (cartridge.GetCartridgeType()) {
                case Cartridge.CartridgeType.Lab:
                    Debug.Log("[Cartridge] Lab " + cartridge.GetInfo().Name + " loaded!");
                    break;
                case Cartridge.CartridgeType.Sandbox:
                    Debug.Log("[Cartridge] Sandbox loaded!");
                    break;
                default:
                    break;
            }
        }

        #endregion // Handlers
    }
}