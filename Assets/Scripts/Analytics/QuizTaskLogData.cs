using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Analytics
{
    [System.Serializable]
    public class QuizTaskLogData : TaskLogData
    {
        public List<string> Options;
        public List<string> SelectedOptions;
        public List<string> Answer; // Correct Answer
    }
}