using System.Collections;
using System.Collections.Generic;
using TMPro;
using UnityEngine;
using UnityEngine.UI;

namespace ThermoVR.Lab
{
    [RequireComponent(typeof(ExtractableData))]
    public class SensorReadout : MonoBehaviour, IExtractableData
    {
        private ExtractableData _extractableData; // extractable data component

        #region Inspector

        public TMP_Text TMP;

        [SerializeField] private Slider m_slider;

        #endregion // Inspector

        #region Unity Callbacks

        private void Awake() {
            _extractableData = this.GetComponent<ExtractableData>();
        }

        #endregion // Unity Callbacks

        #region IExtractableData

        public void ExtractData() {
            _extractableData.ExtractData();
        }

        #endregion // ExtractableData

    }
}