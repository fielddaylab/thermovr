using System.Collections;
using System.Collections.Generic;
using System.Linq;
using TMPro;
using UnityEngine;

namespace ThermoVR.Lab
{
    [RequireComponent(typeof(Touchable))]
    public class Cartridge : MonoBehaviour
    {
        private Touchable m_touchable;

        [SerializeField] private string m_initialLabel;

        #region Inspector 

        [SerializeField] private TextMeshPro[] m_labels;

        #endregion // Inspector

        private bool m_beingGrabbed; // if object is currently being grabbed

        #region Unity Callbacks

        private void Awake() {
            m_touchable = GetComponent<Touchable>();
            m_beingGrabbed = false;
        }

        private void Start() {
            GameMgr.Events.Dispatch(GameEvents.RegisterMovable, this.GetComponent<Touchable>());

            if (m_initialLabel != "Label") {
                SetLabels(m_initialLabel);
            }
        }

        private void Update() {
            bool grabChange = false;

            if (m_beingGrabbed) {
                if (!m_touchable.grabbed) {
                    m_beingGrabbed = false;
                    grabChange = true;
                }
            }
            else {
                if (m_touchable.grabbed) {
                    m_beingGrabbed = true;
                    grabChange = true;
                }
            }

            if (grabChange) {
                // update tablet ghost accordingly
                // pass grab state and specific cartridge
            }
        }
        
        #endregion // Unity Callbacks

        public void SetInfo(LabInfo newInfo) {
            SetLabels(newInfo.Name);
        }

        private void SetLabels(string newLabel) {
            for (int i = 0; i < m_labels.Length; i++) {
                m_labels[i].SetText(newLabel);
            }
        }
    }
}
