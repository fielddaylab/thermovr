using BeauRoutine;
using BeauUtil;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Physics
{
    public class FixedAnchor : MonoBehaviour
    {
        [Required] [SerializeField] private GameObject m_AnchorTo;
        [SerializeField] private Axis m_Axes;

        public Vector3 GetAnchorPoint() {
            return AdjustPos();
        }

        private void Update() {
            this.transform.position = AdjustPos();
        }

        #region Helpers

        private Vector3 AdjustPos() {
            Vector3 currPos = this.transform.position;
            Vector3 anchorPos = m_AnchorTo.transform.position;
            Vector3 newPos = currPos;

            if ((m_Axes & Axis.X) != 0) {
                newPos.x = anchorPos.x;
            }
            if ((m_Axes & Axis.Y) != 0) {
                newPos.y = anchorPos.y;
            }
            if ((m_Axes & Axis.Z) != 0) {
                newPos.z = anchorPos.z;
            }

            return newPos;
        }

        #endregion // Helpers
    }
}