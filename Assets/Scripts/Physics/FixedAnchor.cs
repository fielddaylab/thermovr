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

        private void Update() {
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

            this.transform.position = newPos;
        }
    }
}
