using System.Collections;
using System.Collections.Generic;
using Unity.XR.CoreUtils;
using UnityEngine;
using UnityEngine.XR;

namespace ThermoVR.Controls
{
    public class XRHeightController : MonoBehaviour
    {
        [SerializeField] private float m_InitialHeight = 0.325f;
        [SerializeField] private Transform m_CenterEye;
        [SerializeField] private XROrigin m_Origin;

        private bool initialSet = false;

        private void Update() {
            InputDevices.GetDeviceAtXRNode(XRNode.Head).TryGetFeatureValue(UnityEngine.XR.CommonUsages.userPresence, out bool isPresent);

            if (!initialSet && isPresent) {
                SetHeight();

                initialSet = true;
            }
        }

        private void SetHeight() {
            Vector3 currPos = this.transform.position;

            this.transform.position = new Vector3(currPos.x, m_InitialHeight - (m_CenterEye.transform.localPosition.y - m_Origin.CameraYOffset), currPos.z);
        }
    }

}