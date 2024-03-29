using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Controls
{
    public class DesktopCameraPositioner : MonoBehaviour
    {
        [SerializeField] private Camera m_cam;
        [SerializeField] private Transform m_mainPos;
        [SerializeField] private Transform m_sliderPos;

        public void SetCameraPosition(Transform newPos)
        {
            m_cam.transform.position = newPos.position;
            m_cam.transform.rotation = newPos.rotation;
        }

        #region Editor

#if UNITY_EDITOR

        [ContextMenu("Apply Main Pos")]
        private void ApplyMainPos()
        {
            SetCameraPosition(m_mainPos);
        }

        [ContextMenu("Apply Slider Pos")]
        private void ApplySliderPos()
        {
            SetCameraPosition(m_sliderPos);
        }


#endif // UNITY_EDITOR

        #endregion // Editor
    }
}
