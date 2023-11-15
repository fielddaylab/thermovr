using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Controls {
    /// <summary>
    /// Interfaces with physical buttons by triggering them via raycasts
    /// </summary>
    public class DesktopVRInterfacer : MonoBehaviour
    {
        private GameObject m_Dragging; // the object being grabbed
        private Vector3 m_PrevWorldPos; // previous mouse position

        private void Update() {
            if (Input.GetMouseButtonDown(0)) {
                // left button clicked
                if (RaycastFromMouse(out GameObject objHit)) {
                    // handle dial knobs
                    // start dragging
                    m_Dragging = objHit;
                    m_PrevWorldPos = Camera.main.ScreenToWorldPoint(Input.mousePosition);
                }
            }
            else if (Input.GetMouseButtonUp(0)) {
                // left button released
                if (RaycastFromMouse(out GameObject objHit)) {
                    // handle buttons
                    Pressable btnPressable = objHit.GetComponent<Pressable>();
                    if (btnPressable) {
                        // trigger press
                        btnPressable.Press(false);
                    }

                    // handle dial knobs
                    // end dragging
                    m_Dragging = null;
                }
            }

            if (m_Dragging) {
                Debug.Log("Dragging");
                Vector3 currPos = Camera.main.ScreenToWorldPoint(Input.mousePosition);
                World.Instance.TryInteractable(m_Dragging, currPos, ref m_PrevWorldPos);

                m_PrevWorldPos = currPos;
            }
        }

        #region Helpers

        private bool RaycastFromMouse(out GameObject hitObj) {
            Ray ray = Camera.main.ScreenPointToRay(Input.mousePosition);
            RaycastHit hit;
            if (UnityEngine.Physics.Raycast(ray, out hit)) {
                hitObj = hit.collider.gameObject;
                return true;
            }

            hitObj = null;
            return false;
        }

        #endregion // Helpers
    }
}
