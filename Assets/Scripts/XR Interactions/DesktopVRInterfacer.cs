using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR.Dials;
using UnityEngine;

namespace ThermoVR.Controls {
    /// <summary>
    /// Interfaces with physical buttons by triggering them via raycasts
    /// </summary>
    public class DesktopVRInterfacer : MonoBehaviour
    {
        private string CLICKABLE_LAYER = "Clickable";

        private GameObject m_Dragging; // the object being grabbed
        private Vector3 m_PrevWorldPos; // previous mouse position

        private void Update() {
            if (Input.GetMouseButtonDown(0)) {
                // left button clicked
                if (RaycastFromMouse(CLICKABLE_LAYER, out GameObject objHit)) {
                    // handle dial knobs
                    if (objHit.GetComponent<Dial>()) {
                        // start dragging
                        m_Dragging = objHit;
                        m_PrevWorldPos = Camera.main.ScreenToWorldPoint(Input.mousePosition + new Vector3(0, 0, Vector3.Distance(Camera.main.transform.position, objHit.transform.position)));

                        GameMgr.Events.Dispatch(GameEvents.ObjectGrabbed, m_Dragging);
                    }

                }
            }
            else if (Input.GetMouseButtonUp(0)) {
                // left button released
                if (RaycastFromMouse(CLICKABLE_LAYER, out GameObject objHit)) {
                    // handle buttons
                    Pressable btnPressable = objHit.GetComponent<Pressable>();
                    if (btnPressable) {
                        // trigger press
                        btnPressable.Press(false, true, false);
                    }

                    // handle dial knobs

                }

                if (m_Dragging) {
                    // end dragging
                    GameMgr.Events.Dispatch(GameEvents.ObjectReleased, m_Dragging);
                    m_Dragging = null;
                }
            }

            if (m_Dragging) {
                Vector3 currPos = Camera.main.ScreenToWorldPoint(Input.mousePosition + new Vector3(0, 0, Vector3.Distance(Camera.main.transform.position, m_Dragging.transform.position)));
                World.Instance.TryInteractable(ref m_Dragging, m_PrevWorldPos, ref currPos, null, true);

                m_PrevWorldPos = currPos;
            }
        }

        #region Helpers

        private bool RaycastFromMouse(string layer, out GameObject hitObj) {
            Ray ray = Camera.main.ScreenPointToRay(Input.mousePosition);
            RaycastHit hit;
            if (UnityEngine.Physics.Raycast(ray, out hit, Mathf.Infinity, 1 << LayerMask.NameToLayer(layer))) {
                hitObj = hit.collider.gameObject;
                return true;
            }

            hitObj = null;
            return false;
        }

        #endregion // Helpers
    }
}
