using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR.Dials;
using UnityEngine;

namespace ThermoVR.Controls {
    public enum Hand
    {
        LEFT,
        RIGHT,
        MOUSE // desktop
    }

    /// <summary>
    /// Interfaces with physical buttons by triggering them via raycasts
    /// </summary>
    public class DesktopVRInterfacer : MonoBehaviour
    {
        [SerializeField] private PlacementDotInteractions pdInteractions;

        private string CLICKABLE_LAYER = "Clickable";

        private GameObject m_Dragging; // the object being grabbed
        private Vector3 m_PrevWorldPos; // previous mouse position

        private bool m_nudgeMode;

        private void Update() {
            m_nudgeMode = Input.GetKey(KeyCode.LeftShift);

            if (Input.GetMouseButtonDown(0)) {
                // left button clicked
                if (RaycastFromMouse(CLICKABLE_LAYER, out GameObject objHit)) {
                    // handle dial knobs
                    if (objHit.GetComponent<Dial>()) {
                        // start dragging
                        m_Dragging = objHit;
                        m_PrevWorldPos = Camera.main.ScreenToWorldPoint(Input.mousePosition + new Vector3(0, 0, Vector3.Distance(Camera.main.transform.position, objHit.transform.position)));

                        Dial dd = objHit.GetComponent<Dial>();
                        Hand grabType = Hand.MOUSE;
                        dd.touchable.SetGrabbed(true, grabType);
                        World.Instance.GrabDial(dd, grabType);
                        EventMgr.Events.Dispatch(GameEvents.ObjectGrabbed, m_Dragging);
                    }
                    else if (objHit.GetComponent<PlacementDotInteractions>())
                    {
                        if (World.Instance.ModMgr.GraphBallInteractable())
                        {
                            pdInteractions.BeginInteract(Hand.MOUSE);
                            m_Dragging = pdInteractions.graph;
                        }
                    }
                    else if (objHit.GetComponent<InputProxy>())
                    {
                        var proxy = objHit.GetComponent<InputProxy>();
                        EventMgr.Events.Dispatch(GameEvents.InputProxySelected, proxy);
                        EventMgr.Events.Dispatch(GameEvents.EditToolValStarted, proxy.ToolType());
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
                        btnPressable.Press(false, Hand.MOUSE);
                    }
                }

                if (m_Dragging) {
                    // end dragging
                    Dial dd = m_Dragging.GetComponent<Dial>();
                    if (dd)
                    {
                        World.Instance.ReleaseDial(dd, Hand.MOUSE, true);

                        // stop grabbing
                        dd.touchable.SetGrabbed(false, Hand.MOUSE);
                    }

                    if (m_Dragging == pdInteractions.graph)
                    {
                        pdInteractions.FinishInteract(Hand.MOUSE);
                    }

                    EventMgr.Events.Dispatch(GameEvents.ObjectReleased, m_Dragging);
                    m_Dragging = null;
                }
            }

            if (m_Dragging) {
                Vector3 currPos;

                // graph ball
                if (m_Dragging == pdInteractions.graph)
                {
                    if (WorldPointFromRaycastFromMouse(CLICKABLE_LAYER, out Vector3 worldPos))
                    {
                        currPos = worldPos;
                    }
                    else
                    {
                        currPos = Vector3.zero;
                    }
                }
                else
                {
                    // dials, other
                    currPos = Camera.main.ScreenToWorldPoint(Input.mousePosition + new Vector3(0, 0, Vector3.Distance(Camera.main.transform.position, m_Dragging.transform.position)));
                }

                World.Instance.TryInteractable(ref m_Dragging, m_PrevWorldPos, ref currPos, null, Hand.MOUSE, m_nudgeMode, 0);

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

        private bool WorldPointFromRaycastFromMouse(string layer, out Vector3 worldPos)
        {
            Ray ray = Camera.main.ScreenPointToRay(Input.mousePosition);
            RaycastHit hit;
            if (UnityEngine.Physics.Raycast(ray, out hit, Mathf.Infinity, 1 << LayerMask.NameToLayer(layer)))
            {
                worldPos = hit.point;
                return true;
            }

            worldPos = Vector3.zero;
            return false;
        }

        #endregion // Helpers
    }
}
