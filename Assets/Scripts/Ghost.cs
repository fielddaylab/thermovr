using System.Collections;
using System.Collections.Generic;
using ThermoVR;
using ThermoVR.Lab;
using ThermoVR.Tools;
using UnityEngine;

public class Ghost : MonoBehaviour
{
    private enum GhostType
    {
        Tool,
        Tag
    }

    [SerializeField] private GhostType m_ghostType;
    [SerializeField] private string triggerTag;
    [HideInInspector] public GameObject tool;

    [System.NonSerialized]
    public bool tintersect = false;

    public GameObject obj;

    Collider tool_c;

    private MeshRenderer[] m_mrs;
    private BoxCollider[] m_bcs;
    private Collider m_intersector;
    private bool m_occupied;

    #region Unity Callbacks

    private void Awake() {
        m_mrs = GetComponentsInChildren<MeshRenderer>();
        m_bcs = GetComponentsInChildren<BoxCollider>();

        if (m_ghostType == GhostType.Tag) {
            // Disable initially by default, until a relevant tool is picked up
            SetEnabled(false);
            m_occupied = false;
        }
    }

    private void Start() {
        GameMgr.Events?.Register<Collider>(GameEvents.ColliderReleased, HandleColliderReleased);
        GameMgr.Events?.Register<Collider>(GameEvents.ColliderGrabbed, HandleColliderGrabbed);
    }

    public void set_tool(Tool tool) {
        this.tool = tool.gameObject;
        tool_c = tool.GetComponent<Collider>();
    }

    private void Update() {
        if (m_ghostType == GhostType.Tag) {
            if (tintersect) {
                if (CurrMat() != GameDB.Instance.SnapMat) {
                    SetMat(GameDB.Instance.SnapMat);
                }
            }
            else {
                if (CurrMat() != GameDB.Instance.AvailableMat) {
                    SetMat(GameDB.Instance.AvailableMat);
                }
            }
        }
    }

    void OnTriggerEnter(Collider c) {
        switch (m_ghostType) {
            case GhostType.Tool:
                if (c == tool_c) tintersect = true;
                break;
            case GhostType.Tag:
                if (c.gameObject.tag == triggerTag) {
                    tintersect = true;
                    m_intersector = c;
                }
                break;
            default:
                break;
        }

    }

    void OnTriggerExit(Collider c) {
        switch (m_ghostType) {
            case GhostType.Tool:
                if (c == tool_c) tintersect = false;
                break;
            case GhostType.Tag:
                if (c.gameObject.tag == triggerTag) {
                    tintersect = false;
                    m_intersector = null;
                }
                break;
            default:
                break;
        }
    }

    #endregion // Unity Callbacks

    #region Handlers

    private void HandleColliderReleased(Collider c) {
        switch (m_ghostType) {
            case GhostType.Tool:
                return;
            case GhostType.Tag:
                if (c != null) {
                    if (c.gameObject.tag == triggerTag) {
                        if (c == m_intersector) {
                            // released
                            // cartridge was released onto cartridge slot

                            // Make cartridge a child of cartridge slot
                            c.transform.SetParent(this.transform);
                            c.transform.localPosition = new Vector3(0f, 0f, 0f);
                            c.transform.localRotation = Quaternion.identity;
                            c.transform.localScale = new Vector3(1f, 1f, 1f);

                            // Hide this ghost
                            SetEnabled(false);
                            m_occupied = true;

                            // Activate relevant lab
                            Cartridge cartridge = c.GetComponent<Cartridge>();
                            GameMgr.Events.Dispatch(GameEvents.ActivateCartridge, cartridge);
                        }
                        else {
                            // Hide this ghost
                            SetEnabled(false);
                        }
                    }
                }
                break;
            default:
                break;
        }
    }

    private void HandleColliderGrabbed(Collider c) {
        switch (m_ghostType) {
            case GhostType.Tool:
                return;
            case GhostType.Tag:
                if (c != null) {
                    if (c.gameObject.tag == triggerTag) {
                        if (c == m_intersector) {
                            // Previously inserted cartridge has been removed

                            Debug.Log("[Cartridge] Cartridge removed");

                            // Re-enable active ghost
                            SetEnabled(true);
                            m_occupied = false;

                            Cartridge cartridge = c.GetComponent<Cartridge>();
                            GameMgr.Events.Dispatch(GameEvents.DeactivateCartridge, cartridge);
                        }
                        else {
                            // Highlight active ghost if not occupied
                            if (!m_occupied) {
                                SetEnabled(true);
                            }
                        }
                    }
                }
                break;
            default:
                break;
        }
    }

    #endregion // Handlers

    /// <summary>
    /// Sets the mesh renderers and box colliders of ghost to enabled/disabled state
    /// </summary>
    /// <param name="enabled"></param>
    private void SetEnabled(bool enabled) {
        for (int i = 0; i < m_mrs.Length; i++) {
            m_mrs[i].enabled = enabled;
        }
        for (int i = 0; i < m_bcs.Length; i++) {
            m_bcs[i].enabled = enabled;
        }
    }

    private Material CurrMat() {
        return m_mrs[0].material;
    }

    private void SetMat(Material mat) {
        for (int i = 0; i < m_mrs.Length; i++) {
            m_mrs[i].material = mat;
        }
    }
}

