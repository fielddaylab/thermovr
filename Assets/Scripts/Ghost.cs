using System.Collections;
using System.Collections.Generic;
using ThermoVR;
using ThermoVR.Lab;
using ThermoVR.Tools;
using UnityEngine;

public class Ghost : MonoBehaviour
{
    private enum GhostType {
        Tool,
        Tag
    }

    [SerializeField] private GhostType m_ghostType;
    [SerializeField] private string triggerTag;
    [HideInInspector] public GameObject tool;

    public GameObject obj;

    Collider tool_c;

    private MeshRenderer m_mr;
    private Collider m_intersector;

    private void Awake() {
        m_mr = GetComponent<MeshRenderer>();

    }

    private void Start() {
        GameMgr.Events?.Register<Collider>(GameEvents.ColliderReleased, HandleColliderReleased);
    }

    public void set_tool(Tool tool) {
        this.tool = tool.gameObject;
        tool_c = tool.GetComponent<Collider>();
    }

    [System.NonSerialized]
    public bool tintersect = false;

    private void Update() {
        if (m_ghostType == GhostType.Tag) {
            if (tintersect) {
                if (m_mr.material != GameDB.Instance.SnapMat) {
                    m_mr.material = GameDB.Instance.SnapMat;
                }
            }
            else {
                if (m_mr.material != GameDB.Instance.AvailableMat) {
                    m_mr.material = GameDB.Instance.AvailableMat;
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

    #region Handlers

    /// <summary>
    /// 
    /// </summary>
    /// <param name="c"></param>
    private void HandleColliderReleased(Collider c) {
        switch (m_ghostType) {
            case GhostType.Tool:
                return;
            case GhostType.Tag:
                if (c != null && c == m_intersector) {
                    // released

                    if (c.gameObject.tag == "Cartridge") {
                        // cartridge was released onto cartridge slot

                        // Make cartridge a child of cartridge slot

                        // Activate relevant lab
                        Cartridge cartridge = c.GetComponent<Cartridge>();
                        GameMgr.Events.Dispatch(GameEvents.ActivateCartridge, cartridge);
                    }
                }
                break;
            default:
                break;
        }
    }

    #endregion // Handlers
}

