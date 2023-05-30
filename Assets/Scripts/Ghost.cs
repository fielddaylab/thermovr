using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Ghost : MonoBehaviour
{
    [HideInInspector] public GameObject tool;
    public GameObject obj;

    Collider tool_c;

    public void set_tool(Tool tool) {
        this.tool = tool.gameObject;
        tool_c = tool.GetComponent<Collider>();
    }

    [System.NonSerialized]
    public bool tintersect = false;
    void OnTriggerEnter(Collider c) {
        if (c == tool_c) tintersect = true;
    }

    void OnTriggerExit(Collider c) {
        if (c == tool_c) tintersect = false;
    }

}

