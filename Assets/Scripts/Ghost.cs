using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Ghost : MonoBehaviour
{
  public GameObject tool;
  public GameObject available;
  public GameObject snap;
  Collider tool_c;

  void Awake()
  {
    tool_c = tool.GetComponent<Collider>();
  }

  // Start is called before the first frame update
  void Start()
  {

  }

  // Update is called once per frame
  void Update()
  {

  }

  [System.NonSerialized]
  public bool tintersect = false;
  void OnTriggerEnter(Collider c)
  {
    if(c == tool_c) tintersect = true;
  }

  void OnTriggerExit(Collider c)
  {
    if(c == tool_c) tintersect = false;
  }

}

