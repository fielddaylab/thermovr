using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Grabbable : MonoBehaviour
{
  GameObject lhand;
  GameObject rhand;
  Collider lhand_c;
  Collider rhand_c;
  public Transform og_parent;

  // Start is called before the first frame update
  void Start()
  {
    lhand = GameObject.Find("LHand");
    rhand = GameObject.Find("RHand");
    lhand_c = lhand.GetComponent<Collider>();
    rhand_c = rhand.GetComponent<Collider>();
    og_parent = gameObject.transform.parent;
  }

  // Update is called once per frame
  void Update()
  {

  }

  public bool lintersect = false;
  public bool rintersect = false;
  void OnTriggerEnter(Collider c)
  {
    if(c == lhand_c) lintersect = true;
    if(c == rhand_c) rintersect = true;
  }

  void OnTriggerExit(Collider c)
  {
    if(c == lhand_c) lintersect = false;
    if(c == rhand_c) rintersect = false;
  }

}

