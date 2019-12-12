using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Touchable : MonoBehaviour
{
  GameObject lhand;
  GameObject rhand;
  [System.NonSerialized]
  public bool grabbed = false;
  Collider lhand_c;
  Collider rhand_c;
  [System.NonSerialized]
  public Transform og_parent;

  void Awake()
  {
    lhand = GameObject.Find("LeftControllerAnchor");
    rhand = GameObject.Find("RightControllerAnchor");
    lhand_c = lhand.GetComponent<Collider>();
    rhand_c = rhand.GetComponent<Collider>();
    og_parent = gameObject.transform.parent;
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
  public bool ltouch = false;
  [System.NonSerialized]
  public bool rtouch = false;
  [System.NonSerialized]
  public bool touch = false;
  void OnTriggerEnter(Collider c)
  {
    if(c == lhand_c) ltouch = true;
    if(c == rhand_c) rtouch = true;
    touch = (ltouch || rtouch);
  }

  void OnTriggerExit(Collider c)
  {
    if(c == lhand_c) ltouch = false;
    if(c == rhand_c) rtouch = false;
    touch = (ltouch || rtouch);
  }

}

