using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Pressable : Touchable
{
  const float BUTTONS_BASE_HEIGHT = 0.0115f;

  //GameObject lhand;
  //GameObject rhand;
  //[System.NonSerialized]
  //public bool grabbed = false;
  //protected Collider[] lhand_c;
  //protected Collider[] rhand_c;
  //[System.NonSerialized]
  //public Transform og_parent;

  void Awake()
  {
    lhand = GameObject.Find("LeftControllerAnchor");
    rhand = GameObject.Find("RightControllerAnchor");
    lhand_c = lhand.GetComponentsInChildren<Collider>();
    rhand_c = rhand.GetComponentsInChildren<Collider>();
    og_parent = gameObject.transform.parent;
  }

  private bool cInHandList(Collider c, bool left)
  {
    Collider[] list = left ? lhand_c : rhand_c;
    foreach (Collider candidate in list)
    {
      if (candidate == c) return true;
    }
    return false;
  }

  //[System.NonSerialized]
  //public bool ltouch = false;
  //[System.NonSerialized]
  //public bool rtouch = false;
  //[System.NonSerialized]
  //public bool touch = false;
  protected override void OnTriggerEnter(Collider c)
  {
    if (OVRInput.Get(OVRInput.RawAxis1D.LHandTrigger) > 0f && OVRInput.Get(OVRInput.RawAxis1D.LIndexTrigger) == 0f && cInHandList(c, true)) { ltouch = true; }
    if (OVRInput.Get(OVRInput.RawAxis1D.RHandTrigger) > 0f && OVRInput.Get(OVRInput.RawAxis1D.RIndexTrigger) == 0f && cInHandList(c, false)) { rtouch = true; }
    touch = (ltouch || rtouch);
    if (touch)
    {
      gameObject.GetComponentInChildren<MeshRenderer>().gameObject.transform.localPosition = new Vector3(0.0f, 0.5f * BUTTONS_BASE_HEIGHT, 0.0f);
      //var trans = gameObject.GetComponentInChildren<MeshRenderer>().gameObject.transform;
      //trans.localPosition = new Vector3(trans.position.x, 0.5f * BUTTONS_BASE_HEIGHT, trans.position.z);
    }
    else
    {
      gameObject.GetComponentInChildren<MeshRenderer>().gameObject.transform.localPosition = new Vector3(0.0f, BUTTONS_BASE_HEIGHT, 0.0f);
      //var trans = gameObject.GetComponentInChildren<MeshRenderer>().gameObject.transform;
      //trans.position = new Vector3(trans.position.x, BUTTONS_BASE_HEIGHT, trans.position.z);
    }
  }

  protected override void OnTriggerExit(Collider c)
  {
    if ((OVRInput.Get(OVRInput.Axis1D.PrimaryHandTrigger) == 0f || OVRInput.Get(OVRInput.Axis1D.PrimaryIndexTrigger) > 0f) && cInHandList(c, true)) ltouch = false;
    if ((OVRInput.Get(OVRInput.Axis1D.SecondaryHandTrigger) == 0f || OVRInput.Get(OVRInput.Axis1D.SecondaryIndexTrigger) > 0f) && cInHandList(c, false)) rtouch = false;
    touch = (ltouch || rtouch);
    if (!touch)
    {
      gameObject.GetComponentInChildren<MeshRenderer>().gameObject.transform.localPosition = new Vector3(0.0f, BUTTONS_BASE_HEIGHT, 0.0f);
      //var trans = gameObject.GetComponentInChildren<MeshRenderer>().gameObject.transform;
      //trans.position = new Vector3(trans.position.x, BUTTONS_BASE_HEIGHT, trans.position.z);
    }
  }
  //GameObject l_finger;
  //GameObject r_finger;
  //[System.NonSerialized]
  //Collider lfinger_c;
  //Collider rfinger_c;

  //void Awake()
  //{
  //  l_finger = GameObject.Find("LFinger");
  //  r_finger = GameObject.Find("RFinger");
  //  lfinger_c = l_finger.GetComponent<Collider>();
  //  rfinger_c = r_finger.GetComponent<Collider>();
  //}

  //// Start is called before the first frame update
  //void Start()
  //{

  //}

  //// Update is called once per frame
  //void Update()
  //{

  //}

  //[System.NonSerialized]
  //public bool lfinger = false;
  //[System.NonSerialized]
  //public bool rfinger = false;
  //[System.NonSerialized]
  //public bool pressed = false;
  //void OnTriggerEnter(Collider c)
  //{
  //  // grip squeezed                                           and finger extended                                         and collision enter //ie "the player is pointing and touched the button"
  //  if(OVRInput.Get(OVRInput.Axis1D.PrimaryHandTrigger)   > 0f && OVRInput.Get(OVRInput.Axis1D.PrimaryIndexTrigger)   == 0f && c == lfinger_c) { lfinger = true; }
  //  if(OVRInput.Get(OVRInput.Axis1D.SecondaryHandTrigger) > 0f && OVRInput.Get(OVRInput.Axis1D.SecondaryIndexTrigger) == 0f && c == rfinger_c) { rfinger = true; }
  //  pressed = (lfinger || rfinger);
  //  if (pressed)
  //  {
  //    gameObject.GetComponentInChildren<MeshRenderer>().gameObject.transform.localPosition = new Vector3(0.0f, 0.5f * BUTTONS_BASE_HEIGHT, 0.0f);
  //    //var trans = gameObject.GetComponentInChildren<MeshRenderer>().gameObject.transform;
  //    //trans.localPosition = new Vector3(trans.position.x, 0.5f * BUTTONS_BASE_HEIGHT, trans.position.z);
  //  }
  //  else
  //  {
  //    gameObject.GetComponentInChildren<MeshRenderer>().gameObject.transform.localPosition = new Vector3(0.0f, BUTTONS_BASE_HEIGHT, 0.0f);
  //    //var trans = gameObject.GetComponentInChildren<MeshRenderer>().gameObject.transform;
  //    //trans.position = new Vector3(trans.position.x, BUTTONS_BASE_HEIGHT, trans.position.z);
  //  }
  //}

  //void OnTriggerExit(Collider c)
  //{
  //  // grip released                                            or finger squeezed                                          or collision exit //ie "the player isn't pointing or isn't touching the button"
  //  if(OVRInput.Get(OVRInput.Axis1D.PrimaryHandTrigger)   == 0f || OVRInput.Get(OVRInput.Axis1D.PrimaryIndexTrigger)   > 0f || c == lfinger_c) lfinger = false;
  //  if(OVRInput.Get(OVRInput.Axis1D.SecondaryHandTrigger) == 0f || OVRInput.Get(OVRInput.Axis1D.SecondaryIndexTrigger) > 0f || c == rfinger_c) rfinger = false;
  //  pressed = (lfinger || rfinger);
  //  if (!pressed)
  //  {
  //    gameObject.GetComponentInChildren<MeshRenderer>().gameObject.transform.localPosition = new Vector3(0.0f, BUTTONS_BASE_HEIGHT, 0.0f);
  //    //var trans = gameObject.GetComponentInChildren<MeshRenderer>().gameObject.transform;
  //    //trans.position = new Vector3(trans.position.x, BUTTONS_BASE_HEIGHT, trans.position.z);
  //  }
  //}
}
