using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Touchable : MonoBehaviour
{
  protected GameObject lhand;
  protected GameObject rhand;
  [System.NonSerialized]
  public bool grabbed = false;
  protected Collider[] lhand_c;
  protected Collider[] rhand_c;
  [System.NonSerialized]
  public Transform og_parent;

  void Awake()
  {
    lhand = GameObject.Find("LeftControllerAnchor");
    rhand = GameObject.Find("RightControllerAnchor");
    lhand_c = lhand.GetComponentsInChildren<Collider>();
    rhand_c = rhand.GetComponentsInChildren<Collider>();
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

  private bool cInHandList(Collider c, bool left)
  {
    Collider[] list = left ? lhand_c : rhand_c;
    foreach (Collider candidate in list)
    {
      if (candidate == c) return true;
    }
    return false;
  }

  [System.NonSerialized]
  public bool ltouch = false;
  [System.NonSerialized]
  public bool rtouch = false;
  [System.NonSerialized]
  public bool touch = false;
  protected virtual void OnTriggerEnter(Collider c)
  {
    if (cInHandList(c, true))  ltouch = true;
    if (cInHandList(c, false)) rtouch = true;

    touch = (ltouch || rtouch);

    if (touch)
    {
      var light_list = gameObject.GetComponentsInChildren<Lightable>();
      foreach (Lightable light in light_list)
      {
        light.SetLit(touch); // set to lit if touched, or unlit if not.
      }
    }
  }

  protected virtual void OnTriggerExit(Collider c)
  {
    if(cInHandList(c, true))  ltouch = false;
    if(cInHandList(c, false)) rtouch = false;

    touch = (ltouch || rtouch);

    if (!touch)
    {
      var light_list = gameObject.GetComponentsInChildren<Lightable>();
      foreach (Lightable light in light_list)
      {
        light.SetLit(touch); // set to lit if touched, or unlit if not.
      }
    }
  }

}

