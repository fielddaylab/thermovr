/*
DOCUMENTATION- phil, 12/16/19 [intended to be a description as of a point in time, NOT nec a prescription for how it should be evolved- feel free to uproot]

This component gets added to objects which 1. are not "Tool"s, yet 2. still act as "visual aids", generally characterised as:
"things that can be stored, grabbed, moved around, and if left disengaged, will fly back to their stored position"
*/

using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class VisAid : MonoBehaviour
{
  [System.NonSerialized]
  public bool stored = false;
  [System.NonSerialized]
  public Touchable touchable;
  [System.NonSerialized]
  public Rigidbody rigidbody;
  [System.NonSerialized]
  public Transform og_transform;
  [System.NonSerialized]
  public Vector3 og_localPosition;
  [System.NonSerialized]
  public Quaternion og_localRotation;
  [System.NonSerialized]
  public float t_free = 0.0f;

  void Awake()
  {
    touchable = gameObject.GetComponent<Touchable>();
    rigidbody = gameObject.GetComponent<Rigidbody>();
    og_transform = gameObject.transform;
    og_localPosition = gameObject.transform.localPosition;
    og_localRotation = gameObject.transform.localRotation;
  }

  // Start is called before the first frame update
  void Start()
  {

  }

  // Update is called once per frame
  void Update()
  {
    t_free += Time.deltaTime;
    if(!stored && !touchable.IsGrabbed() && t_free > 1.0f)
    {
      rigidbody.isKinematic = true;
      Vector3 og_position = touchable.og_parent.position+og_localPosition;
      gameObject.transform.position = Vector3.Lerp(gameObject.transform.position,og_position,0.1f);
      gameObject.transform.localRotation = Quaternion.Lerp(gameObject.transform.localRotation,og_localRotation,0.1f);
      if(Vector3.Distance(gameObject.transform.position, og_position) < 0.01f)
      {
        stored = true;
        gameObject.transform.parent = touchable.og_parent;
        //gameObject.transform.localPosition = og_localPosition;
        gameObject.transform.position = og_position;
        gameObject.transform.localScale = new Vector3(1f,1f,1f);
        gameObject.transform.localRotation = og_localRotation;
      }
    }
    else if(stored || touchable.IsGrabbed())
      t_free = 0.0f;
  }

}

