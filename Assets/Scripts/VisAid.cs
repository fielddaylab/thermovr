using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class VisAid : MonoBehaviour
{
  [System.NonSerialized]
  public bool stored = false;
  [System.NonSerialized]
  public Grabbable grabbable;
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
    grabbable = gameObject.GetComponent<Grabbable>();
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
    if(!stored && !grabbable.grabbed && t_free > 1.0f)
    {
      rigidbody.isKinematic = true;
      Vector3 og_position = grabbable.og_parent.position+og_localPosition;
      gameObject.transform.position = Vector3.Lerp(gameObject.transform.position,og_position,0.1f);
      if(Vector3.Distance(gameObject.transform.position, og_position) < 0.1f)
      {
        stored = true;
        gameObject.transform.localPosition = og_localPosition;
        gameObject.transform.localRotation = og_localRotation;
      }
    }
    else if(stored || grabbable.grabbed)
      t_free = 0.0f;
  }

}

