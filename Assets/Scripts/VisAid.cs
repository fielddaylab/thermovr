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
  public Transform og_transform;
  public Vector3 og_localPosition;
  public Quaternion og_localRotation;

  void Awake()
  {
    grabbable = gameObject.GetComponent<Grabbable>();
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
    if(!stored && !grabbable.grabbed)
    {
      Vector3 og_position = grabbable.og_parent.position+og_localPosition;
      gameObject.transform.position = Vector3.Lerp(gameObject.transform.position,og_position,0.1f);
      if(Vector3.Distance(gameObject.transform.position, og_position) < 0.1f)
      {
        stored = true;
        gameObject.transform.localPosition = og_localPosition;
        gameObject.transform.localRotation = og_localRotation;
      }
    }

  }
}

