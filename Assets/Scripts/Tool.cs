using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Tool : MonoBehaviour
{
  [System.NonSerialized]
  public bool engaged = false;
  [System.NonSerialized]
  public bool stored = false;
  [System.NonSerialized]
  public Grabbable grabbable;
  [System.NonSerialized]
  public Rigidbody rigidbody;
  [System.NonSerialized]
  public float t_free = 0.0f;

  public GameObject storage;
  [System.NonSerialized]
  public Ghost storage_ghost;
  [System.NonSerialized]
  public Grabbable storage_grabbable;
  [System.NonSerialized]
  public GameObject storage_available;
  [System.NonSerialized]
  public MeshRenderer storage_available_meshrenderer;
  [System.NonSerialized]
  public MeshRenderer storage_snap_meshrenderer;
  [System.NonSerialized]
  public GameObject storage_snap;

  public GameObject active;
  [System.NonSerialized]
  public Ghost active_ghost;
  [System.NonSerialized]
  public Grabbable active_grabbable;
  [System.NonSerialized]
  public GameObject active_available;
  [System.NonSerialized]
  public MeshRenderer active_available_meshrenderer;
  [System.NonSerialized]
  public MeshRenderer active_snap_meshrenderer;
  [System.NonSerialized]
  public GameObject active_snap;

  public GameObject dial;
  [System.NonSerialized]
  public Dial dial_dial;
  [System.NonSerialized]
  public Grabbable dial_grabbable;

  void Awake()
  {
    grabbable = gameObject.GetComponent<Grabbable>();
    rigidbody = gameObject.GetComponent<Rigidbody>();

    storage_ghost = storage.GetComponent<Ghost>();
    storage_grabbable = storage.GetComponent<Grabbable>();
    storage_available = storage_ghost.available;
    storage_available_meshrenderer = storage_available.GetComponent<MeshRenderer>();
    storage_snap = storage_ghost.snap;
    storage_snap_meshrenderer = storage_snap.GetComponent<MeshRenderer>();

    active_ghost = active.GetComponent<Ghost>();
    active_grabbable = active.GetComponent<Grabbable>();
    active_available = active_ghost.available;
    active_available_meshrenderer = active_available.GetComponent<MeshRenderer>();
    active_snap = active_ghost.snap;
    active_snap_meshrenderer = active_snap.GetComponent<MeshRenderer>();

    dial_dial = dial.GetComponent<Dial>();
    dial_grabbable = dial.GetComponent<Grabbable>();
  }

  // Start is called before the first frame update
  void Start()
  {

  }

  // Update is called once per frame
  void Update()
  {
    t_free += Time.deltaTime;
    if(!engaged && !stored && !grabbable.grabbed && t_free > 1.0f)
    {
      rigidbody.isKinematic = true;
      gameObject.transform.position = Vector3.Lerp(gameObject.transform.position,storage.transform.position,0.1f);
      if(Vector3.Distance(gameObject.transform.position, storage.transform.position) < 0.1f)
      {
        gameObject.transform.SetParent(storage.transform);
        stored = true;
        gameObject.transform.localPosition = new Vector3(0f,0f,0f);
        gameObject.transform.localRotation = Quaternion.identity;
      }
    }
    else if(engaged || stored || grabbable.grabbed)
      t_free = 0.0f;
  }

}

