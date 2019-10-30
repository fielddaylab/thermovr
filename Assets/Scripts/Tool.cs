using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Tool : MonoBehaviour
{
  [System.NonSerialized]
  public bool engaged = false;
  [System.NonSerialized]
  public Grabbable grabbable;

  public GameObject storage;
  [System.NonSerialized]
  public Grabbable storage_grabbable;
  [System.NonSerialized]
  public MeshRenderer storage_available_meshrenderer;
  [System.NonSerialized]
  public MeshRenderer storage_snap_meshrenderer;

  public GameObject active;
  [System.NonSerialized]
  public Grabbable active_grabbable;
  [System.NonSerialized]
  public MeshRenderer active_available_meshrenderer;
  [System.NonSerialized]
  public MeshRenderer active_snap_meshrenderer;

  public GameObject dial;
  [System.NonSerialized]
  public Dial dial_dial;
  [System.NonSerialized]
  public Grabbable dial_grabbable;

  void Awake()
  {
    grabbable = gameObject.GetComponent<Grabbable>();

    storage_grabbable = storage.GetComponent<Grabbable>();
    storage_available_meshrenderer = storage.GetComponent<Ghost>().available.GetComponent<MeshRenderer>();
    storage_snap_meshrenderer = storage.GetComponent<Ghost>().snap.GetComponent<MeshRenderer>();

    active_grabbable = active.GetComponent<Grabbable>();
    active_available_meshrenderer = active.GetComponent<Ghost>().available.GetComponent<MeshRenderer>();
    active_snap_meshrenderer = active.GetComponent<Ghost>().snap.GetComponent<MeshRenderer>();

    dial_dial = dial.GetComponent<Dial>();;
    dial_grabbable = dial.GetComponent<Grabbable>();
  }

  // Start is called before the first frame update
  void Start()
  {

  }

  // Update is called once per frame
  void Update()
  {

  }
}
