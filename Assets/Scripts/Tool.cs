/*
DOCUMENTATION- phil, 12/16/19 [intended to be a description as of a point in time, NOT nec a prescription for how it should be evolved- feel free to uproot]

The various tools which can be variably engaged to the container, thrown, or stored.
*/

using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using TMPro;

public class Tool : MonoBehaviour
{
  [System.NonSerialized]
  public bool engaged = false;
  [System.NonSerialized]
  public bool stored = false;
  [System.NonSerialized]
  public Touchable touchable;
  [System.NonSerialized]
  public BoxCollider boxcollider;
  [System.NonSerialized]
  public Rigidbody rigidbody;
  [System.NonSerialized]
  public float t_free = 0.0f;

  public GameObject mesh;

  public GameObject storage;
  [System.NonSerialized]
  public float default_storage_scale;
  [System.NonSerialized]
  public Ghost storage_ghost;
  [System.NonSerialized]
  public Touchable storage_touchable;
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
  public Touchable active_touchable;
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
  public Touchable dial_touchable;

  public GameObject text;
  public GameObject textv;
  [System.NonSerialized]
  public TextMeshPro textv_tmpro;
  [System.NonSerialized]
  public MeshRenderer textv_meshrenderer;

  public GameObject textl;
  [System.NonSerialized]
  public TextMeshPro textl_tmpro;
  [System.NonSerialized]
  public MeshRenderer textl_meshrenderer;

  public GameObject textd;
  [System.NonSerialized]
  public TextMeshPro textd_tmpro;
  [System.NonSerialized]
  public MeshRenderer textd_meshrenderer;

  public GameObject textn;

  [System.NonSerialized]
  public Fadable text_fadable;
  [System.NonSerialized]
  public bool disabled;

  void Awake()
  {
    touchable = gameObject.GetComponent<Touchable>();
    boxcollider = gameObject.GetComponent<BoxCollider>();
    rigidbody = gameObject.GetComponent<Rigidbody>();

    default_storage_scale = storage.transform.localScale.x; //could grab any dimension
    storage_ghost = storage.GetComponent<Ghost>();
    storage_touchable = storage.GetComponent<Touchable>();
    storage_available = storage_ghost.available;
    storage_available_meshrenderer = storage_available.GetComponent<MeshRenderer>();
    storage_snap = storage_ghost.snap;
    storage_snap_meshrenderer = storage_snap.GetComponent<MeshRenderer>();

    active_ghost = active.GetComponent<Ghost>();
    active_touchable = active.GetComponent<Touchable>();
    active_available = active_ghost.available;
    active_available_meshrenderer = active_available.GetComponent<MeshRenderer>();
    active_snap = active_ghost.snap;
    active_snap_meshrenderer = active_snap.GetComponent<MeshRenderer>();

    dial_dial = dial.GetComponent<Dial>();
    dial_touchable = dial.GetComponent<Touchable>();

    textv_tmpro = textv.GetComponent<TextMeshPro>();
    textv_meshrenderer = textv.GetComponent<MeshRenderer>();
    textl_tmpro = textl.GetComponent<TextMeshPro>();
    textl_meshrenderer = textl.GetComponent<MeshRenderer>();
    text_fadable = GetComponent<Fadable>();
    if(textd)
    {
      textd_tmpro = textd.GetComponent<TextMeshPro>();
      textd_meshrenderer = textd.GetComponent<MeshRenderer>();
    }

    disabled = false;
  }

  // Start is called before the first frame update
  void Start()
  {

  }

  // Update is called once per frame
  void Update()
  {
    t_free += Time.deltaTime;

    if (text_fadable) text_fadable.set_factive(touchable.touch || dial_dial.examined);

    // If tool has been sitting unused and unmoved for long enough, start it moving back to storage.
    if(!engaged && !stored && !touchable.grabbed && t_free > 1.0f)
    {
      rigidbody.isKinematic = true;
      gameObject.transform.position = Vector3.Lerp(gameObject.transform.position,storage.transform.position,0.1f);
      gameObject.transform.localRotation = Quaternion.Lerp(gameObject.transform.localRotation,storage.transform.localRotation,0.1f);
      if(Vector3.Distance(gameObject.transform.position, storage.transform.position) < 0.01f)
      {
        gameObject.transform.SetParent(storage.transform);
        stored = true;
        gameObject.transform.localPosition = new Vector3(0f,0f,0f);
        gameObject.transform.localScale = new Vector3(1f,1f,1f);
        gameObject.transform.localRotation = Quaternion.identity;
        float v = storage.transform.localScale.x; //can grab any dimension
        Vector3 invscale = new Vector3(1f/v,1f/v,1f/v);
        text.transform.localScale = invscale;
      }
    }
    else if(engaged || stored || touchable.grabbed)
      t_free = 0.0f;
  }

}

