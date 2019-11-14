using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using TMPro;

public class Dial : MonoBehaviour
{
  public float val = 0.0f;
  [System.NonSerialized]
  public float prev_val = 0.0f;
  [System.NonSerialized]
  public float map = 0.0f;
  [System.NonSerialized]
  public float min_map = 0.0f;
  [System.NonSerialized]
  public float max_map = 1.0f;
  [System.NonSerialized]
  public string unit = "";
  public GameObject tool;

  [System.NonSerialized]
  public bool examined = false;
  GameObject meter;

  void Awake()
  {
    meter = gameObject.transform.GetChild(1).gameObject;
  }

  // Start is called before the first frame update
  void Start()
  {
  }

  // Update is called once per frame
  void Update()
  {
    Vector3 lp = meter.transform.localPosition;
    lp.z = -0.25f+val/2.0f;
    meter.transform.localPosition = lp;
    map = min_map+(max_map-min_map)*val;
  }

}
