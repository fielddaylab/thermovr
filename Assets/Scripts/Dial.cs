using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using TMPro;

/**
 * Runs a "dial" object, which is really a slider.
 * This controls the amount of heat I/O or weight from a tool
 **/
public class Dial : MonoBehaviour
{
  public int response_power; //sometimes, we want log mapping
  [System.NonSerialized]
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
    meter = gameObject.transform.GetChild(0).gameObject;
  }

  // Start is called before the first frame update
  void Start()
  {
  }

  // Update is called once per frame
  void Update()
  {
    Vector3 lp = meter.transform.localPosition;
    lp.x = 0.05f-val*0.1f;
    meter.transform.localPosition = lp;
    //map = min_map+(max_map-min_map)*val;
    map = response_power > 1 ? mapSharp() : mapLinear();
  }

  public void Reset()
  {
    val = 0.0f;
    Update();
  }

  /*
   * Standard way to map from 0-1 slder "val" range to min-max "tool" range.
   */
  private float mapLinear()
  { return min_map + (max_map - min_map) * val;  }

  /*
   * Non-linear mapping, which raises the 0-1 "val" to given power.
   * This gives us super-linear behavior, so for a 0-1 range we have "slow" growth of mapped value initially,
   * accelerating as we approach 1.
   * Useful for sliders whose tool's influence tends to apply too rapidly at small values.
   */
  private float mapSharp()
  { return min_map + (max_map - min_map) * Mathf.Pow(val, response_power);  }
}
