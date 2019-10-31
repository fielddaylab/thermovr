using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Dial : MonoBehaviour
{
  public float val = 0.0f;
  public GameObject tool;
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
    lp.z = -0.5f+val;
    meter.transform.localPosition = lp;
  }

}
