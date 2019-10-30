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
    meter = gameObject.transform.GetChild(0).gameObject;
  }

  // Start is called before the first frame update
  void Start()
  {

  }

  // Update is called once per frame
  void Update()
  {
    meter.transform.localRotation = Quaternion.Euler(-90.0f+val*180.0f,0.0f,0.0f);
  }

}
