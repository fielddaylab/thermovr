using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Halfable : MonoBehaviour
{
  public GameObject full;
  public GameObject half;
  [System.NonSerialized]
  public MeshRenderer full_meshrenderer;
  [System.NonSerialized]
  public MeshRenderer half_meshrenderer;
  [System.NonSerialized]
  public bool halved = false;

  public void setHalf(bool h)
  {
    if(h == halved) return;
    halved = h;
    if(halved)
    {
      full_meshrenderer.enabled = false;
      half_meshrenderer.enabled = true;
    }
    else
    {
      full_meshrenderer.enabled = true;
      half_meshrenderer.enabled = false;
    }
  }

  // Start is called before the first frame update
  void Start()
  {
    full_meshrenderer = full.GetComponent<MeshRenderer>();
    half_meshrenderer = half.GetComponent<MeshRenderer>();
    half_meshrenderer.enabled = false;
    halved = false;
  }

  // Update is called once per frame
  void Update()
  {

  }

}

