using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Fadable : MonoBehaviour
{
  //config
  public float t_in = 0.2f;
  public float t_start_out = 3.0f;
  public float t_end_out = 3.5f;
  public float t_stale = 0.5f;

  //state
  private bool factive = false;
  private bool next_factive = false;
  [System.NonSerialized]
  public float t_factive = 0.0f;
  [System.NonSerialized]
  public float t_infactive = 3.5f; //manually set to t_end_out (NOT t_end_out+t_stale- needs to be not stale for buffer window)

  //derived
  [System.NonSerialized]
  public float alpha = 0f;
  [System.NonSerialized]
  public bool stale = false;

  public void set_factive(bool a)
  {
    next_factive = a;
  }

  // Start is called before the first frame update
  void Start()
  {

  }

  // Update is called once per frame
  void Update()
  {
    if(factive != next_factive) stale = false;
    factive = next_factive;

    t_factive   += Time.deltaTime;
    t_infactive += Time.deltaTime;
    if(factive) t_infactive = 0.0f;
    else        t_factive   = 0.0f;

    if(factive)
    {
      if(t_factive < t_in)
      {
        alpha = t_factive/t_in;
        stale = false;
      }
      else
      {
        alpha = 1.0f;
        if(t_factive > t_in+t_stale) stale = true;
      }
    }
    else
    {
           if(t_infactive < t_start_out) alpha = 1.0f;
      else if(t_infactive > t_end_out)
      {
        alpha = 0.0f;
        if(t_infactive > t_end_out+t_stale) stale = true;
      }
      else
      {
        alpha = (t_end_out-t_infactive)/(t_end_out-t_start_out);
        stale = false;
      }
    }
  }
}

