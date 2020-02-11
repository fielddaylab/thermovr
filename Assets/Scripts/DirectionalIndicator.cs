using System.Collections;
using System.Collections.Generic;
using UnityEngine;

/**
 * Simple class to handle the arrow indicators for heat in/out of the system.
 * Not sure it will be useful for whatever we use in the future, but good enough for now.
 * -Luke
 **/
public class DirectionalIndicator : MonoBehaviour
{
  const int MAX_FRAME = 60;
  public Material indicator_material;
  public Transform T_left;
  public Transform T_right;
  [System.NonSerialized]
  public bool running;
  private int frame;

    // Start is called before the first frame update
  void Start()
  {
    frame = 0;
    running = true;
  }

  public void Stop()
  {
    running = false;
    MeshRenderer[] mr_list = gameObject.GetComponentsInChildren<MeshRenderer>();
    foreach (MeshRenderer mr in mr_list)
    {
      mr.enabled = false;
    }
  }

  // Update is called once per frame
  void Update()
  {
    if (running)
    {
      indicator_material.SetTextureOffset("_MainTex", new Vector2((float)(frame % MAX_FRAME) / (float)MAX_FRAME, 0.0f));
      frame++;
    }
  }

  public void FlowDirection(bool outward)
  {
    float y_rot = outward ? 90f : 0f;
    T_left.rotation = Quaternion.Euler(90f, y_rot, 0f);
    T_right.rotation = Quaternion.Euler(90f, y_rot, 0f);
  }

  public void Go()
  {
    frame = 0;
    running = true;
    MeshRenderer[] mr_list = gameObject.GetComponentsInChildren<MeshRenderer>();
    foreach (MeshRenderer mr in mr_list)
    {
      mr.enabled = true;
    }
  }
}
