using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Flasher : MonoBehaviour
{
  private const int PERIOD = 15;

  public Lightable flash_light;

  private int frame;
  private bool flashing;

  // Start is called before the first frame update
  void Start()
  {
    frame = 0;
  }

  // Update is called once per frame
  void Update()
  {
    if (flashing)
    {
      if (frame == 0)
      {
        flash_light.SetLit(true);
      }
      else if (frame == PERIOD)
      {
        flash_light.SetLit(false);
      }
      frame = (frame + 1) % (2 * PERIOD);
    }
  }

  public void Stop()
  {
    if (flashing)
    {
      flashing = false;
    }
  }

  public void Flash()
  {
    if (!flashing)
    {
      flashing = true;
      frame = 0;
    }
  }
}
