using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Lightable : MonoBehaviour
{
  public bool use_custom_mats;
  public Material base_mat;
  public Material lit_mat;
  private MeshRenderer mesh_renderer;

  // Start is called before the first frame update
  void Start()
  {
    mesh_renderer = gameObject.GetComponent<MeshRenderer>();
  }

  public void SetLit(bool on)
  {
    if (on)
    {
      if (use_custom_mats)
      {
        mesh_renderer.material = lit_mat;
      }
      else
      {
        mesh_renderer.material.EnableKeyword("_Emission");
        Color color = (mesh_renderer.material.mainTexture != null) ? new Color(0.5f, 0.5f, 0.5f) : mesh_renderer.material.color;
        mesh_renderer.material.SetColor("_EmissionColor", color);
      }
    }
    else // Not On
    {
      if (use_custom_mats)
      {
        mesh_renderer.material = base_mat;
      }
      else
      {
        mesh_renderer.material.SetColor("_EmissionColor", Color.black);
      }
    }
  }
}
