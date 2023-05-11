using System.Collections;
using System.Collections.Generic;
using UnityEngine;

/**
 * Class to allow "lighting" of an object.
 * By default, just adds emission uniformly to the object, based on material/texture color.
 * Optionally, one can provide "lit" and "unlit" materials, instead of using a naive "emission" approach to lighting.
 **/
[RequireComponent(typeof(MeshRenderer))]
public class Lightable : MonoBehaviour
{
    public bool use_custom_mats;
    public Material base_mat;
    public Material lit_mat;
    private MeshRenderer mesh_renderer;

    private bool is_on;

    // Start is called before the first frame update
    void Start()
    {
        mesh_renderer = gameObject.GetComponent<MeshRenderer>();
        is_on = false;
    }

    public void SetLit(bool on)
    {
        if (mesh_renderer == null)
        {
            mesh_renderer = gameObject.GetComponent<MeshRenderer>();
        }

        if (on && !is_on) // turn on
        {
            if (use_custom_mats)
            {
                mesh_renderer.materials[0] = lit_mat;
            }
            else
            {
                mesh_renderer.materials[0].EnableKeyword("_Emission");
                // if there is a texture, use grey light, which will be modified by texture to match color (it appears).
                // else, emit light matching exactly to the material color (which wouldn't work so well if we had a black object, I suppose...)
                Color color = (mesh_renderer.materials[0].mainTexture != null) ? new Color(0.5f, 0.5f, 0.5f) : mesh_renderer.materials[0].color;
                mesh_renderer.materials[0].SetColor("_EmissionColor", color);
            }
            is_on = true;
        }
        else if (!on && is_on) // turn off
        {
            if (use_custom_mats)
            {
                Debug.Log(this.gameObject.name);
                mesh_renderer.materials[0] = base_mat;
            }
            else
            {
                mesh_renderer.materials[0].SetColor("_EmissionColor", Color.black);
            }
            is_on = false;
        }
    }
}
