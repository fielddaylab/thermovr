using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using TMPro;

/**
 * Presumably used for tabs on the clipboard.
 * Doesn't seem to be in use.
 **/
public class Tab : MonoBehaviour
{

  public TextMeshPro tmp;
  public GameObject backing;

  [System.NonSerialized]
  public FingerToggleable fingertoggleable;
  [System.NonSerialized]
  public MeshRenderer backing_meshrenderer;

  void Awake()
  {
    fingertoggleable = GetComponent<FingerToggleable>();
    backing_meshrenderer = backing.GetComponent<MeshRenderer>();
  }

  // Start is called before the first frame update
  void Start()
  {
  }

  // Update is called once per frame
  void Update()
  {
  }

}

