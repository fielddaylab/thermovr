using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using TMPro;

public class Quizo : MonoBehaviour
{
  [System.NonSerialized]
  public Lazerable lazerable;

  public TextMeshPro tmp;
  [System.NonSerialized]
  public TextContainer tc;

  void Awake()
  {
    lazerable = GetComponent<Lazerable>();
    tc = tmp.GetComponent<TextContainer>();
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

