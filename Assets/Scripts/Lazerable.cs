using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Lazerable : MonoBehaviour
{
  GameObject llazer;
  GameObject rlazer;
  Collider llazer_c;
  Collider rlazer_c;

  void Awake()
  {
    llazer = GameObject.Find("LLazer");
    rlazer = GameObject.Find("RLazer");
    llazer_c = llazer.GetComponent<Collider>();
    rlazer_c = rlazer.GetComponent<Collider>();
  }

  // Start is called before the first frame update
  void Start()
  {

  }

  // Update is called once per frame
  void Update()
  {

  }

  [System.NonSerialized]
  public bool lhit = false;
  [System.NonSerialized]
  public bool rhit = false;
  [System.NonSerialized]
  public bool hit = false;
  void OnTriggerEnter(Collider c)
  {
    if(c == llazer_c) lhit = true;
    if(c == rlazer_c) rhit = true;
    hit = (lhit || rhit);
  }

  void OnTriggerExit(Collider c)
  {
    if(c == llazer_c) lhit = false;
    if(c == rlazer_c) rhit = false;
    hit = (lhit || rhit);
  }

}

