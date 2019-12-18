using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class ChallengeBall : MonoBehaviour
{
  GameObject state_dot;
  Collider state_c;

  // Start is called before the first frame update
  void Start()
  {
    state_dot = GameObject.Find("gstate");
    state_c = state_dot.GetComponent<Collider>();
  }

  // Update is called once per frame
  void Update()
  {

  }

  [System.NonSerialized]
  public bool win = false;
  void OnTriggerEnter(Collider c)
  {
    if(c == state_c) win = true;
  }

  void OnTriggerExit(Collider c)
  {
    if(c == state_c) win = false;
  }

}
