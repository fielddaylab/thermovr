using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class FingerToggleable : MonoBehaviour
{
  GameObject l_finger;
  GameObject r_finger;
  [System.NonSerialized]
  public bool on = false;
  Collider lfinger_c;
  Collider rfinger_c;

  void Awake()
  {
    l_finger = GameObject.Find("LFinger");
    r_finger = GameObject.Find("RFinger");
    lfinger_c = l_finger.GetComponent<Collider>();
    rfinger_c = r_finger.GetComponent<Collider>();
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
  public bool lfinger = false;
  [System.NonSerialized]
  public bool rfinger = false;
  [System.NonSerialized]
  public bool finger = false;
  void OnTriggerEnter(Collider c)
  {
    if(c == lfinger_c) lfinger = true;
    if(c == rfinger_c) rfinger = true;
    finger = (lfinger || rfinger);
    if(c == lfinger_c || c == rfinger_c) on = !on;
  }

  void OnTriggerExit(Collider c)
  {
    if(c == lfinger_c) lfinger = false;
    if(c == rfinger_c) rfinger = false;
    finger = (lfinger || rfinger);
  }

}

