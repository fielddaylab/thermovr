using System.Collections;
using System.Collections.Generic;
using UnityEngine;

/**
 * This seems to be some starting-point code for pressable buttons and things.
 * I believe this is attached to some objects on the clipboard, for quizzes or something.
 * - Luke
 **/
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

  [System.NonSerialized]
  public bool lfinger = false;
  [System.NonSerialized]
  public bool rfinger = false;
  [System.NonSerialized]
  public bool finger = false;
  void OnTriggerEnter(Collider c)
  {
    // grip squeezed                                           and finger extended                                         and collision enter //ie "the player is pointing and touched the button"
    if(OVRInput.Get(OVRInput.Axis1D.PrimaryHandTrigger)   > 0f && OVRInput.Get(OVRInput.Axis1D.PrimaryIndexTrigger)   == 0f && c == lfinger_c) { lfinger = true; on = !on; }
    if(OVRInput.Get(OVRInput.Axis1D.SecondaryHandTrigger) > 0f && OVRInput.Get(OVRInput.Axis1D.SecondaryIndexTrigger) == 0f && c == rfinger_c) { rfinger = true; on = !on; }
    finger = (lfinger || rfinger);
  }

  void OnTriggerExit(Collider c)
  {
    // grip released                                            or finger squeezed                                          or collision exit //ie "the player isn't pointing or isn't touching the button"
    if(OVRInput.Get(OVRInput.Axis1D.PrimaryHandTrigger)   == 0f || OVRInput.Get(OVRInput.Axis1D.PrimaryIndexTrigger)   > 0f || c == lfinger_c) lfinger = false;
    if(OVRInput.Get(OVRInput.Axis1D.SecondaryHandTrigger) == 0f || OVRInput.Get(OVRInput.Axis1D.SecondaryIndexTrigger) > 0f || c == rfinger_c) rfinger = false;
    finger = (lfinger || rfinger);
  }

}

