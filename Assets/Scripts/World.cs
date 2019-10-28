using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class World : MonoBehaviour
{
  List<TextMesh> DEBUGTEXTS;

  ThermoMath thermo;
  GameObject lhand;
  GameObject rhand;

  public Material hand_empty;
  public Material hand_intersecting;
  public Material hand_grabbing;

  List<Grabbable> movables;
  GameObject lgrabbed = null;
  GameObject rgrabbed = null;
  float ltwist = 0.0f;
  float rtwist = 0.0f;

  List<Tool> tools;
  List<Grabbable> dials;

  bool lhtrigger = false;
  bool rhtrigger = false;
  bool litrigger = false;
  bool ritrigger = false;

  //buttons
  GameObject TInc;
  GameObject TDec;
  GameObject PInc;
  GameObject PDec;
  GameObject VInc;
  GameObject VDec;

  // Start is called before the first frame update
  void Start()
  {
    DEBUGTEXTS = new List<TextMesh>();
    GameObject dtexts = GameObject.Find("DEBUGTEXTS");
    foreach(Transform child in dtexts.transform)
      DEBUGTEXTS.Add(child.gameObject.GetComponent<TextMesh>());

    thermo = GameObject.Find("Oracle").GetComponent<ThermoMath>();

    lhand  = GameObject.Find("LHand");
    rhand  = GameObject.Find("RHand");

    lhand.GetComponent<MeshRenderer>().material = hand_empty;
    rhand.GetComponent<MeshRenderer>().material = hand_empty;

    tools = new List<Tool>();
    tools.Add(GameObject.Find("Tool_Insulator").GetComponent<Tool>());
    tools.Add(GameObject.Find("Tool_Clamp").GetComponent<Tool>());
    tools.Add(GameObject.Find("Tool_Burner").GetComponent<Tool>());
    tools.Add(GameObject.Find("Tool_Coil").GetComponent<Tool>());
    tools.Add(GameObject.Find("Tool_Weight").GetComponent<Tool>());
    tools.Add(GameObject.Find("Tool_Balloon").GetComponent<Tool>());

    dials = new List<Grabbable>();
    for(int i = 0; i < tools.Count; i++) dials.Add(tools[i].dial.GetComponent<Grabbable>());

    movables = new List<Grabbable>();
    for(int i = 0; i < tools.Count; i++) movables.Add(tools[i].gameObject.GetComponent<Grabbable>()); //important that tools take priority, so they can be grabbed and removed
    movables.Add(GameObject.Find("Graph").GetComponent<Grabbable>());
    movables.Add(GameObject.Find("Vessel").GetComponent<Grabbable>());

    //buttons
    TInc = GameObject.Find("TInc");
    movables.Add(TInc.GetComponent<Grabbable>());
    TDec = GameObject.Find("TDec");
    movables.Add(TDec.GetComponent<Grabbable>());
    PInc = GameObject.Find("PInc");
    movables.Add(PInc.GetComponent<Grabbable>());
    PDec = GameObject.Find("PDec");
    movables.Add(PDec.GetComponent<Grabbable>());
    VInc = GameObject.Find("VInc");
    movables.Add(VInc.GetComponent<Grabbable>());
    VDec = GameObject.Find("VDec");
    movables.Add(VDec.GetComponent<Grabbable>());

  }

  void TryAct(GameObject actable, float twist_val, ref float r_twist)
  {
         if(actable == TInc) thermo.inc_t();
    else if(actable == TDec) thermo.dec_t();
    else if(actable == PInc) thermo.inc_p();
    else if(actable == PDec) thermo.dec_p();
    else if(actable == VInc) thermo.inc_v();
    else if(actable == VDec) thermo.dec_v();
    else
    {
      Dial d = actable.GetComponent<Dial>();
      if(d != null)
      {
        float dtwist = r_twist-twist_val;
        d.val = Mathf.Clamp(d.val-dtwist,0.0f,1.0f);
        DEBUGTEXTS[0].text = twist_val.ToString();
        DEBUGTEXTS[1].text = dtwist.ToString();
        DEBUGTEXTS[2].text = d.val.ToString();
      }
    }

    r_twist = twist_val;
  }

  void TryGrab(bool which, float htrigger_val, float itrigger_val, float twist_val, ref bool r_htrigger, ref bool r_itrigger, ref float r_twist, ref GameObject r_hand, ref GameObject r_grabbed, ref GameObject r_ohand, ref GameObject r_ograbbed)
  {
    float htrigger_threshhold = 0.1f;
    float itrigger_threshhold = 0.1f;

    int htrigger_delta = 0;
    if(!r_htrigger && htrigger_val > htrigger_threshhold)
    {
      htrigger_delta = 1;
      r_htrigger = true;
    }
    else if(r_htrigger && htrigger_val <= htrigger_threshhold)
    {
      htrigger_delta = -1;
      r_htrigger = false;
    }

    int itrigger_delta = 0;
    if(!r_itrigger && itrigger_val > itrigger_threshhold)
    {
      itrigger_delta = 1;
      r_itrigger = true;
      r_twist = twist_val;
    }
    else if(r_itrigger && itrigger_val <= itrigger_threshhold)
    {
      itrigger_delta = -1;
      r_itrigger = false;
    }

    if(r_grabbed == null && htrigger_delta == 1)
    {
      for(int i = 0; r_grabbed == null && i < movables.Count; i++)
      {
        if(
           ( which && movables[i].lintersect) ||
           (!which && movables[i].rintersect)
          )
        {
          r_grabbed = movables[i].gameObject;
          r_grabbed.transform.SetParent(r_hand.transform);
          r_hand.GetComponent<MeshRenderer>().material = hand_grabbing;
          if(r_grabbed == r_ograbbed)
          {
            r_ograbbed = null;
            r_ohand.GetComponent<MeshRenderer>().material = hand_intersecting;
          }
        }
      }
      for(int i = 0; r_grabbed == null && i < dials.Count; i++)
      {
        if(
           ( which && dials[i].lintersect) ||
           (!which && dials[i].rintersect)
          )
        {
          r_grabbed = dials[i].gameObject;
          r_hand.GetComponent<MeshRenderer>().material = hand_grabbing;
          if(r_grabbed == r_ograbbed)
          {
            r_ograbbed = null;
            r_ohand.GetComponent<MeshRenderer>().material = hand_intersecting;
          }
        }
      }
    }
    else if(r_grabbed && htrigger_delta == -1)
    {
      Tool t = r_grabbed.GetComponent<Tool>();
      if(t)
      {
        if(
           ( which && t.active.GetComponent<Grabbable>().lintersect) ||
           (!which && t.active.GetComponent<Grabbable>().rintersect)
          )
        {
          r_grabbed.transform.SetParent(t.active.transform);
          r_grabbed.transform.localPosition = new Vector3(0f,0f,0f);
          r_grabbed.transform.localRotation = Quaternion.identity;
        }
        else if(
           ( which && t.storage.GetComponent<Grabbable>().lintersect) ||
           (!which && t.storage.GetComponent<Grabbable>().rintersect)
          )
        {
          r_grabbed.transform.SetParent(t.storage.transform);
          r_grabbed.transform.localPosition = new Vector3(0f,0f,0f);
          r_grabbed.transform.localRotation = Quaternion.identity;
        }
        else t = null;
      }
      if(t == null) r_grabbed.transform.SetParent(r_grabbed.GetComponent<Grabbable>().og_parent); //ok to do, even with a dial

      r_grabbed = null;
    }
    if(r_grabbed == null)
    {
      r_hand.GetComponent<MeshRenderer>().material = hand_empty;
      for(int i = 0; i < movables.Count; i++)
      {
        if(
           ( which && movables[i].lintersect) ||
           (!which && movables[i].rintersect)
          )
        {
          r_hand.GetComponent<MeshRenderer>().material = hand_intersecting;
          break;
        }
      }
    }

    if(r_grabbed && r_itrigger) TryAct(r_grabbed, twist_val, ref r_twist);
  }

  // Update is called once per frame
  void Update()
  {
    TryGrab(true,  OVRInput.Get(OVRInput.Axis1D.PrimaryHandTrigger),   OVRInput.Get(OVRInput.Axis1D.PrimaryIndexTrigger),   lhand.transform.rotation.z, ref lhtrigger, ref litrigger, ref ltwist, ref lhand, ref lgrabbed, ref rhand, ref rgrabbed); //left hand
    TryGrab(false, OVRInput.Get(OVRInput.Axis1D.SecondaryHandTrigger), OVRInput.Get(OVRInput.Axis1D.SecondaryIndexTrigger), rhand.transform.rotation.z, ref rhtrigger, ref ritrigger, ref rtwist, ref rhand, ref rgrabbed, ref lhand, ref lgrabbed); //right hand
  }

}

