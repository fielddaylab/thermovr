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
  float ly = 0.0f;
  float ry = 0.0f;

  List<Tool> tools;
  List<Grabbable> dials;
  GameObject dial_insulator;
  GameObject dial_clamp;
  GameObject dial_burner;
  GameObject dial_coil;
  GameObject dial_weight;
  GameObject dial_balloon;

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

    Tool t;
    tools = new List<Tool>();
    t = GameObject.Find("Tool_Insulator").GetComponent<Tool>(); dial_insulator = t.dial; tools.Add(t);
    t = GameObject.Find("Tool_Clamp"    ).GetComponent<Tool>(); dial_clamp     = t.dial; tools.Add(t);
    t = GameObject.Find("Tool_Burner"   ).GetComponent<Tool>(); dial_burner    = t.dial; tools.Add(t);
    t = GameObject.Find("Tool_Coil"     ).GetComponent<Tool>(); dial_coil      = t.dial; tools.Add(t);
    t = GameObject.Find("Tool_Weight"   ).GetComponent<Tool>(); dial_weight    = t.dial; tools.Add(t);
    t = GameObject.Find("Tool_Balloon"  ).GetComponent<Tool>(); dial_balloon   = t.dial; tools.Add(t);

    for(int i = 0; i < tools.Count; i++)
    {
      t = tools[i];
      GameObject g = t.gameObject;
      g.transform.SetParent(t.storage.gameObject.transform);
      g.transform.localPosition = new Vector3(0.0f,0.0f,0.0f);
    }

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

  void TryAct(GameObject actable, float y_val, ref float r_y)
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
        float dy = (r_y-y_val);
        d.val = Mathf.Clamp(d.val-dy,0.0f,1.0f);

        if(actable == dial_insulator) thermo.set_tp(d.val);
        if(actable == dial_clamp)     thermo.set_vp(d.val);
        if(actable == dial_burner)    ; //do nothing; passive
        if(actable == dial_coil)      ; //do nothing; passive
        if(actable == dial_weight)    thermo.set_pp(d.val);
        if(actable == dial_balloon)   thermo.set_pp(d.val);
      }
    }

    r_y = y_val;
  }

  void TryGrab(bool which, float htrigger_val, float itrigger_val, float y_val, ref bool r_htrigger, ref bool r_itrigger, ref float r_y, ref GameObject r_hand, ref GameObject r_grabbed, ref GameObject r_ohand, ref GameObject r_ograbbed)
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
      r_y = y_val;
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

    if(r_grabbed && r_itrigger) TryAct(r_grabbed, y_val, ref r_y);
  }

  // Update is called once per frame
  void Update()
  {
    TryGrab(true,  OVRInput.Get(OVRInput.Axis1D.PrimaryHandTrigger),   OVRInput.Get(OVRInput.Axis1D.PrimaryIndexTrigger),   lhand.transform.position.y, ref lhtrigger, ref litrigger, ref ly, ref lhand, ref lgrabbed, ref rhand, ref rgrabbed); //left hand
    TryGrab(false, OVRInput.Get(OVRInput.Axis1D.SecondaryHandTrigger), OVRInput.Get(OVRInput.Axis1D.SecondaryIndexTrigger), rhand.transform.position.y, ref rhtrigger, ref ritrigger, ref ry, ref rhand, ref rgrabbed, ref lhand, ref lgrabbed); //right hand
    /*
    DEBUGTEXTS[0].text = lhand.transform.eulerAngles.x.ToString();
    DEBUGTEXTS[1].text = lhand.transform.eulerAngles.y.ToString();
    DEBUGTEXTS[2].text = lhand.transform.eulerAngles.z.ToString();
    */
  }

}

