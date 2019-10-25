using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class World : MonoBehaviour
{

  ThermoMath thermo;
  GameObject lhand;
  GameObject rhand;

  public Material hand_empty;
  public Material hand_intersecting;
  public Material hand_grabbing;

  List<Touchable> grabbables;
  GameObject lgrabbed = null;
  GameObject rgrabbed = null;

  List<Tool> tools;

  bool ltrigger = false;
  bool rtrigger = false;

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
    thermo = GameObject.Find("Oracle").GetComponent<ThermoMath>();

    lhand  = GameObject.Find("LHand");
    rhand  = GameObject.Find("RHand");

    lhand.GetComponent<MeshRenderer>().material = hand_empty;
    rhand.GetComponent<MeshRenderer>().material = hand_empty;

    tools = new List<Tool>();
    tools.Add(GameObject.Find("Tool_Insulator").GetComponent<Tool>());

    grabbables = new List<Touchable>();
    for(int i = 0; i < tools.Count; i++) grabbables.Add(tools[i].gameObject.GetComponent<Touchable>()); //important that tools take priority, so they can be grabbed and removed
    grabbables.Add(GameObject.Find("Graph").GetComponent<Touchable>());
    grabbables.Add(GameObject.Find("Vessel").GetComponent<Touchable>());

    //buttons
    TInc = GameObject.Find("TInc");
    grabbables.Add(TInc.GetComponent<Touchable>());
    TDec = GameObject.Find("TDec");
    grabbables.Add(TDec.GetComponent<Touchable>());
    PInc = GameObject.Find("PInc");
    grabbables.Add(PInc.GetComponent<Touchable>());
    PDec = GameObject.Find("PDec");
    grabbables.Add(PDec.GetComponent<Touchable>());
    VInc = GameObject.Find("VInc");
    grabbables.Add(VInc.GetComponent<Touchable>());
    VDec = GameObject.Find("VDec");
    grabbables.Add(VDec.GetComponent<Touchable>());

  }

  void TryAct(GameObject actable)
  {
         if(actable == TInc) thermo.inc_t();
    else if(actable == TDec) thermo.dec_t();
    else if(actable == PInc) thermo.inc_p();
    else if(actable == PDec) thermo.dec_p();
    else if(actable == VInc) thermo.inc_v();
    else if(actable == VDec) thermo.dec_v();
  }

  // Update is called once per frame
  void Update()
  {
    float trigger_threshhold = 0.1f;

    //left hand
    int ltrigger_delta = 0;
    if(!ltrigger && OVRInput.Get(OVRInput.Axis1D.PrimaryHandTrigger) > trigger_threshhold)
    {
      ltrigger_delta = 1;
      ltrigger = true;
    }
    else if(ltrigger && OVRInput.Get(OVRInput.Axis1D.PrimaryHandTrigger) <= trigger_threshhold)
    {
      ltrigger_delta = -1;
      ltrigger = false;
    }

    if(lgrabbed == null && ltrigger_delta == 1)
    {
      for(int i = 0; lgrabbed == null && i < grabbables.Count; i++)
      {
        if(grabbables[i].lintersect)
        {
          lgrabbed = grabbables[i].gameObject;
          lgrabbed.transform.SetParent(lhand.transform);
          lhand.GetComponent<MeshRenderer>().material = hand_grabbing;
          if(lgrabbed == rgrabbed)
          {
            rgrabbed = null;
            rhand.GetComponent<MeshRenderer>().material = hand_intersecting;
          }
        }
      }
    }
    else if(lgrabbed && ltrigger_delta == -1)
    {
      Tool t = lgrabbed.GetComponent<Tool>();
      if(t)
      {
        if(t.active.GetComponent<Touchable>().lintersect)
        {
          lgrabbed.transform.SetParent(t.active.transform);
          lgrabbed.transform.localPosition = new Vector3(0f,0f,0f);
          lgrabbed.transform.localRotation = new Vector3(0f,0f,0f);
        }
        else if(t.storage.GetComponent<Touchable>().lintersect)
        {
          lgrabbed.transform.SetParent(t.storage.transform);
          lgrabbed.transform.localPosition = new Vector3(0f,0f,0f);
          lgrabbed.transform.localRotation = new Vector3(0f,0f,0f);
        }
        else t = null;
      }
      if(t == null) lgrabbed.transform.SetParent(lgrabbed.GetComponent<Touchable>().og_parent);

      lgrabbed = null;
    }
    if(lgrabbed == null)
    {
      lhand.GetComponent<MeshRenderer>().material = hand_empty;
      for(int i = 0; i < grabbables.Count; i++)
      {
        if(grabbables[i].lintersect)
        {
          lhand.GetComponent<MeshRenderer>().material = hand_intersecting;
          break;
        }
      }
    }

    //right hand
    int rtrigger_delta = 0;
    if(!rtrigger && OVRInput.Get(OVRInput.Axis1D.SecondaryHandTrigger) > trigger_threshhold)
    {
      rtrigger_delta = 1;
      rtrigger = true;
    }
    else if(rtrigger && OVRInput.Get(OVRInput.Axis1D.SecondaryHandTrigger) <= trigger_threshhold)
    {
      rtrigger_delta = -1;
      rtrigger = false;
    }
    if(rgrabbed == null && rtrigger_delta == 1)
    {
      for(int i = 0; rgrabbed == null && i < grabbables.Count; i++)
      {
        if(grabbables[i].rintersect)
        {
          rgrabbed = grabbables[i].gameObject;
          rgrabbed.transform.SetParent(rhand.transform);
          rhand.GetComponent<MeshRenderer>().material = hand_grabbing;
          if(rgrabbed == lgrabbed)
          {
            lgrabbed = null;
            lhand.GetComponent<MeshRenderer>().material = hand_intersecting;
          }
        }
      }
    }
    else if(rgrabbed && rtrigger_delta == -1)
    {
      Tool t = rgrabbed.GetComponent<Tool>();
      if(t)
      {
        if(t.active.GetComponent<Touchable>().rintersect)
        {
          rgrabbed.transform.SetParent(t.active.transform);
          rgrabbed.transform.localPosition = new Vector3(0f,0f,0f);
          rgrabbed.transform.localRotation = new Vector3(0f,0f,0f);
        }
        else if(t.storage.GetComponent<Touchable>().rintersect)
        {
          rgrabbed.transform.SetParent(t.storage.transform);
          rgrabbed.transform.localPosition = new Vector3(0f,0f,0f);
          rgrabbed.transform.localRotation = new Vector3(0f,0f,0f);
        }
        else t = null;
      }
      if(t == null) rgrabbed.transform.SetParent(rgrabbed.GetComponent<Touchable>().og_parent);

      rgrabbed = null;
    }
    if(rgrabbed == null)
    {
      rhand.GetComponent<MeshRenderer>().material = hand_empty;
      for(int i = 0; i < grabbables.Count; i++)
      {
        if(grabbables[i].rintersect)
        {
          rhand.GetComponent<MeshRenderer>().material = hand_intersecting;
          break;
        }
      }
    }

    //buttons
    if(lgrabbed && OVRInput.Get(OVRInput.Axis1D.PrimaryIndexTrigger)   > trigger_threshhold) TryAct(lgrabbed);
    if(rgrabbed && OVRInput.Get(OVRInput.Axis1D.SecondaryIndexTrigger) > trigger_threshhold) TryAct(rgrabbed);

  }

}

