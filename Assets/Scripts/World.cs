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
  GameObject workspace;
  GameObject handle_workspace;
  GameObject lgrabbed = null;
  GameObject rgrabbed = null;
  float ly = 0.0f;
  float ry = 0.0f;

  List<Tool> tools;
  Tool tool_insulator;
  Tool tool_clamp;
  Tool tool_burner;
  Tool tool_coil;
  Tool tool_weight;
  Tool tool_balloon;

  bool lhtrigger = false;
  bool rhtrigger = false;
  bool litrigger = false;
  bool ritrigger = false;

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
    t = GameObject.Find("Tool_Insulator").GetComponent<Tool>(); tool_insulator = t; tools.Add(t);
    t = GameObject.Find("Tool_Clamp"    ).GetComponent<Tool>(); tool_clamp     = t; tools.Add(t);
    t = GameObject.Find("Tool_Burner"   ).GetComponent<Tool>(); tool_burner    = t; tools.Add(t);
    t = GameObject.Find("Tool_Coil"     ).GetComponent<Tool>(); tool_coil      = t; tools.Add(t);
    t = GameObject.Find("Tool_Weight"   ).GetComponent<Tool>(); tool_weight    = t; tools.Add(t);
    t = GameObject.Find("Tool_Balloon"  ).GetComponent<Tool>(); tool_balloon   = t; tools.Add(t);

    workspace = GameObject.Find("Workspace");
    handle_workspace = GameObject.Find("Handle_Workspace");

    for(int i = 0; i < tools.Count; i++)
    {
      t = tools[i];
      t.active_available_meshrenderer.enabled = false;
      t.active_snap_meshrenderer.enabled = false;
      t.storage_available_meshrenderer.enabled = false;
      t.storage_snap_meshrenderer.enabled = false;
      GameObject g = t.gameObject;
      g.transform.SetParent(t.storage.gameObject.transform);
      g.transform.localPosition = new Vector3(0.0f,0.0f,0.0f);
    }

    movables = new List<Grabbable>();
    for(int i = 0; i < tools.Count; i++) movables.Add(tools[i].grabbable); //important that tools take priority, so they can be grabbed and removed
    movables.Add(GameObject.Find("Graph").GetComponent<Grabbable>());
    movables.Add(GameObject.Find("Vessel").GetComponent<Grabbable>());
  }

  void TryApplyTool(Tool t)
  {
    if(!t.engaged) return;
    if(t == tool_insulator)
    {
      if(tool_clamp.engaged) thermo.tp_get_p(t.dial_dial.val);
      else                   thermo.tp_get_v(t.dial_dial.val);
    }
    else if(t == tool_clamp)
    {
      if(tool_insulator.engaged) thermo.vp_get_p(t.dial_dial.val);
      else                       thermo.vp_get_t(t.dial_dial.val);
    }
    else if(t == tool_weight || t == tool_balloon)
    {
      if(tool_clamp.engaged) return; //weight does nothing!
      float weight = 0;
      if(tool_weight.engaged)  weight += tool_weight.dial_dial.val;
      if(tool_balloon.engaged) weight -= tool_balloon.dial_dial.val;
      weight = (weight+1)/2; //normalize 0-1

      if(tool_insulator.engaged) thermo.pp_get_v(weight);
      else                       thermo.pp_get_t(weight);
    }
    else if(t == tool_burner) ; //do nothing; passive
    else if(t == tool_coil)   ; //do nothing; passive
  }

  void TryAct(GameObject actable, float y_val, ref float r_y)
  {
    if(actable == handle_workspace)
    {
      float dy = (r_y-y_val);
      workspace.transform.position = new Vector3(workspace.transform.position.x,workspace.transform.position.y-dy,workspace.transform.position.z);
    }
    else
    {
      Dial d = actable.GetComponent<Dial>();
      if(d != null)
      {
        Tool t = d.tool.GetComponent<Tool>();
        float dy = (r_y-y_val);
        d.val = Mathf.Clamp(d.val-dy,0.0f,1.0f);

        TryApplyTool(t);
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
      //first try movables
      for(int i = 0; r_grabbed == null && i < movables.Count; i++)
      {
        if(
           ( which && movables[i].lintersect) ||
           (!which && movables[i].rintersect)
          )
        {
          r_grabbed = movables[i].gameObject;
          Tool t = r_grabbed.GetComponent<Tool>();
          if(t) t.engaged = false;
          r_grabbed.transform.SetParent(r_hand.transform);
          if(r_grabbed == r_ograbbed) r_ograbbed = null;
        }
      }
      //then dials
      if(r_grabbed == null)
      {
        for(int i = 0; i < tools.Count; i++)
        {
          if(
             ( which && tools[i].dial_grabbable.lintersect) ||
             (!which && tools[i].dial_grabbable.rintersect)
            )
          {
            r_grabbed = tools[i].dial;
            if(r_grabbed == r_ograbbed) r_ograbbed = null;
          }
        }
      }
      //then extraaneous
      if(r_grabbed == null)
      {
        Grabbable g = handle_workspace.GetComponent<Grabbable>();
        if(
          ( which && g.lintersect) ||
          (!which && g.rintersect)
        )
        {
          r_grabbed = handle_workspace;
          if(r_grabbed == r_ograbbed) r_ograbbed = null;
        }
      }
    }
    else if(r_grabbed && htrigger_delta == -1)
    {
      Tool t = r_grabbed.GetComponent<Tool>();
      if(t)
      {
        if(
           ( which && t.active_grabbable.lintersect) ||
           (!which && t.active_grabbable.rintersect)
          )
        {
          r_grabbed.transform.SetParent(t.active.transform);
          t.engaged = true;
          TryApplyTool(t);
          r_grabbed.transform.localPosition = new Vector3(0f,0f,0f);
          r_grabbed.transform.localRotation = Quaternion.identity;
        }
        else if(
           ( which && t.storage_grabbable.lintersect) ||
           (!which && t.storage_grabbable.rintersect)
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

    if(r_grabbed && r_itrigger) TryAct(r_grabbed, y_val, ref r_y);
  }

  void UpdateGrabVis()
  {
    for(int i = 0; i < tools.Count; i++)
    {
      Tool t = tools[i];
      GameObject g = t.gameObject;

      //active
      if(
        (lgrabbed == g && t.active_grabbable.lintersect) ||
        (rgrabbed == g && t.active_grabbable.rintersect)
      )
      {
        t.active_available_meshrenderer.enabled = false;
        t.active_snap_meshrenderer.enabled =      true;
      }
      else if(lgrabbed == g || rgrabbed == g)
      {
        t.active_available_meshrenderer.enabled = true;
        t.active_snap_meshrenderer.enabled =      false;
      }
      else
      {
        t.active_snap_meshrenderer.enabled =      false;
        t.active_available_meshrenderer.enabled = false;
      }
      //storage
      if(
        (lgrabbed == g && t.storage_grabbable.lintersect) ||
        (rgrabbed == g && t.storage_grabbable.rintersect)
      )
      {
        t.storage_available_meshrenderer.enabled = false;
        t.storage_snap_meshrenderer.enabled =      true;
      }
      else if(lgrabbed == g || rgrabbed == g)
      {
        t.storage_available_meshrenderer.enabled = true;
        t.storage_snap_meshrenderer.enabled =      false;
      }
      else
      {
        t.storage_snap_meshrenderer.enabled =      false;
        t.storage_available_meshrenderer.enabled = false;
      }
    }

    Grabbable gr;
    bool lintersect = false;
    bool rintersect = false;
    for(int i = 0; i < movables.Count; i++)
    {
      gr = movables[i];
      if(gr.lintersect) lintersect = true;
      if(gr.rintersect) rintersect = true;
    }
    gr = handle_workspace.GetComponent<Grabbable>();
    if(gr.lintersect) lintersect = true;
    if(gr.rintersect) rintersect = true;

         if(lgrabbed)   lhand.GetComponent<MeshRenderer>().material = hand_grabbing;
    else if(lintersect) lhand.GetComponent<MeshRenderer>().material = hand_intersecting;
    else                lhand.GetComponent<MeshRenderer>().material = hand_empty;

         if(rgrabbed)   rhand.GetComponent<MeshRenderer>().material = hand_grabbing;
    else if(rintersect) rhand.GetComponent<MeshRenderer>().material = hand_intersecting;
    else                rhand.GetComponent<MeshRenderer>().material = hand_empty;
  }

  // Update is called once per frame
  void Update()
  {
    //passive effects
    float d_heat = 0;
    if(tool_burner.engaged) d_heat += tool_burner.dial_dial.val;
    if(tool_coil.engaged)   d_heat -= tool_coil.dial_dial.val;
    if(d_heat != 0.0f)
    {
      thermo.h_get_t(d_heat);
    }

    TryGrab(true,  OVRInput.Get(OVRInput.Axis1D.PrimaryHandTrigger),   OVRInput.Get(OVRInput.Axis1D.PrimaryIndexTrigger),   lhand.transform.position.y, ref lhtrigger, ref litrigger, ref ly, ref lhand, ref lgrabbed, ref rhand, ref rgrabbed); //left hand
    TryGrab(false, OVRInput.Get(OVRInput.Axis1D.SecondaryHandTrigger), OVRInput.Get(OVRInput.Axis1D.SecondaryIndexTrigger), rhand.transform.position.y, ref rhtrigger, ref ritrigger, ref ry, ref rhand, ref rgrabbed, ref lhand, ref lgrabbed); //right hand

    UpdateGrabVis();
    /*
    DEBUGTEXTS[0].text = lhand.transform.eulerAngles.x.ToString();
    DEBUGTEXTS[1].text = lhand.transform.eulerAngles.y.ToString();
    DEBUGTEXTS[2].text = lhand.transform.eulerAngles.z.ToString();
    */
  }

}

