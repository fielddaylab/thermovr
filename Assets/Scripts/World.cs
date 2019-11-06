using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class World : MonoBehaviour
{
  List<TextMesh> DEBUGTEXTS;

  ThermoMath thermo;
  GameObject lhand;
  Vector3 lhand_pos;
  Vector3 lhand_vel;
  MeshRenderer lhand_meshrenderer;
  GameObject rhand;
  Vector3 rhand_pos;
  Vector3 rhand_vel;
  MeshRenderer rhand_meshrenderer;

  public Material hand_empty;
  public Material hand_intersecting;
  public Material hand_grabbing;

  List<Grabbable> movables;
  GameObject workspace;
  GameObject handle_workspace;
  GameObject lgrabbed = null;
  GameObject rgrabbed = null;
  float lz = 0.0f;
  float rz = 0.0f;
  float ly = 0.0f;
  float ry = 0.0f;

  List<Tool> tools;
  Tool tool_insulator;
  Tool tool_clamp;
  Tool tool_burner;
  Tool tool_coil;
  Tool tool_weight;
  Tool tool_balloon;

  ParticleSystem flame; //special case

  GameObject vessel;
  GameObject graph;

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

    lhand  = GameObject.Find("LeftControllerAnchor");
    lhand_pos = lhand.transform.position;
    lhand_vel = new Vector3(0.0f,0.0f,0.0f);
    lhand_meshrenderer = lhand.transform.GetChild(0).gameObject.GetComponent<MeshRenderer>();
    rhand  = GameObject.Find("RightControllerAnchor");
    rhand_pos = rhand.transform.position;
    rhand_vel = new Vector3(0.0f,0.0f,0.0f);
    rhand_meshrenderer = rhand.transform.GetChild(0).gameObject.GetComponent<MeshRenderer>();

    lhand_meshrenderer.material = hand_empty;
    rhand_meshrenderer.material = hand_empty;

    Tool t;
    tools = new List<Tool>();
    t = GameObject.Find("Tool_Insulator").GetComponent<Tool>(); tool_insulator = t; tools.Add(t);
    t = GameObject.Find("Tool_Clamp"    ).GetComponent<Tool>(); tool_clamp     = t; tools.Add(t);
    t = GameObject.Find("Tool_Burner"   ).GetComponent<Tool>(); tool_burner    = t; tools.Add(t);
    t = GameObject.Find("Tool_Coil"     ).GetComponent<Tool>(); tool_coil      = t; tools.Add(t);
    t = GameObject.Find("Tool_Weight"   ).GetComponent<Tool>(); tool_weight    = t; tools.Add(t);
    t = GameObject.Find("Tool_Balloon"  ).GetComponent<Tool>(); tool_balloon   = t; tools.Add(t);

    flame = GameObject.Find("Flame").GetComponent<ParticleSystem>();

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
      t.stored = true;
      g.transform.localPosition = new Vector3(0.0f,0.0f,0.0f);
    }

    vessel = GameObject.Find("Vessel");
    graph = GameObject.Find("Graph");

    movables = new List<Grabbable>();
    for(int i = 0; i < tools.Count; i++) movables.Add(tools[i].grabbable); //important that tools take priority, so they can be grabbed and removed
    movables.Add(graph.GetComponent<Grabbable>());
    movables.Add(vessel.GetComponent<Grabbable>());
  }

  void TryApplyTool(Tool t)
  {
    if(t == tool_insulator)
    {
      //visual
      //math
      if(!t.engaged) return;
      if(tool_clamp.engaged) thermo.tp_get_p(t.dial_dial.val);
      else                   thermo.tp_get_v(t.dial_dial.val);
    }
    else if(t == tool_clamp)
    {
      //visual change handled by thermomath (alters piston height)
      //math
      if(!t.engaged) return;
      if(tool_insulator.engaged) thermo.vp_get_p(t.dial_dial.val);
      else                       thermo.vp_get_t(t.dial_dial.val);
    }
    else if(t == tool_weight || t == tool_balloon)
    {
      //visual
      float v = 1.0f;
           if(t == tool_weight)  v += tool_weight.dial_dial.val;
      else if(t == tool_balloon) v += tool_balloon.dial_dial.val;
      Vector3 scale = new Vector3(v,v,v);
      t.gameObject.transform.localScale = scale;
      if(t == tool_balloon) v *= 0.125f;
      scale = new Vector3(v,v,v);
      t.active_available.transform.localScale = scale;
      t.active_snap.transform.localScale = scale;
      t.storage_available.transform.localScale = scale;
      t.storage_snap.transform.localScale = scale;

      //math
      if(!t.engaged || tool_clamp.engaged) return; //weight does nothing!
      float weight = 0;
      if(tool_weight.engaged)  weight += tool_weight.dial_dial.val;
      if(tool_balloon.engaged) weight -= tool_balloon.dial_dial.val;
      weight = (weight+1)/2; //normalize 0-1

      if(tool_insulator.engaged) thermo.pp_get_v(weight);
      else                       thermo.pp_get_t(weight);
    }
    else if(t == tool_burner)
    {
      //visual
      var vel = flame.velocityOverLifetime;
      vel.speedModifierMultiplier = Mathf.Lerp(0.1f,0.5f,t.dial_dial.val);
      //math change happens passively- only need to update visuals
    }
    else if(t == tool_coil)
    {
      //visual
      ;
      //math change happens passively- only need to update visuals
    }

  }

  void TryAct(GameObject actable, float z_val, float y_val, ref float r_z, ref float r_y)
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
        float dx = (r_z-z_val)*2.0f;
        d.val = Mathf.Clamp(d.val-dx,0.0f,1.0f);

        TryApplyTool(t);
      }
    }
  }

  void TryGrab(bool which, float htrigger_val, float itrigger_val, float z_val, float y_val, Vector3 hand_vel, ref bool r_htrigger, ref bool r_itrigger, ref float r_z, ref float r_y, ref GameObject r_hand, ref GameObject r_grabbed, ref GameObject r_ohand, ref GameObject r_ograbbed)
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
      r_z = z_val;
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
          movables[i].grabbed = true;
          Tool t = r_grabbed.GetComponent<Tool>();
          if(t)
          {
            t.engaged = false;
            t.stored = false;
            t.rigidbody.isKinematic = true;
            t.boxcollider.isTrigger = false;
          }
          VisAid v = r_grabbed.GetComponent<VisAid>();
          if(v)
          {
            v.stored = false;
            v.rigidbody.isKinematic = true;
          }
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
            tools[i].dial_grabbable.grabbed = true;
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
          g.grabbed = true;
          if(r_grabbed == r_ograbbed) r_ograbbed = null;
        }
      }
    }
    else if(r_grabbed && htrigger_delta == -1)
    {
      Tool t = r_grabbed.GetComponent<Tool>();
      if(t)
      {
        if(t.active_ghost.tintersect)
        {
          r_grabbed.transform.SetParent(t.active.transform);
          t.engaged = true;
          t.boxcollider.isTrigger = true;
          if(
            (t == tool_insulator && tool_clamp.engaged) ||
            (t == tool_clamp && tool_insulator.engaged)
          )
          {
            Tool ot;
            if(t == tool_insulator) ot = tool_clamp;
            else                    ot = tool_insulator;
            ot.engaged = false;
            ot.boxcollider.isTrigger = false;
            ot.gameObject.transform.SetParent(ot.grabbable.og_parent);
            ot.rigidbody.isKinematic = false;
            ot.rigidbody.velocity = new Vector3(0f,0f,0f);
          }
          TryApplyTool(t);
          r_grabbed.transform.localPosition = new Vector3(0f,0f,0f);
          r_grabbed.transform.localRotation = Quaternion.identity;
        }
        else if(t.storage_ghost.tintersect)
        {
          r_grabbed.transform.SetParent(t.storage.transform);
          t.stored = true;
          r_grabbed.transform.localPosition = new Vector3(0f,0f,0f);
          r_grabbed.transform.localRotation = Quaternion.identity;
        }
        else
        {
          r_grabbed.transform.SetParent(r_grabbed.GetComponent<Grabbable>().og_parent);
          t.rigidbody.isKinematic = false;
          t.rigidbody.velocity = hand_vel;
        }
      }
      else
      {
        r_grabbed.transform.SetParent(r_grabbed.GetComponent<Grabbable>().og_parent); //ok to do, even with a dial
        VisAid v = r_grabbed.GetComponent<VisAid>();
        if(v)
        {
          v.rigidbody.isKinematic = false;
          v.rigidbody.velocity = hand_vel;
        }
      }

      r_grabbed.GetComponent<Grabbable>().grabbed = false;
      r_grabbed = null;
    }

    if(r_grabbed) TryAct(r_grabbed, z_val, y_val, ref r_z, ref r_y);

    r_z = z_val;
    r_y = y_val;
  }

  void UpdateGrabVis()
  {
    for(int i = 0; i < tools.Count; i++)
    {
      Tool t = tools[i];
      GameObject g = t.gameObject;

      if(lgrabbed == g || rgrabbed == g)
      {
        //active
        if(t.active_ghost.tintersect)
        {
          t.active_available_meshrenderer.enabled = false;
          t.active_snap_meshrenderer.enabled =      true;
        }
        else
        {
          t.active_available_meshrenderer.enabled = true;
          t.active_snap_meshrenderer.enabled =      false;
        }
        //storage
        if(t.storage_ghost.tintersect)
        {
          t.storage_available_meshrenderer.enabled = false;
          t.storage_snap_meshrenderer.enabled =      true;
        }
        else
        {
          t.storage_available_meshrenderer.enabled = true;
          t.storage_snap_meshrenderer.enabled =      false;
        }
      }
      else
      {
        t.active_snap_meshrenderer.enabled =      false;
        t.active_available_meshrenderer.enabled = false;
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
    for(int i = 0; i < tools.Count; i++)
    {
      gr = tools[i].dial_grabbable;
      if(gr.lintersect) lintersect = true;
      if(gr.rintersect) rintersect = true;
    }
    gr = handle_workspace.GetComponent<Grabbable>();
    if(gr.lintersect) lintersect = true;
    if(gr.rintersect) rintersect = true;

         if(lgrabbed)   lhand_meshrenderer.material = hand_grabbing;
    else if(lintersect) lhand_meshrenderer.material = hand_intersecting;
    else                lhand_meshrenderer.material = hand_empty;

         if(rgrabbed)   rhand_meshrenderer.material = hand_grabbing;
    else if(rintersect) rhand_meshrenderer.material = hand_intersecting;
    else                rhand_meshrenderer.material = hand_empty;
  }

  // Update is called once per frame
  void Update()
  {
    //running blended average
    lhand_vel += (lhand.transform.position-lhand_pos)/Time.deltaTime;
    lhand_vel *= 0.5f;
    lhand_pos = lhand.transform.position;

    rhand_vel += (rhand.transform.position-rhand_pos)/Time.deltaTime;
    rhand_vel *= 0.5f;
    rhand_pos = rhand.transform.position;

    /*
    //snap from delay
    lhand_pos = Vector3.Lerp(lhand_pos,lhand.transform.position,0.01f);
    lhand_vel = (lhand.transform.position-lhand_pos)/Time.deltaTime;

    rhand_pos = Vector3.Lerp(rhand_pos,rhand.transform.position,0.01f);
    rhand_vel = (rhand.transform.position-rhand_pos)/Time.deltaTime;
    */

    //passive effects
    float d_heat = 0;
    if(tool_burner.engaged) d_heat += tool_burner.dial_dial.val;
    if(tool_coil.engaged)   d_heat -= tool_coil.dial_dial.val;
    if(d_heat != 0.0f)
    {
      thermo.h_get_t(d_heat);
    }

    float lhandt  = OVRInput.Get(OVRInput.Axis1D.PrimaryHandTrigger);
    float lindext = OVRInput.Get(OVRInput.Axis1D.PrimaryIndexTrigger);
    float rhandt  = OVRInput.Get(OVRInput.Axis1D.SecondaryHandTrigger);
    float rindext = OVRInput.Get(OVRInput.Axis1D.SecondaryIndexTrigger);
    //index compatibility
    lhandt += lindext;
    rhandt += rindext;
    TryGrab(true,  lhandt, lindext, lhand.transform.position.z, lhand.transform.position.y, lhand_vel, ref lhtrigger, ref litrigger, ref lz, ref ly, ref lhand, ref lgrabbed, ref rhand, ref rgrabbed); //left hand
    TryGrab(false, rhandt, rindext, rhand.transform.position.z, rhand.transform.position.y, rhand_vel, ref rhtrigger, ref ritrigger, ref rz, ref ry, ref rhand, ref rgrabbed, ref lhand, ref lgrabbed); //right hand

    UpdateGrabVis();

//    DEBUGTEXTS[0].text = OVRInput.Get(OVRInput.Axis1D.PrimaryIndexTrigger).ToString();
//    DEBUGTEXTS[1].text = OVRInput.Get(OVRInput.Axis1D.SecondaryIndexTrigger).ToString();
//    DEBUGTEXTS[2].text = "NA";
  }

}

