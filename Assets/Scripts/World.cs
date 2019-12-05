using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using TMPro;

public class World : MonoBehaviour
{
  public Material hand_empty;
  public Material hand_touching;
  public Material hand_grabbing;
  public Material quiz_default;
  public Material quiz_hi;
  public Material quiz_sel;
  public Material quiz_hisel;

  ThermoState thermo;
  GameObject cam_offset;
  GameObject ceye;
  GameObject lhand;
  GameObject llazer;
  GameObject llazer_viz;
  Vector3 lhand_pos;
  Vector3 lhand_vel;
  MeshRenderer lhand_meshrenderer;
  MeshRenderer llazer_meshrenderer;
  Fadable llazer_fadable;
  Fadable rlazer_fadable;
  GameObject rhand;
  GameObject rlazer;
  GameObject rlazer_viz;
  Vector3 rhand_pos;
  Vector3 rhand_vel;
  MeshRenderer rhand_meshrenderer;
  MeshRenderer rlazer_meshrenderer;

  GameObject vrcenter;
  Lazerable vrcenter_lazerable;
  MeshRenderer vrcenter_backing_meshrenderer;

  List<Touchable> movables;
  GameObject workspace;
  GameObject handle_workspace;
  Touchable handle_workspace_touchable;
  GameObject lgrabbed = null;
  GameObject rgrabbed = null;
  int lhtrigger_delta = 0;
  int litrigger_delta = 0;
  int rhtrigger_delta = 0;
  int ritrigger_delta = 0;
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
  Tool tool_clipboard;

  ParticleSystem flame; //special case

  bool halfed = false;
  List<Halfable> halfables;
  GameObject halfer;
  Touchable halfer_touchable;

  GameObject vessel;
  GameObject graph;

  bool lhtrigger = false;
  bool rhtrigger = false;
  bool litrigger = false;
  bool ritrigger = false;

  int question = 0;
  List<string> questions;
  List<string> options;
  List<int> answers;
  List<int> givens;
  List<Quizo> option_quizos;
  Quizo qconfirm_quizo;
  GameObject qboard;
  Lazerable qboard_lazerable;
  TextMeshPro qtext_tmp;
  int qselected = -1;

  double applied_weight = 0;
  double applied_heat = 0;

  // Start is called before the first frame update
  void Start()
  {
    thermo = GameObject.Find("Oracle").GetComponent<ThermoState>();

    cam_offset = GameObject.Find("CamOffset");
    ceye = GameObject.Find("CenterEyeAnchor");

    lhand  = GameObject.Find("LeftControllerAnchor");
    llazer  = GameObject.Find("LLazer");
    llazer_viz  = llazer.transform.GetChild(0).gameObject;
    lhand_pos = lhand.transform.position;
    lhand_vel = new Vector3(0.0f,0.0f,0.0f);
    lhand_meshrenderer = lhand.transform.GetChild(0).gameObject.GetComponent<MeshRenderer>();
    llazer_meshrenderer = llazer.transform.GetChild(0).gameObject.GetComponent<MeshRenderer>();
    llazer_meshrenderer.enabled = false;
    llazer_fadable = llazer.GetComponent<Fadable>();

    rhand  = GameObject.Find("RightControllerAnchor");
    rlazer  = GameObject.Find("RLazer");
    rlazer_viz  = rlazer.transform.GetChild(0).gameObject;
    rhand_pos = rhand.transform.position;
    rhand_vel = new Vector3(0.0f,0.0f,0.0f);
    rhand_meshrenderer = rhand.transform.GetChild(0).gameObject.GetComponent<MeshRenderer>();
    rlazer_meshrenderer = rlazer.transform.GetChild(0).gameObject.GetComponent<MeshRenderer>();
    rlazer_meshrenderer.enabled = false;
    rlazer_fadable = rlazer.GetComponent<Fadable>();

    llazer_fadable.t_in = 0.1f;
    rlazer_fadable.t_in = 0.1f;
    llazer_fadable.t_start_out = 0.1f;
    rlazer_fadable.t_start_out = 0.1f;
    llazer_fadable.t_end_out = 0.2f;
    llazer_fadable.t_end_out = 0.2f;
    lhand_meshrenderer.material = hand_empty;
    rhand_meshrenderer.material = hand_empty;

    Tool t;
    tools = new List<Tool>();
    t = GameObject.Find("Tool_Insulator").GetComponent<Tool>(); tool_insulator = t; tools.Add(t); t.dial_dial.min_map = 0.0f; t.dial_dial.max_map = 1.0f; t.dial_dial.unit = "n";
    t = GameObject.Find("Tool_Clamp"    ).GetComponent<Tool>(); tool_clamp     = t; tools.Add(t); t.dial_dial.min_map = 0.0f; t.dial_dial.max_map = 1.0f; t.dial_dial.unit = "h";
    t = GameObject.Find("Tool_Burner"   ).GetComponent<Tool>(); tool_burner    = t; tools.Add(t); t.dial_dial.min_map = 0.0f; t.dial_dial.max_map = 1.0f; t.dial_dial.unit = "J/s";
    t = GameObject.Find("Tool_Coil"     ).GetComponent<Tool>(); tool_coil      = t; tools.Add(t); t.dial_dial.min_map = 0.0f; t.dial_dial.max_map = 1.0f; t.dial_dial.unit = "J/s";
    t = GameObject.Find("Tool_Weight"   ).GetComponent<Tool>(); tool_weight    = t; tools.Add(t); t.dial_dial.min_map = 0.0f; t.dial_dial.max_map = 1.0f; t.dial_dial.unit = "kg";
    t = GameObject.Find("Tool_Balloon"  ).GetComponent<Tool>(); tool_balloon   = t; tools.Add(t); t.dial_dial.min_map = 0.0f; t.dial_dial.max_map = 1.0f; t.dial_dial.unit = "kg";
    t = GameObject.Find("Tool_Clipboard").GetComponent<Tool>(); tool_clipboard = t; tools.Add(t); t.dial_dial.min_map = 0.0f; t.dial_dial.max_map = 1.0f; t.dial_dial.unit = "N/A";

    flame = GameObject.Find("Flame").GetComponent<ParticleSystem>();

    workspace = GameObject.Find("Workspace");
    handle_workspace = GameObject.Find("Handle_Workspace");
    handle_workspace_touchable = handle_workspace.GetComponent<Touchable>();

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
      t.textv_tmp.SetText("{0:3}"+t.dial_dial.unit,(float)t.dial_dial.map);
    }

    vessel = GameObject.Find("Vessel");
    graph = GameObject.Find("Graph");

    vrcenter = GameObject.Find("VRCenter");
    vrcenter_lazerable = vrcenter.GetComponent<Lazerable>();
    vrcenter_backing_meshrenderer = GameObject.Find("VRCenter").transform.GetChild(1).GetComponent<MeshRenderer>();

    movables = new List<Touchable>();
    for(int i = 0; i < tools.Count; i++) movables.Add(tools[i].touchable); //important that tools take priority, so they can be grabbed and removed
    movables.Add(graph.GetComponent<Touchable>());
    movables.Add(vessel.GetComponent<Touchable>());

    halfer = GameObject.Find("Halfer");
    halfer_touchable = halfer.GetComponent<Touchable>();
    halfables = new List<Halfable>();
    halfables.Add(GameObject.Find("Container").GetComponent<Halfable>());
    halfables.Add(GameObject.Find("Piston").GetComponent<Halfable>());
    halfables.Add(GameObject.Find("Contents").GetComponent<Halfable>());

    questions = new List<string>();
    options = new List<string>();
    answers = new List<int>();
    givens = new List<int>();

    questions.Add("What's the q?");
    options.Add("A. WHAAA");
    options.Add("B. WHOO");
    options.Add("C. bababa");
    options.Add("D. dddddd");
    answers.Add(2);
    givens.Add(-1);

    questions.Add("What's the next q?");
    options.Add("A. WHAAA");
    options.Add("B. WHOO");
    options.Add("C. bababa");
    options.Add("D. dddddd");
    answers.Add(2);
    givens.Add(-1);

    questions.Add("What's the q?");
    options.Add("A. WHAAA");
    options.Add("B. WHOO");
    options.Add("C. bababa");
    options.Add("D. dddddd");
    answers.Add(2);
    givens.Add(-1);

    questions.Add("What's the q?");
    options.Add("A. WHAAA");
    options.Add("B. WHOO");
    options.Add("C. bababa");
    options.Add("D. dddddd");
    answers.Add(2);
    givens.Add(-1);

    questions.Add("What's the q?");
    options.Add("A. WHAAA");
    options.Add("B. WHOO");
    options.Add("C. bababa");
    options.Add("D. dddddd");
    answers.Add(2);
    givens.Add(-1);

    option_quizos = new List<Quizo>();
    option_quizos.Add(GameObject.Find("QA").GetComponent<Quizo>());
    option_quizos.Add(GameObject.Find("QB").GetComponent<Quizo>());
    option_quizos.Add(GameObject.Find("QC").GetComponent<Quizo>());
    option_quizos.Add(GameObject.Find("QD").GetComponent<Quizo>());
    qconfirm_quizo = GameObject.Find("QConfirm").GetComponent<Quizo>();
    qboard = GameObject.Find("Qboard");
    qboard_lazerable = qboard.GetComponent<Lazerable>();
    qtext_tmp = GameObject.Find("Qtext").GetComponent<TextMeshPro>();
    SetQuizText();
  }

  void SetQuizText()
  {
    qtext_tmp.SetText(questions[question]);
    for(int i = 0; i < 4; i++) option_quizos[i].tmp.SetText(options[question*4+i]);
  }

  /*
  tried during:
  - newly snapped on
  - dial altered
  */
  void TryApplyTool(Tool t)
  {
    if(t == tool_insulator)
    {
      applied_heat = 0;
      if(t.engaged) return;
      if(tool_burner.engaged) applied_heat += tool_burner.dial_dial.map;
      if(tool_coil.engaged)   applied_heat -= tool_coil.dial_dial.map;
    }
    else if(t == tool_burner || t == tool_coil) //NO IMMEDIATE EFFECT
    {
      //visual
      if(t == tool_burner)
      {
        var vel = flame.velocityOverLifetime;
        vel.speedModifierMultiplier = Mathf.Lerp(0.1f,0.5f,t.dial_dial.val);
      }
      else if(t == tool_coil)
      {
        //visually update coil
      }

      applied_heat = 0;
      if(tool_insulator.engaged) return; //heat does nothing!
      if(tool_burner.engaged) applied_heat += tool_burner.dial_dial.map;
      if(tool_coil.engaged)   applied_heat -= tool_coil.dial_dial.map;
      //math change happens passively- only need to update visuals
    }
    else if(t == tool_clamp)
    {
      applied_weight = 0;
      if(t.engaged) return;
      if(tool_weight.engaged)  applied_weight += tool_weight.dial_dial.map;
      if(tool_balloon.engaged) applied_weight -= tool_balloon.dial_dial.map;
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
      t.active_available.transform.localScale  = scale;
      t.active_snap.transform.localScale       = scale;
      t.storage_available.transform.localScale = scale;
      t.storage_snap.transform.localScale      = scale;

      //math
      applied_weight = 0;
      if(tool_clamp.engaged) return; //weight does nothing!
      if(tool_weight.engaged)  applied_weight += tool_weight.dial_dial.map;
      if(tool_balloon.engaged) applied_weight -= tool_balloon.dial_dial.map;
      //math change happens passively- only need to update visuals
    }
    else if(t == tool_clipboard)
    {
      ; //do nothing!
    }

  }

  void TryAct(GameObject actable, float z_val, float y_val, ref float r_z, ref float r_y)
  {
    //grabbing handle
    if(actable == handle_workspace)
    {
      float dy = (r_y-y_val);
      workspace.transform.position = new Vector3(workspace.transform.position.x,workspace.transform.position.y-dy,workspace.transform.position.z);
    }
    else
    {
      Dial d = actable.GetComponent<Dial>();
      //grabbing dial
      if(d != null)
      {
        Tool t = d.tool.GetComponent<Tool>();
        float dx = (r_z-z_val)*2.0f;
        d.val = Mathf.Clamp(d.val-dx,0.0f,1.0f);

        TryApplyTool(t);
      }
    }
  }

  void TryGrab(bool which, float htrigger_val, float itrigger_val, float z_val, float y_val, Vector3 hand_vel, ref bool r_htrigger, ref bool r_itrigger, ref int r_htrigger_delta, ref int r_itrigger_delta, ref float r_z, ref float r_y, ref GameObject r_hand, ref GameObject r_grabbed, ref GameObject r_ohand, ref GameObject r_ograbbed)
  {
    float htrigger_threshhold = 0.1f;
    float itrigger_threshhold = 0.1f;

    //find deltas
    r_htrigger_delta = 0;
    if(!r_htrigger && htrigger_val > htrigger_threshhold)
    {
      r_htrigger_delta = 1;
      r_htrigger = true;
    }
    else if(r_htrigger && htrigger_val <= htrigger_threshhold)
    {
      r_htrigger_delta = -1;
      r_htrigger = false;
    }

    r_itrigger_delta = 0;
    if(!r_itrigger && itrigger_val > itrigger_threshhold)
    {
      r_itrigger_delta = 1;
      r_itrigger = true;
      r_z = z_val;
      r_y = y_val;
    }
    else if(r_itrigger && itrigger_val <= itrigger_threshhold)
    {
      r_itrigger_delta = -1;
      r_itrigger = false;
    }

    //find new grabs
    if(r_grabbed == null && r_htrigger_delta == 1)
    {
      //first try movables
      for(int i = 0; r_grabbed == null && i < movables.Count; i++)
      {
        if(
           ( which && movables[i].ltouch) ||
           (!which && movables[i].rtouch)
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
            TryApplyTool(t);
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
             ( which && tools[i].dial_touchable.ltouch) ||
             (!which && tools[i].dial_touchable.rtouch)
            )
          {
            r_grabbed = tools[i].dial;
            tools[i].dial_touchable.grabbed = true;
            if(r_grabbed == r_ograbbed) r_ograbbed = null;
          }
        }
      }
      //then extraaneous
      if(r_grabbed == null)
      {
        Touchable g = handle_workspace_touchable;
        if(
          ( which && g.ltouch) ||
          (!which && g.rtouch)
        )
        {
          r_grabbed = handle_workspace;
          g.grabbed = true;
          if(r_grabbed == r_ograbbed) r_ograbbed = null;
        }
      }
      if(r_grabbed == null)
      {
        Touchable g = halfer_touchable;
        if(
          ( which && g.ltouch) ||
          (!which && g.rtouch)
        )
        {
          r_grabbed = halfer;
          g.grabbed = true;
          if(r_grabbed == r_ograbbed) r_ograbbed = null;
        }
      }
    }
    //find new releases
    else if(r_grabbed && r_htrigger_delta == -1)
    {
      Tool t = r_grabbed.GetComponent<Tool>();
      if(t)
      {
        //tool newly engaged
        if(t.active_ghost.tintersect)
        {
          r_grabbed.transform.SetParent(t.active.transform);
          t.engaged = true;
          t.boxcollider.isTrigger = true;
               if(t == tool_insulator) t.dial_dial.val = (float)ThermoMath.percent_given_t(thermo.temperature);
          else if(t == tool_clamp)     t.dial_dial.val = (float)ThermoMath.percent_given_v(thermo.volume);
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
          r_grabbed.transform.SetParent(r_grabbed.GetComponent<Touchable>().og_parent);
          t.rigidbody.isKinematic = false;
          t.rigidbody.velocity = hand_vel;
        }
      }
      else
      {
        r_grabbed.transform.SetParent(r_grabbed.GetComponent<Touchable>().og_parent); //ok to do, even with a dial
        VisAid v = r_grabbed.GetComponent<VisAid>();
        if(v)
        {
          v.rigidbody.isKinematic = false;
          v.rigidbody.velocity = hand_vel;
        }
      }

      r_grabbed.GetComponent<Touchable>().grabbed = false;
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

    Touchable gr;
    bool ltouch = false;
    bool rtouch = false;
    for(int i = 0; i < movables.Count; i++)
    {
      gr = movables[i];
      if(gr.ltouch) ltouch = true;
      if(gr.rtouch) rtouch = true;
    }
    for(int i = 0; i < tools.Count; i++)
    {
      gr = tools[i].dial_touchable;
      if(gr.ltouch) ltouch = true;
      if(gr.rtouch) rtouch = true;
    }
    gr = handle_workspace_touchable;
    if(gr.ltouch) ltouch = true;
    if(gr.rtouch) rtouch = true;
    gr = halfer_touchable;
    if(gr.ltouch) ltouch = true;
    if(gr.rtouch) rtouch = true;

         if(lgrabbed) lhand_meshrenderer.material = hand_grabbing;
    else if(ltouch)   lhand_meshrenderer.material = hand_touching;
    else              lhand_meshrenderer.material = hand_empty;

         if(rgrabbed) rhand_meshrenderer.material = hand_grabbing;
    else if(rtouch)   rhand_meshrenderer.material = hand_touching;
    else              rhand_meshrenderer.material = hand_empty;
  }

  // Update is called once per frame
  void Update()
  {
    //passive effects
    if(applied_weight != 0) thermo.add_weight(applied_weight);
    if(applied_heat   != 0)
    {
      if(tool_clamp.engaged) thermo.add_heat_constant_v(applied_heat);
      else                   thermo.add_heat_constant_p(applied_heat);
    }

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

    float lhandt  = OVRInput.Get(OVRInput.Axis1D.PrimaryHandTrigger);
    float lindext = OVRInput.Get(OVRInput.Axis1D.PrimaryIndexTrigger);
    float rhandt  = OVRInput.Get(OVRInput.Axis1D.SecondaryHandTrigger);
    float rindext = OVRInput.Get(OVRInput.Axis1D.SecondaryIndexTrigger);
    //index compatibility
    lhandt += lindext;
    rhandt += rindext;
    TryGrab(true,  lhandt, lindext, lhand.transform.position.z, lhand.transform.position.y, lhand_vel, ref lhtrigger, ref litrigger, ref lhtrigger_delta, ref litrigger_delta, ref lz, ref ly, ref lhand, ref lgrabbed, ref rhand, ref rgrabbed); //left hand
    TryGrab(false, rhandt, rindext, rhand.transform.position.z, rhand.transform.position.y, rhand_vel, ref rhtrigger, ref ritrigger, ref rhtrigger_delta, ref ritrigger_delta, ref rz, ref ry, ref rhand, ref rgrabbed, ref lhand, ref lgrabbed); //right hand

    if(
      (litrigger_delta > 0 && halfer_touchable.ltouch) ||
      (ritrigger_delta > 0 && halfer_touchable.rtouch)
      )
    {
      halfed = !halfed;
      for(int i = 0; i < halfables.Count; i++)
        halfables[i].setHalf(halfed);
    }

    UpdateGrabVis();

    //quiz
    if(qboard_lazerable.lhit || vrcenter_lazerable.lhit) llazer_fadable.set_factive(true);
    else                                                 llazer_fadable.set_factive(false);
    if(qboard_lazerable.rhit || vrcenter_lazerable.rhit) rlazer_fadable.set_factive(true);
    else                                                 rlazer_fadable.set_factive(false);

    if(vrcenter_lazerable.hit)
    {
      if(
        (litrigger_delta == 1 && vrcenter_lazerable.lhit) ||
        (ritrigger_delta == 1 && vrcenter_lazerable.rhit)
      )
      {
        vrcenter_backing_meshrenderer.material = quiz_hisel;
        UnityEngine.XR.InputTracking.Recenter();
        OVRManager.display.RecenterPose();
        Vector3 pos = ceye.transform.localPosition*-1.0f;
        pos.y = 0.0f;
        cam_offset.transform.position = pos;
      }
      else vrcenter_backing_meshrenderer.material = quiz_hi;
    }
    else vrcenter_backing_meshrenderer.material = quiz_default;

    if(qboard_lazerable.lhit || vrcenter_lazerable.lhit)
    {
      float d = 0f;
           if(qboard_lazerable.lhit)   d = Vector3.Distance(lhand.transform.position,qboard.transform.position);
      else if(vrcenter_lazerable.lhit) d = Vector3.Distance(lhand.transform.position,vrcenter.transform.position);
      llazer_viz.transform.localPosition = new Vector3(llazer_viz.transform.localPosition.x,llazer_viz.transform.localPosition.y,d/2);
      llazer_viz.transform.localScale    = new Vector3(llazer_viz.transform.localScale.x,d/2,llazer_viz.transform.localScale.z);
    }
    if(qboard_lazerable.rhit || vrcenter_lazerable.rhit)
    {
      float d = 0f;
           if(qboard_lazerable.rhit)   d = Vector3.Distance(lhand.transform.position,qboard.transform.position);
      else if(vrcenter_lazerable.rhit) d = Vector3.Distance(lhand.transform.position,vrcenter.transform.position);
      rlazer_viz.transform.localPosition = new Vector3(rlazer_viz.transform.localPosition.x,rlazer_viz.transform.localPosition.y,d/2);
      rlazer_viz.transform.localScale    = new Vector3(rlazer_viz.transform.localScale.x,d/2,rlazer_viz.transform.localScale.z);
    }

    if(!llazer_fadable.stale)
    {
      if(llazer_fadable.alpha == 0.0f) llazer_meshrenderer.enabled = false;
      else
      {
        llazer_meshrenderer.enabled = true;
        Color c = llazer_meshrenderer.material.color;
        llazer_meshrenderer.material.color = new Color(c.r,c.g,c.b,llazer_fadable.alpha);
      }
    }
    if(!rlazer_fadable.stale)
    {
      if(rlazer_fadable.alpha == 0.0f) rlazer_meshrenderer.enabled = false;
      else
      {
        rlazer_meshrenderer.enabled = true;
        Color c = rlazer_meshrenderer.material.color;
        rlazer_meshrenderer.material.color = new Color(c.r,c.g,c.b,rlazer_fadable.alpha);
      }
    }

    Quizo q;
    for(int i = 0; i < option_quizos.Count; i++)
    {
      q = option_quizos[i];
      if(
        (litrigger_delta == 1 && q.lazerable.lhit) ||
        (ritrigger_delta == 1 && q.lazerable.rhit)
      )
      {
        if(i == qselected) qselected = -1;
        else               qselected = i;
      }
      if(i == qselected)
      {
        if(q.lazerable.hit) q.backing_meshrenderer.material = quiz_hisel;
        else                q.backing_meshrenderer.material = quiz_sel;
      }
      else
      {
        if(q.lazerable.hit) q.backing_meshrenderer.material = quiz_hi;
        else                q.backing_meshrenderer.material = quiz_default;
      }
    }
    if(qselected != -1)
    {
      if(qconfirm_quizo.lazerable.hit) qconfirm_quizo.backing_meshrenderer.material = quiz_hisel;
      else                             qconfirm_quizo.backing_meshrenderer.material = quiz_sel;
    }
    else
    {
      if(qconfirm_quizo.lazerable.hit) qconfirm_quizo.backing_meshrenderer.material = quiz_hi;
      else                             qconfirm_quizo.backing_meshrenderer.material = quiz_default;
    }

    //tooltext
    Tool t;
    for(int i = 0; i < tools.Count; i++)
    {
      t = tools[i];
      t.dial_dial.examined = false;
      if(t.dial == lgrabbed || t.dial == rgrabbed) t.dial_dial.examined = true;
      if(t.dial_dial.val != t.dial_dial.prev_val)
      {
        t.textv_tmp.SetText("{0:3}"+t.dial_dial.unit,(float)t.dial_dial.map);
        t.dial_dial.examined = true;
      }
      t.dial_dial.prev_val = t.dial_dial.val;
    }

    for(int i = 0; i < tools.Count; i++)
    {
      t = tools[i];
      if(!t.text_fadable.stale)
      {
        if(t.text_fadable.alpha == 0.0f)
        {
          t.textv_meshrenderer.enabled = false;
          t.textl_meshrenderer.enabled = false;
        }
        else
        {
          t.textv_meshrenderer.enabled = true;
          t.textl_meshrenderer.enabled = true;
          Color32 c = new Color32(0,0,0,(byte)(t.text_fadable.alpha*255));
          t.textv_tmp.faceColor = c;
          t.textl_tmp.faceColor = c;
        }
      }
    }

/*
//overwrite insulator's text for debugging purposes
{
    Tool to = tools[0];
    to.textv_tmp.SetText("{0}",litrigger_delta);
    to.dial_dial.examined = true;
    to.dial_dial.prev_val = to.dial_dial.val*-1;
    to.textv_meshrenderer.enabled = true;
    Color32 c = new Color32(0,0,0,255);
    to.textv_tmp.faceColor = c;
    to.textl_tmp.faceColor = c;
}
*/

  }

}

