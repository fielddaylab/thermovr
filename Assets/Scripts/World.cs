/*
DOCUMENTATION- phil, 12/16/19

This class manages all the interaction in the scene.
It relies on ThermoState to keep track of any thermodynamic-centric state, but other than that, this is responsible for everything moving about the scene.
It should be instantiated as a game object "Oracle" at the root of the scene heirarchy.
There are unfortunately somewhat inconsistent patterns of what variables are defined publicly via the editor inspector, and which are set in code, though I tried to err toward the latter where possible.
*/

using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using TMPro;

public class World : MonoBehaviour
{
  public Material hand_empty;
  Material[] hand_emptys;
  public Material hand_touching;
  Material[] hand_touchings;
  public Material hand_grabbing;
  Material[] hand_grabbings;
  public Material tab_default;
  public Material tab_hi;
  public Material tab_sel;
  public Material tab_hisel;

  ThermoState thermo;
  GameObject cam_offset;
  GameObject ceye;
  GameObject lhand;
  GameObject lactualhand;
  Vector3 lhand_pos;
  Vector3 lhand_vel;
  SkinnedMeshRenderer lhand_meshrenderer;
  GameObject rhand;
  GameObject ractualhand;
  Vector3 rhand_pos;
  Vector3 rhand_vel;
  SkinnedMeshRenderer rhand_meshrenderer;

  GameObject vrcenter;
  FingerToggleable vrcenter_fingertoggleable;
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
  float lz = 0f;
  float rz = 0f;
  float ly = 0f;
  float ry = 0f;

  List<Tool> tools;
  Tool tool_insulator;
  Tool tool_clamp;
  Tool tool_burner;
  Tool tool_coil;
  Tool tool_weight;
  Tool tool_balloon;

  ParticleSystem flame; //special case

  bool halfed = false;
  List<Halfable> halfables;
  GameObject halfer;
  Touchable halfer_touchable;

  GameObject vessel;
  GameObject graph;
  GameObject state_dot;
  GameObject challenge_dot;
  GameObject clipboard;

  ChallengeBall challenge_ball_collide;

  bool lhtrigger = false;
  bool rhtrigger = false;
  bool litrigger = false;
  bool ritrigger = false;

  GameObject instructions_parent;
  GameObject challenge_parent;
  GameObject quiz_parent;
  List<Tab> mode_tabs;
  int board_mode = 0; //0- instructions, 1- challenge, 2- quiz, 3- congrats
  int question = 0;
  List<string> questions;
  List<string> options;
  List<int> answers;
  List<int> givens;
  List<Tab> option_tabs;
  Tab qconfirm_tab;
  GameObject board;
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

    hand_emptys    = new Material[]{hand_empty};
    hand_touchings = new Material[]{hand_touching};
    hand_grabbings = new Material[]{hand_grabbing};

    lhand  = GameObject.Find("LeftControllerAnchor");
    lactualhand = lhand.transform.GetChild(0).gameObject;
    lhand_pos = lhand.transform.position;
    lhand_vel = new Vector3(0f,0f,0f);
    lhand_meshrenderer = GameObject.Find("hands:Lhand").GetComponent<SkinnedMeshRenderer>();

    rhand  = GameObject.Find("RightControllerAnchor");
    ractualhand = rhand.transform.GetChild(0).gameObject;
    rhand_pos = rhand.transform.position;
    rhand_vel = new Vector3(0f,0f,0f);
    rhand_meshrenderer = GameObject.Find("hands:Rhand").GetComponent<SkinnedMeshRenderer>();

    lhand_meshrenderer.materials = hand_emptys;
    rhand_meshrenderer.materials = hand_emptys;

    Tool t;
    tools = new List<Tool>();
    t = GameObject.Find("Tool_Insulator").GetComponent<Tool>(); tool_insulator = t; tools.Add(t); t.dial_dial.min_map =  0f; t.dial_dial.max_map = 1f; t.dial_dial.unit = "n";
    t = GameObject.Find("Tool_Clamp"    ).GetComponent<Tool>(); tool_clamp     = t; tools.Add(t); t.dial_dial.min_map =  0f; t.dial_dial.max_map = 1f; t.dial_dial.unit = "h";
    t = GameObject.Find("Tool_Burner"   ).GetComponent<Tool>(); tool_burner    = t; tools.Add(t); t.dial_dial.min_map =  1f; t.dial_dial.max_map =  1000f*100f; t.dial_dial.unit = "J/s";
    t = GameObject.Find("Tool_Coil"     ).GetComponent<Tool>(); tool_coil      = t; tools.Add(t); t.dial_dial.min_map = -1f; t.dial_dial.max_map = -1000f*100f; t.dial_dial.unit = "J/s";
    double kg_corresponding_to_10mpa = thermo.surfacearea*1550/*M^2->in^2*/*(10*1453.8/*MPa->psi*/)*0.453592/*lb->kg*/;
    t = GameObject.Find("Tool_Weight"   ).GetComponent<Tool>(); tool_weight    = t; tools.Add(t); t.dial_dial.min_map =  0f; t.dial_dial.max_map =  (float)kg_corresponding_to_10mpa; t.dial_dial.unit = "kg";
    t = GameObject.Find("Tool_Balloon"  ).GetComponent<Tool>(); tool_balloon   = t; tools.Add(t); t.dial_dial.min_map =  0f; t.dial_dial.max_map = -(float)kg_corresponding_to_10mpa; t.dial_dial.unit = "kg";

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
      g.transform.localPosition = new Vector3(0f,0f,0f);
      g.transform.localScale = new Vector3(1f,1f,1f);
      g.transform.localRotation = Quaternion.identity;
      float v = t.storage.transform.localScale.x; //can grab any dimension
      Vector3 invscale = new Vector3(1f/v,1f/v,1f/v);
      t.text.transform.localScale = invscale;
      t.textv_tmp.SetText("{0:3}"+t.dial_dial.unit,(float)t.dial_dial.map);
    }

    vessel = GameObject.Find("Vessel");
    graph = GameObject.Find("Graph");
    state_dot = GameObject.Find("gstate");
    challenge_dot = GameObject.Find("cstate");
    clipboard = GameObject.Find("Clipboard");

    challenge_ball_collide = challenge_dot.GetComponent<ChallengeBall>();

    vrcenter = GameObject.Find("VRCenter");
    vrcenter_fingertoggleable = vrcenter.GetComponent<FingerToggleable>();
    vrcenter_backing_meshrenderer = vrcenter.transform.GetChild(1).GetComponent<MeshRenderer>();

    movables = new List<Touchable>();
    for(int i = 0; i < tools.Count; i++) movables.Add(tools[i].touchable); //important that tools take priority, so they can be grabbed and removed
    movables.Add(graph.GetComponent<Touchable>());
    movables.Add(clipboard.GetComponent<Touchable>());

    halfer = GameObject.Find("Halfer");
    halfer_touchable = halfer.GetComponent<Touchable>();
    halfables = new List<Halfable>();
    halfables.Add(GameObject.Find("Container"     ).GetComponent<Halfable>());
    halfables.Add(GameObject.Find("Tool_Insulator").GetComponent<Halfable>());
    halfables.Add(GameObject.Find("Tool_Coil"     ).GetComponent<Halfable>());

    instructions_parent = GameObject.Find("Instructions");
    challenge_parent = GameObject.Find("Challenge");
    quiz_parent = GameObject.Find("Quiz");

    mode_tabs = new List<Tab>();
    mode_tabs.Add(GameObject.Find("ModeInstructions").GetComponent<Tab>());
    mode_tabs.Add(GameObject.Find("ModeChallenge").GetComponent<Tab>());
    mode_tabs.Add(GameObject.Find("ModeQuiz").GetComponent<Tab>());

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

    option_tabs = new List<Tab>();
    option_tabs.Add(GameObject.Find("QA").GetComponent<Tab>());
    option_tabs.Add(GameObject.Find("QB").GetComponent<Tab>());
    option_tabs.Add(GameObject.Find("QC").GetComponent<Tab>());
    option_tabs.Add(GameObject.Find("QD").GetComponent<Tab>());
    qconfirm_tab = GameObject.Find("QConfirm").GetComponent<Tab>();
    board = GameObject.Find("Board");
    qtext_tmp = GameObject.Find("Qtext").GetComponent<TextMeshPro>();
    SetQuizText();
    SetChallengeBall();
    SetAllHalfed(true);
  }

  void SetQuizText()
  {
    qtext_tmp.SetText(questions[question]);
    for(int i = 0; i < 4; i++) option_tabs[i].tmp.SetText(options[question*4+i]);
  }

  void SetChallengeBall()
  {
    double volume      = ThermoMath.v_given_percent(Random.Range(0.1f,0.9f));
    double temperature = ThermoMath.t_given_percent(Random.Range(0.1f,0.9f));
    double pressure    = ThermoMath.p_given_vt(volume, temperature);
    challenge_dot.transform.localPosition = thermo.plot(pressure, volume, temperature);
  }

  void SetAllHalfed(bool h)
  {
    halfed = h;
    for(int i = 0; i < halfables.Count; i++)
      halfables[i].setHalf(halfed);
    //special case, only halfed when engaged
    if(!tool_coil.engaged) tool_coil.gameObject.GetComponent<Halfable>().setHalf(false);
    if(!tool_insulator.engaged) tool_insulator.gameObject.GetComponent<Halfable>().setHalf(false);
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
      float v = 1f;
           if(t == tool_weight)  v += tool_weight.dial_dial.val;
      else if(t == tool_balloon) v += tool_balloon.dial_dial.val;
      Vector3 scale;
      Vector3 invscale;

      scale = new Vector3(v,v,v);
      invscale = new Vector3(1f/v,1f/v,1f/v);
      t.active.transform.localScale = scale;
      if(t.engaged) t.text.transform.localScale = invscale;

      v *= t.default_storage_scale;
      scale = new Vector3(v,v,v);
      invscale = new Vector3(1f/v,1f/v,1f/v);
      t.storage.transform.localScale = scale;
      if(t.stored) t.text.transform.localScale = invscale;

      //math
      applied_weight = 0;
      if(tool_clamp.engaged) return; //weight does nothing!
      if(tool_weight.engaged)  applied_weight += tool_weight.dial_dial.map;
      if(tool_balloon.engaged) applied_weight -= tool_balloon.dial_dial.map;
      //math change happens passively- only need to update visuals
    }

  }

  //safe to call if not interactable, as it will just do nothing
  void TryInteractable(GameObject actable, float x_val, float y_val, ref float r_x, ref float r_y)
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
        float dx = (r_x-x_val)*-10f;
        d.val = Mathf.Clamp(d.val-dx,0f,1f);

        TryApplyTool(t);
      }
    }
  }

  //"which": true -> left, false -> right
  void TryHand(bool which, float htrigger_val, float itrigger_val, float x_val, float y_val, Vector3 hand_vel, ref bool r_htrigger, ref bool r_itrigger, ref int r_htrigger_delta, ref int r_itrigger_delta, ref float r_x, ref float r_y, ref GameObject r_hand, ref GameObject r_grabbed, ref GameObject r_ohand, ref GameObject r_ograbbed)
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
      r_x = x_val;
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
          r_grabbed.transform.SetParent(r_hand.transform);
          if(r_grabbed == r_ograbbed) r_ograbbed = null;
          movables[i].grabbed = true;
          Tool t = r_grabbed.GetComponent<Tool>();
          if(t)
          {
            t.engaged = false;
            t.stored = false;
            r_grabbed.transform.localScale = new Vector3(1f,1f,1f);
            t.text.transform.localScale = new Vector3(1f,1f,1f);
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
      if(r_grabbed != null) //newly grabbed
      {
        Halfable h = r_grabbed.GetComponent<Halfable>();
        if(h != null) h.setHalf(false); //nothing should be halfed while being grabbed
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
          r_grabbed.transform.localScale = new Vector3(1f,1f,1f);
          float v = t.active.transform.localScale.x; //can grab any dimension
          Vector3 invscale = new Vector3(1f/v,1f/v,1f/v);
          t.text.transform.localScale = invscale;
          Halfable h = r_grabbed.GetComponent<Halfable>();
          if(h != null) h.setHalf(halfed); //conform to half-ness while engaged
        }
        //tool newly stored
        else if(t.storage_ghost.tintersect)
        {
          r_grabbed.transform.SetParent(t.storage.transform);
          t.stored = true;
          r_grabbed.transform.localPosition = new Vector3(0f,0f,0f);
          r_grabbed.transform.localRotation = Quaternion.identity;
          r_grabbed.transform.localScale = new Vector3(1f,1f,1f);
          float v = t.storage.transform.localScale.x; //can grab any dimension
          Vector3 invscale = new Vector3(1f/v,1f/v,1f/v);
          t.text.transform.localScale = invscale;
        }
        //tool released
        else
        {
          r_grabbed.transform.SetParent(r_grabbed.GetComponent<Touchable>().og_parent);
          r_grabbed.transform.localScale = new Vector3(1f,1f,1f);
          t.text.transform.localScale = new Vector3(1f,1f,1f);
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

    if(r_grabbed) TryInteractable(r_grabbed, x_val, y_val, ref r_x, ref r_y);

    r_x = x_val;
    r_y = y_val;


    //halfer
    if(
      r_itrigger_delta > 0 &&
      (
        ( which && halfer_touchable.ltouch) ||
        (!which && halfer_touchable.rtouch)
      )
    )
      SetAllHalfed(!halfed);

    //centerer
    if(vrcenter_fingertoggleable.finger)
    {
      if(
        ( which && vrcenter_fingertoggleable.lfinger) ||
        (!which && vrcenter_fingertoggleable.rfinger)
      )
      {
        vrcenter_backing_meshrenderer.material = tab_hisel;
        UnityEngine.XR.InputTracking.Recenter();
        OVRManager.display.RecenterPose();
        Vector3 pos = cam_offset.transform.localPosition-(cam_offset.transform.localPosition+ceye.transform.localPosition);
        pos.y = 0f;
        cam_offset.transform.localPosition = pos;
      }
      else vrcenter_backing_meshrenderer.material = tab_hi;
    }
    else vrcenter_backing_meshrenderer.material = tab_default;

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
          t.active_snap_meshrenderer.enabled      = true;
        }
        else
        {
          t.active_available_meshrenderer.enabled = true;
          t.active_snap_meshrenderer.enabled      = false;
        }
        //storage
        if(t.storage_ghost.tintersect)
        {
          t.storage_available_meshrenderer.enabled = false;
          t.storage_snap_meshrenderer.enabled      = true;
        }
        else
        {
          t.storage_available_meshrenderer.enabled = true;
          t.storage_snap_meshrenderer.enabled      = false;
        }
      }
      else
      {
        t.active_snap_meshrenderer.enabled      = false;
        t.active_available_meshrenderer.enabled = false;
        t.storage_snap_meshrenderer.enabled      = false;
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

         if(lgrabbed) lhand_meshrenderer.materials = hand_grabbings;
    else if(ltouch)   lhand_meshrenderer.materials = hand_touchings;
    else              lhand_meshrenderer.materials = hand_emptys;

         if(rgrabbed) rhand_meshrenderer.materials = hand_grabbings;
    else if(rtouch)   rhand_meshrenderer.materials = hand_touchings;
    else              rhand_meshrenderer.materials = hand_emptys;
  }

  //give it a list of fingertoggleables, and it manipulates them to act as a singularly-selectable list
  int reconcileDependentSelectables(int known, List<Tab> list)
  {
    int n_toggled = 0;
    Tab t;
    for(int i = 0; i < list.Count; i++)
    {
      t = list[i];
      if(t.fingertoggleable.on) n_toggled++;
    }

    if(n_toggled <= 1)
    {
      known = -1;
      for(int i = 0; i < list.Count; i++)
      {
        t = list[i];
        if(t.fingertoggleable.on) known = i;
      }
    }
    else //need conflict resolution!
    {
      known = -1;
      for(int i = 0; i < list.Count; i++)
      {
        t = list[i];
        if(t.fingertoggleable.on)
        {
          if(known == -1) known = i;
          else
          {
            //if only t is intersecting, prefer t
            if(t.fingertoggleable.finger && !list[known].fingertoggleable.finger)
            {
              list[known].fingertoggleable.on = false;
              known = i;
            }
            else //prefer previous (ignore t)
              t.fingertoggleable.on = false;
          }
        }
      }
    }

    return known;
  }

  void updateSelectableVis(int known, List<Tab> list)
  {
    Tab t;
    for(int i = 0; i < list.Count; i++)
    {
      t = list[i];
      if(known == i)
      {
        if(t.fingertoggleable.finger) t.backing_meshrenderer.material = tab_hisel;
        else                          t.backing_meshrenderer.material = tab_sel;
      }
      else
      {
        if(t.fingertoggleable.finger) t.backing_meshrenderer.material = tab_hi;
        else                          t.backing_meshrenderer.material = tab_default;
      }
    }
  }

  float hack_timer = 0f;
  void Update()
  {
    hack_timer += Time.deltaTime;
    while(hack_timer > 100f) hack_timer -= 100f;

    //hands keep trying to run away- no idea why (this is a silly way to keep them still)
    lactualhand.transform.localPosition    = new Vector3(0f,0f,0f);
    lactualhand.transform.localEulerAngles = new Vector3(0f,0f,90f);//localRotation = Quaternion.identity;
    ractualhand.transform.localPosition    = new Vector3(0f,0f,0f);
    ractualhand.transform.localEulerAngles = new Vector3(0f,0f,-90f);//localRotation = Quaternion.identity;

    //passive effects
    if(applied_weight != 0) thermo.add_pressure(applied_weight); //TODO: must convert weight to pressure!
    if(applied_heat   != 0)
    {
      if(tool_clamp.engaged) thermo.add_heat_constant_v(applied_heat); //TODO: make sure "applied heat" is first converted to appropriate unit!
      else                   thermo.add_heat_constant_p(applied_heat); //TODO: make sure "applied heat" is first converted to appropriate unit!
    }

    //running blended average of hand velocity, for consistent "throwing"
    lhand_vel += (lhand.transform.position-lhand_pos)/Time.deltaTime;
    lhand_vel *= 0.5f;
    lhand_pos = lhand.transform.position;

    rhand_vel += (rhand.transform.position-rhand_pos)/Time.deltaTime;
    rhand_vel *= 0.5f;
    rhand_pos = rhand.transform.position;

    float lhandt  = OVRInput.Get(OVRInput.Axis1D.PrimaryHandTrigger);
    float lindext = OVRInput.Get(OVRInput.Axis1D.PrimaryIndexTrigger);
    float rhandt  = OVRInput.Get(OVRInput.Axis1D.SecondaryHandTrigger);
    float rindext = OVRInput.Get(OVRInput.Axis1D.SecondaryIndexTrigger);
    //index compatibility
    if(OVRInput.Get(OVRInput.Button.One,OVRInput.Controller.LTouch)) lhandt = 1.0f;
    if(OVRInput.Get(OVRInput.Button.One,OVRInput.Controller.RTouch)) rhandt = 1.0f;
    // update 12/19/19- ovr just doesn't recognize index input. so hacking a timed squeeze/release for testing
    //if((int)hack_timer%2 == 1) { lindext = 1f; rindext = 1f; }
    lhandt += lindext;
    rhandt += rindext;
    TryHand(true,  lhandt, lindext, lhand.transform.position.x, lhand.transform.position.y, lhand_vel, ref lhtrigger, ref litrigger, ref lhtrigger_delta, ref litrigger_delta, ref lz, ref ly, ref lhand, ref lgrabbed, ref rhand, ref rgrabbed); //left hand
    TryHand(false, rhandt, rindext, rhand.transform.position.x, rhand.transform.position.y, rhand_vel, ref rhtrigger, ref ritrigger, ref rhtrigger_delta, ref ritrigger_delta, ref rz, ref ry, ref rhand, ref rgrabbed, ref lhand, ref lgrabbed); //right hand

    UpdateGrabVis();

    //clipboard
    int old_board_mode = board_mode;
    board_mode = reconcileDependentSelectables(board_mode, mode_tabs);
    if(board_mode == -1) board_mode = old_board_mode;
    updateSelectableVis(board_mode, mode_tabs);

    switch(board_mode)
    {
      case 0: //instructions
        if(!instructions_parent.activeSelf) instructions_parent.SetActive(true);
        if( challenge_parent.activeSelf)    challenge_parent.SetActive(   false);
        if( quiz_parent.activeSelf)         quiz_parent.SetActive(        false);
        break;
      case 1: //challenge
        if( instructions_parent.activeSelf) instructions_parent.SetActive(false);
        if(!challenge_parent.activeSelf)    challenge_parent.SetActive(   true);
        if( quiz_parent.activeSelf)         quiz_parent.SetActive(        false);
        break;
      case 2: //quiz
        if( instructions_parent.activeSelf) instructions_parent.SetActive(false);
        if( challenge_parent.activeSelf)    challenge_parent.SetActive(   false);
        if(!quiz_parent.activeSelf)         quiz_parent.SetActive(        true);
        qselected = reconcileDependentSelectables(qselected, option_tabs);
        updateSelectableVis(qselected, option_tabs);
        if(qselected != -1)
        {
          if(qconfirm_tab.fingertoggleable.finger) qconfirm_tab.backing_meshrenderer.material = tab_hisel;
          else                                     qconfirm_tab.backing_meshrenderer.material = tab_sel;
        }
        else
        {
          if(qconfirm_tab.fingertoggleable.finger) qconfirm_tab.backing_meshrenderer.material = tab_hi;
          else                                     qconfirm_tab.backing_meshrenderer.material = tab_default;
        }
        break;
      case 3: //congratulations
        board_mode = 2; //immediately return back to quiz until this section is actually implemented
        break;
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
        if(t.text_fadable.alpha == 0f)
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

  }

}

