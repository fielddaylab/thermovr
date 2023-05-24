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
    const float CONTAINER_INSULATION_COEFFICIENT = 0.1f; // Not really based on a physical material, just a way to roughly simulate imperfect insulation.
    const bool QUIZ_ENABLED = false;
    const bool CHALLENGE_ENABLED = false;
    public const double DELTA_PRESSURE_CUTOFF = 100.0;

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
    public DirectionalIndicator arrows;

    // ThermoState thermo;
    ThermoPresent thermo_present;
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
    Touchable graph_touchable;
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
    Tool tool_ambientPressure;
    Tool tool_roomTemp;
    Tool tool_percentInsulation;
    Tool tool_clampRange;

    ParticleSystem flame; //special case

    bool halfed = false;
    List<Halfable> halfables;
    GameObject halfer;
    Touchable halfer_touchable;
    FingerToggleable halfer_fingertoggleable;
    GameObject reset;
    Touchable reset_touchable;
    FingerToggleable reset_fingertoggleable;

    GameObject vessel;
    GameObject graph;
    GameObject state_dot;
    GameObject placement_dot;
    Vector3 placement_thermo;
    bool placement_thermo_reasonable;
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
    bool qconfirmed = false;
    GameObject board;
    TextMeshPro qtext_tmp;
    int qselected = -1;

    double room_temp = 72; // in F
    double applied_heat = 0;
    double applied_weight = 0;
    double ambient_pressure = 0;

    // Start is called before the first frame update
    void Start() {
        // All this code does, at end of day, is find all the objects to manage,
        // and set initial values and such as needed.
        thermo_present = GameObject.Find("Oracle").GetComponent<ThermoPresent>();

        cam_offset = GameObject.Find("CamOffset");
        ceye = GameObject.Find("CenterEyeAnchor");

        hand_emptys = new Material[] { hand_empty };
        hand_touchings = new Material[] { hand_touching };
        hand_grabbings = new Material[] { hand_grabbing };

        lhand = GameObject.Find("LeftControllerAnchor");
        lactualhand = lhand.transform.GetChild(0).gameObject;
        lhand_pos = lhand.transform.position;
        lhand_vel = new Vector3(0f, 0f, 0f);
        lhand_meshrenderer = GameObject.Find("hands:Lhand").GetComponent<SkinnedMeshRenderer>();

        rhand = GameObject.Find("RightControllerAnchor");
        ractualhand = rhand.transform.GetChild(0).gameObject;
        rhand_pos = rhand.transform.position;
        rhand_vel = new Vector3(0f, 0f, 0f);
        rhand_meshrenderer = GameObject.Find("hands:Rhand").GetComponent<SkinnedMeshRenderer>();

        lhand_meshrenderer.materials = hand_emptys;
        rhand_meshrenderer.materials = hand_emptys;

        // As we grab them, set ranges on tool dials (sliders).
        Tool t;
        tools = new List<Tool>();
        t = GameObject.Find("Tool_Insulator").GetComponent<Tool>(); tool_insulator = t; tools.Add(t); t.dial_dial.min_map = 0f; t.dial_dial.max_map = 1f; t.dial_dial.unit = ""; t.dial_dial.display_unit = t.dial_dial.unit;
        t = GameObject.Find("Tool_Clamp").GetComponent<Tool>(); tool_clamp = t; tools.Add(t); t.dial_dial.min_map = 0f; t.dial_dial.max_map = 1f; t.dial_dial.unit = ""; t.dial_dial.display_unit = t.dial_dial.unit;
        t = GameObject.Find("Tool_Burner").GetComponent<Tool>(); tool_burner = t; tools.Add(t); t.dial_dial.min_map = 0f; t.dial_dial.max_map = 1000f * 100f; t.dial_dial.unit = "J/s"; t.dial_dial.display_unit = "kJ/s"; t.dial_dial.display_mul = 0.001f;
        t = GameObject.Find("Tool_Coil").GetComponent<Tool>(); tool_coil = t; tools.Add(t); t.dial_dial.min_map = 0f; t.dial_dial.max_map = -1000f * 100f; t.dial_dial.unit = "J/s"; t.dial_dial.display_unit = "kJ/s"; t.dial_dial.display_mul = 0.001f;
        double kg_corresponding_to_10mpa = thermo_present.get_surfacearea_insqr() * (10 * 1453.8/*MPa->psi*/) * 0.453592/*lb->kg*/;
        double kg_corresponding_to_2mpa = thermo_present.get_surfacearea_insqr() * (2 * 1453.8/*MPa->psi*/) * 0.453592/*lb->kg*/; // 10 MPa seems way too big, sooooo... we'll just do 2 MPa.
        t = GameObject.Find("Tool_Weight").GetComponent<Tool>(); tool_weight = t; tools.Add(t); t.dial_dial.min_map = 0f; t.dial_dial.max_map = (float)kg_corresponding_to_10mpa; t.dial_dial.unit = "kg"; t.dial_dial.display_unit = t.dial_dial.unit;
        t = GameObject.Find("Tool_Balloon").GetComponent<Tool>(); tool_balloon = t; tools.Add(t); t.dial_dial.min_map = 0f; t.dial_dial.max_map = -(float)kg_corresponding_to_10mpa / 500.0f; t.dial_dial.unit = "kg"; t.dial_dial.display_unit = t.dial_dial.unit;
        // TODO: establish logical bounds and units on the ambient pressure dial
        t = GameObject.Find("Tool_AmbientPressure").GetComponent<Tool>(); tool_ambientPressure = t; tools.Add(t); t.dial_dial.min_map = 0f; t.dial_dial.max_map = (float)kg_corresponding_to_2mpa / 100; t.dial_dial.unit = "kPa"; t.dial_dial.display_unit = t.dial_dial.unit;
        // TODO: allow for non-zero starting values / positions; apply relevant heat transfer between room and piston
        t = GameObject.Find("Tool_RoomTemp").GetComponent<Tool>(); tool_roomTemp = t; tools.Add(t); t.dial_dial.min_map = -100; t.dial_dial.max_map = 200; t.dial_dial.unit = "F"; t.dial_dial.display_unit = "F";
        t = GameObject.Find("Tool_PercentInsulation").GetComponent<Tool>(); tool_percentInsulation = t; tools.Add(t); t.dial_dial.min_map = 0f; t.dial_dial.max_map = 100; t.dial_dial.unit = "%"; t.dial_dial.display_unit = "%";
        // t = GameObject.Find("Tool_ClampRange").GetComponent<Tool>(); tool_burner = t; tools.Add(t); t.dial_dial.min_map = 0f; t.dial_dial.max_map = 1000f * 100f; t.dial_dial.unit = "J/s"; t.dial_dial.display_unit = "kJ/s"; t.dial_dial.display_mul = 0.001f;


        //Phil 06/23/2020: I can't figure out why tool_insulator and tool_clamp are both disabled. In the editor, they appear enabled until I hit play. I can't find anywhere in code that disables them. But they get to here, and are disabled, so I'm just re-enabling them
        tool_insulator.enabled = true;
        tool_clamp.enabled = true;

        flame = GameObject.Find("Flame").GetComponent<ParticleSystem>();

        workspace = GameObject.Find("Workspace");
        handle_workspace = GameObject.Find("Handle_Workspace");
        handle_workspace_touchable = handle_workspace.GetComponent<Touchable>();

        // set initial states of meshrenderers and transforms for our tools.
        for (int i = 0; i < tools.Count; i++) {
            t = tools[i];
            if (t.active_available_meshrenderer != null) {
                t.active_available_meshrenderer.enabled = false;
            }
            if (t.storage_meshrenderer != null) {
                t.storage_meshrenderer.enabled = false;
            }
            if (t.storage != null) {
                GameObject g = t.gameObject;
                g.transform.SetParent(t.storage.gameObject.transform);
                t.stored = true;
                g.transform.localPosition = new Vector3(0f, 0f, 0f);
                g.transform.localScale = new Vector3(1f, 1f, 1f);
                g.transform.localRotation = Quaternion.identity;
                float v = t.storage.transform.localScale.x; //can grab any dimension
                Vector3 invscale = new Vector3(1f / v, 1f / v, 1f / v);
                t.text.transform.localScale = invscale;
            }

            UpdateToolText(t);
        }

        vessel = GameObject.Find("Vessel");
        graph = GameObject.Find("Graph");
        graph_touchable = graph.GetComponent<Touchable>();
        state_dot = GameObject.Find("gstate");
        placement_dot = GameObject.Find("tstate");
        placement_dot.GetComponent<Renderer>().enabled = false;
        placement_thermo_reasonable = false;
        challenge_dot = GameObject.Find("cstate");
        clipboard = GameObject.Find("Clipboard");

        challenge_ball_collide = challenge_dot.GetComponent<ChallengeBall>();

        vrcenter = GameObject.Find("VRCenter");
        vrcenter_fingertoggleable = vrcenter.GetComponent<FingerToggleable>();
        vrcenter_backing_meshrenderer = vrcenter.transform.GetChild(1).GetComponent<MeshRenderer>();

        movables = new List<Touchable>();
        for (int i = 0; i < tools.Count; i++) movables.Add(tools[i].touchable); //important that tools take priority, so they can be grabbed and removed
        movables.Add(clipboard.GetComponent<Touchable>());

        halfer = GameObject.Find("Halfer");
        halfer_touchable = halfer.GetComponent<Touchable>();
        halfer_fingertoggleable = halfer.GetComponent<FingerToggleable>();
        reset = GameObject.Find("Reset");
        reset_touchable = reset.GetComponent<Touchable>();
        reset_fingertoggleable = reset.GetComponent<FingerToggleable>();
        halfables = new List<Halfable>();
        halfables.Add(GameObject.Find("Container").GetComponent<Halfable>());
        halfables.Add(GameObject.Find("Tool_Insulator").GetComponent<Halfable>());
        halfables.Add(GameObject.Find("Tool_Coil").GetComponent<Halfable>());

        // A bunch of stuff related to initializing clipboard.
        instructions_parent = GameObject.Find("Instructions");
        challenge_parent = GameObject.Find("Challenge");
        quiz_parent = GameObject.Find("Quiz");

        mode_tabs = new List<Tab>();
        mode_tabs.Add(GameObject.Find("ModeInstructions").GetComponent<Tab>());
        if (CHALLENGE_ENABLED) mode_tabs.Add(GameObject.Find("ModeChallenge").GetComponent<Tab>());
        if (QUIZ_ENABLED) mode_tabs.Add(GameObject.Find("ModeQuiz").GetComponent<Tab>());

        questions = new List<string>();
        options = new List<string>();
        answers = new List<int>(); //the correct answer
        givens = new List<int>(); //the recorded answer given by the user (default -1)

        questions.Add("Here is an example question- in what region is the water?");
        options.Add("A. Solid");
        options.Add("B. Liquid");
        options.Add("C. Vapor");
        options.Add("D. Two Phase");
        answers.Add(2);
        givens.Add(-1);

        questions.Add("Here is an example question- in what region is the water?");
        options.Add("A. Solid");
        options.Add("B. Liquid");
        options.Add("C. Vapor");
        options.Add("D. Two Phase");
        answers.Add(2);
        givens.Add(-1);

        questions.Add("Here is an example question- in what region is the water?");
        options.Add("A. Solid");
        options.Add("B. Liquid");
        options.Add("C. Vapor");
        options.Add("D. Two Phase");
        answers.Add(2);
        givens.Add(-1);

        questions.Add("Here is an example question- in what region is the water?");
        options.Add("A. Solid");
        options.Add("B. Liquid");
        options.Add("C. Vapor");
        options.Add("D. Two Phase");
        answers.Add(2);
        givens.Add(-1);

        questions.Add("Here is an example question- in what region is the water?");
        options.Add("A. Solid");
        options.Add("B. Liquid");
        options.Add("C. Vapor");
        options.Add("D. Two Phase");
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

    void SetQuizText() {
        qtext_tmp.SetText(questions[question]);
        for (int i = 0; i < 4; i++) option_tabs[i].tmp.SetText(options[question * 4 + i]);
    }

    void SetChallengeBall() {
        double volume = ThermoMath.v_given_percent(Random.Range(0.1f, 0.9f));
        double temperature = ThermoMath.t_given_percent(Random.Range(0.1f, 0.9f));
        double pressure = ThermoMath.p_given_vt(volume, temperature);
        challenge_dot.transform.localPosition = thermo_present.plot(pressure, volume, temperature);
    }

    void SetAllHalfed(bool h) {
        halfed = h;
        for (int i = 0; i < halfables.Count; i++)
            halfables[i].setHalf(halfed);
        //special case, only halfed when engaged
        if (!tool_coil.engaged) tool_coil.gameObject.GetComponent<Halfable>().setHalf(false);
        if (!tool_insulator.engaged) tool_insulator.gameObject.GetComponent<Halfable>().setHalf(false);
    }

    Vector3 popVector() {
        return new Vector3(Random.Range(-1f, 1f), 1f, Random.Range(-1f, 1f));
    }

    // The three functions below are used to manage attach/detach and storage of tools.
    // Generally, they have to set the transforms properly, update state variables,
    // and update text.
    void ActivateTool(Tool t) {
        GameObject o = t.gameObject;
        o.transform.SetParent(t.active.transform);
        t.touchable.grabbed = false;
        t.engaged = true;
        t.stored = false;
        t.boxcollider.isTrigger = true;
        // Not sure the two below are ever used? We don't have dials for these tools in use.
        //     if(t == tool_insulator) t.dial_dial.val = (float)ThermoMath.percent_given_t(thermo.temperature);
        //else if(t == tool_clamp)     t.dial_dial.val = (float)ThermoMath.percent_given_v(thermo.volume);
        t.dial_dial.Reset(); // reset tool when we add it.
        UpdateToolText(t);
        UpdateApplyTool(t);
        o.transform.localPosition = new Vector3(0f, 0f, 0f);
        o.transform.localRotation = Quaternion.identity;
        o.transform.localScale = new Vector3(1f, 1f, 1f);
        float v = t.active.transform.localScale.x; //can grab any dimension
        Vector3 invscale = new Vector3(1f / v, 1f / v, 1f / v);
        t.text.transform.localScale = invscale;
        Halfable h = o.GetComponent<Halfable>();
        if (h != null) h.setHalf(halfed); //conform to half-ness while engaged
    }
    void StoreTool(Tool t) {
        GameObject o = t.gameObject;
        o.transform.SetParent(t.storage.transform);
        t.touchable.grabbed = false;
        t.engaged = false;
        t.stored = true;
        o.transform.localPosition = new Vector3(0f, 0f, 0f);
        o.transform.localRotation = Quaternion.identity;
        o.transform.localScale = new Vector3(1f, 1f, 1f);
        float v = t.storage.transform.localScale.x; //can grab any dimension
        Vector3 invscale = new Vector3(1f / v, 1f / v, 1f / v);
        t.text.transform.localScale = invscale;
        Halfable h = o.GetComponent<Halfable>();
        if (h != null) h.setHalf(false); //Un-half when we store a tool.
        t.dial_dial.Reset(); // definitely need to reset tool when we store it.
        UpdateToolText(t);
        UpdateApplyTool(t);
    }
    void DetachTool(Tool t, Vector3 vel) {
        GameObject o = t.gameObject;
        o.transform.SetParent(t.touchable.og_parent);
        t.touchable.grabbed = false;
        t.engaged = false;
        t.stored = false;
        o.transform.localScale = new Vector3(1f, 1f, 1f);
        t.text.transform.localScale = new Vector3(1f, 1f, 1f);
        t.rigidbody.isKinematic = false;
        t.rigidbody.velocity = vel;
        t.dial_dial.Reset(); // may as well reset tool when we remove it, too.
        UpdateToolText(t);
        UpdateApplyTool(t);
    }

    /*
    tried during:
    - newly grabbed
    - newly snapped on
    - dial altered
    */
    void UpdateApplyTool(Tool t) //alters "applied_x"
    {
        if (t == tool_insulator) {
            //do nothing
            return;
        }
        else if (t == tool_burner || t == tool_coil) {
            if (tool_burner.engaged && tool_coil.engaged) {
                if (t == tool_burner) DetachTool(tool_coil, popVector());
                if (t == tool_coil) DetachTool(tool_burner, popVector());
            }

            if (t == tool_burner) {
                if (t.dial_dial.map == 0.0f) {
                    var e = flame.emission;
                    e.enabled = false;
                }
                else {
                    var e = flame.emission;
                    e.enabled = true;
                }
                var vel = flame.velocityOverLifetime;
                vel.speedModifierMultiplier = Mathf.Lerp(0.1f, 0.5f, t.dial_dial.val);
            }
            else if (t == tool_coil) {
                //TODO: coil visuals?
            }

            applied_heat = 0;
            if (tool_burner.engaged) applied_heat += tool_burner.dial_dial.map;
            if (tool_coil.engaged) applied_heat += tool_coil.dial_dial.map;
        }
        else if (t == tool_clamp) {
            //const float MAX_WIDTH = 2.0f;
            //const float TOP_BASE_Y_POS = 0.049f;
            //const float MID_BASE_Y_POS = 0.0164f;
            //var clamp_mid = GameObject.Find("clamp_mid");
            //var clamp_mid_mesh = clamp_mid.GetComponent<MeshFilter>().mesh;
            //var clamp_top = GameObject.Find("clamp_mid");
            //var clamp_top_mesh = clamp_top.GetComponent<MeshFilter>().mesh;
            //clamp_mid.transform.localScale = new Vector3(1.0f, 1.0f + MAX_WIDTH * t.dial_dial.val, 1.0f);
            //*clamp_mid.transform.localScale.y
            //clamp_mid.transform.position = new Vector3(clamp_mid.transform.position.x, MID_BASE_Y_POS + clamp_mid_mesh.bounds.extents.y, clamp_mid.transform.position.z);
            //clamp_top.transform.position = new Vector3(clamp_mid.transform.position.x, TOP_BASE_Y_POS + clamp_mid_mesh.bounds.extents.y, clamp_mid.transform.position.z);
            applied_weight = 0;
            if (tool_weight.engaged) applied_weight += tool_weight.dial_dial.map;
            if (tool_balloon.engaged) applied_weight += tool_balloon.dial_dial.map;
        }
        else if (t == tool_weight || t == tool_balloon) {
            if (tool_weight.engaged && tool_balloon.engaged) {
                if (t == tool_weight) DetachTool(tool_balloon, popVector());
                else if (t == tool_balloon) DetachTool(tool_weight, popVector());
            }

            float v = 1f;
            if (t == tool_weight) v += tool_weight.dial_dial.val;
            else if (t == tool_balloon) v += tool_balloon.dial_dial.val;
            Vector3 scale;
            Vector3 invscale;

            scale = new Vector3(v, v, v);
            invscale = new Vector3(1f / v, 1f / v, 1f / v);
            t.active.transform.localScale = scale;
            if (t.engaged) t.text.transform.localScale = invscale;

            v *= t.default_storage_scale;
            scale = new Vector3(v, v, v);
            invscale = new Vector3(1f / v, 1f / v, 1f / v);
            t.storage.transform.localScale = scale;
            if (t.stored) t.text.transform.localScale = invscale;

            //math
            applied_weight = 0;
            if (tool_weight.engaged) applied_weight += tool_weight.dial_dial.map;
            if (tool_balloon.engaged) applied_weight += tool_balloon.dial_dial.map;
        }

    }
    void UpdateToolText(Tool t) {
        if (!t.dial_dial) return;
        if (t == tool_insulator || t == tool_clamp) return;

        string txt = string.Format("{0:3} " + t.dial_dial.display_unit, (float)(t.dial_dial.map * t.dial_dial.display_mul));
        //t.textv_tmpro.SetText(txt);
        //if(t.textd_tmpro) t.textd_tmpro.SetText(txt);
        t.textv_tmpro.SetText("{0:3} " + t.dial_dial.display_unit, (float)(t.dial_dial.map * t.dial_dial.display_mul));
        if (t.textd_tmpro) t.textd_tmpro.SetText("{0:3} " + t.dial_dial.display_unit, (float)(t.dial_dial.map * t.dial_dial.display_mul));
    }

    //safe to call if not interactable, as it will just do nothing
    bool floatNumeric(float f) {
        if (double.IsNaN(f)) return false;
        if (double.IsInfinity(f)) return false;
        return true;
    }
    void TryInteractable(GameObject actable, Vector3 pos, ref float r_x, ref float r_y) {
        //grabbing handle
        if (actable == handle_workspace) {
            float dy = (r_y - pos.y);
            workspace.transform.position = new Vector3(workspace.transform.position.x, workspace.transform.position.y - dy, workspace.transform.position.z);
        }
        else if (actable == graph) {
            Vector3 localspace = graph.transform.InverseTransformPoint(pos);
            Vector3 correctedspace = new Vector3(localspace.z, localspace.y, -localspace.x) * 4.0f; //rotate 90, mul by 4 (inverse transform of gmodel)

            //note: thermospace is v,p,t

            //Vector3 thermoguess = thermo.guessPlot(ThermoMath.t_neutral, correctedspace.y, correctedspace.x);
            Vector3 thermoguess = thermo_present.guessMeshPlot(correctedspace.x, correctedspace.y, correctedspace.z);
            Vector3 localguess = thermo_present.plot(thermoguess.y, thermoguess.x, thermoguess.z); //note swizzle!

            if (floatNumeric(localguess.x) && floatNumeric(localguess.y) && floatNumeric(localguess.z)) {
                placement_dot.transform.localPosition = localguess;
                placement_thermo = thermoguess;
                placement_thermo_reasonable = true;
            }
            else placement_thermo_reasonable = false;
        }
        else {
            Dial d = actable.GetComponent<Dial>();
            //grabbing dial
            if (d != null) {
                Tool t = d.tool.GetComponent<Tool>();
                float dx = (r_x - pos.x) * -10f;
                d.val = Mathf.Clamp(d.val - dx, 0f, 1f);
                //if this close to either end, assume user wants min/max
                if (d.val < 0.005) d.val = 0f;
                if (d.val > 0.995) d.val = 1f;
                d.forceMap();

                UpdateApplyTool(t);
            }
        }
    }

    /*
     * This function seems to handle all possible interactions between the hand and other objects.
     * Honestly, I haven't quite got a full understanding of this ~200-line behemoth.
     */
    //"left_hand": true -> left, false -> right
    void TryHand(bool left_hand, float htrigger_val, float itrigger_val, Vector3 hand_pos, Vector3 hand_vel, ref bool ref_htrigger, ref bool ref_itrigger, ref int ref_htrigger_delta, ref int ref_itrigger_delta, ref float ref_x, ref float ref_y, ref GameObject ref_hand, ref GameObject ref_grabbed, ref GameObject ref_ohand, ref GameObject ref_ograbbed) {
        float htrigger_threshhold = 0.1f;
        float itrigger_threshhold = 0.1f;

        //find deltas
        ref_htrigger_delta = 0;
        if (!ref_htrigger && htrigger_val > htrigger_threshhold) {
            ref_htrigger_delta = 1;
            ref_htrigger = true;
        }
        else if (ref_htrigger && htrigger_val <= htrigger_threshhold) {
            ref_htrigger_delta = -1;
            ref_htrigger = false;
        }

        ref_itrigger_delta = 0;
        if (!ref_itrigger && itrigger_val > itrigger_threshhold) {
            ref_itrigger_delta = 1;
            ref_itrigger = true;
            ref_x = hand_pos.x;
            ref_y = hand_pos.y;
        }
        else if (ref_itrigger && itrigger_val <= itrigger_threshhold) {
            ref_itrigger_delta = -1;
            ref_itrigger = false;
        }

        //find new grabs
        if (ref_grabbed == null && ref_htrigger_delta == 1) {
            //first try movables
            for (int i = 0; ref_grabbed == null && i < movables.Count; i++) {
                if ( //object newly grabbed
                   (left_hand && movables[i].ltouch) ||
                   (!left_hand && movables[i].rtouch)
                  ) {
                    ref_grabbed = movables[i].gameObject;
                    ref_grabbed.transform.SetParent(ref_hand.transform);
                    if (ref_grabbed == ref_ograbbed) ref_ograbbed = null;
                    movables[i].grabbed = true;
                    Tool t = ref_grabbed.GetComponent<Tool>();
                    if (t) //newly grabbed object is a tool
                    {
                        t.audioS.Play();
                        t.engaged = false;
                        t.stored = false;
                        ref_grabbed.transform.localScale = new Vector3(1f, 1f, 1f);
                        t.text.transform.localScale = new Vector3(1f, 1f, 1f);
                        t.rigidbody.isKinematic = true;
                        t.boxcollider.isTrigger = false;
                        UpdateApplyTool(t);
                    }
                    VisAid v = ref_grabbed.GetComponent<VisAid>();
                    if (v) //newly grabbed object is a visaid
                    {
                        v.stored = false;
                        v.rigidbody.isKinematic = true;
                    }
                }
            }
            //then dials
            if (ref_grabbed == null) {
                for (int i = 0; i < tools.Count; i++) {
                    if ( //dial newly grabbed
                       (left_hand && tools[i].dial_touchable.ltouch) ||
                       (!left_hand && tools[i].dial_touchable.rtouch)
                      ) {
                        ref_grabbed = tools[i].dial;
                        tools[i].dial_touchable.grabbed = true;
                        if (ref_grabbed == ref_ograbbed) ref_ograbbed = null;
                    }
                }
            }

            //then extraaneous
            if (ref_grabbed == null) //still not holding anything
            {
                Touchable g = handle_workspace_touchable;
                if ( //handle newly grabbed
                  (left_hand && g.ltouch) ||
                  (!left_hand && g.rtouch)
                ) {
                    ref_grabbed = handle_workspace;
                    g.grabbed = true;
                    if (ref_grabbed == ref_ograbbed) ref_ograbbed = null;
                }
            }

            if (ref_grabbed == null) //still not holding anything
            {
                Touchable g = graph_touchable;
                if ( //graph newly grabbed
                  (left_hand && g.ltouch) ||
                  (!left_hand && g.rtouch)
                ) {
                    ref_grabbed = graph;
                    state_dot.GetComponent<Renderer>().enabled = false;
                    placement_dot.GetComponent<Renderer>().enabled = true;
                    g.grabbed = true;
                    if (ref_grabbed == ref_ograbbed) ref_ograbbed = null;
                }
            }

            if (ref_grabbed != null) //something newly grabbed
            {
                Halfable h = ref_grabbed.GetComponent<Halfable>();
                if (h != null) h.setHalf(false); //nothing should be halfed while being grabbed
            }
        }
        //find new releases
        else if (ref_grabbed && ref_htrigger_delta == -1) //something newly released
        {
            Tool t = ref_grabbed.GetComponent<Tool>();
            if (t) //tool newly released
            {
                t.audioS.Play();
                if (t.active_ghost.tintersect) //tool released making it active
                {
                    ActivateTool(t);
                }
                else if (t.storage_ghost.tintersect) //tool released making it stored
                {
                    StoreTool(t);
                }
                else //tool released nowhere special
                {
                    DetachTool(t, hand_vel);
                }
            }
            else //newly released object is NOT a tool
            {
                ref_grabbed.transform.SetParent(ref_grabbed.GetComponent<Touchable>().og_parent); //ok to do, even with a dial
                VisAid v = ref_grabbed.GetComponent<VisAid>();
                if (v) //visaid newly released
                {
                    v.rigidbody.isKinematic = false;
                    v.rigidbody.velocity = hand_vel;
                }
                if (ref_grabbed == graph) {
                    placement_dot.GetComponent<Renderer>().enabled = false;
                    if (placement_thermo_reasonable) thermo_present.warp_pv(placement_thermo.y, placement_thermo.x, placement_thermo.z);
                    state_dot.GetComponent<Renderer>().enabled = true;
                }
            }

            ref_grabbed.GetComponent<Touchable>().grabbed = false;
            ref_grabbed = null;
        }

        if (ref_grabbed) TryInteractable(ref_grabbed, hand_pos, ref ref_x, ref ref_y);

        ref_x = hand_pos.x;
        ref_y = hand_pos.y;

        if (ref_grabbed == null) //still not holding anything
        {
            if (halfer_fingertoggleable.finger) //finger hitting halfer object
            {
                if ( //we're currently checking the correct hand
                  (left_hand && halfer_fingertoggleable.lfinger) ||
                  (!left_hand && halfer_fingertoggleable.rfinger)
                ) { //halfer state
                    if (halfer_fingertoggleable.on) {
                        halfer.GetComponent<AudioSource>().Play();
                        SetAllHalfed(!halfed);
                        halfer_fingertoggleable.on = false;
                    }
                }
            }
            else halfer_fingertoggleable.on = false;

            if (reset_fingertoggleable.finger) //finger hitting reset object
            {
                if ( //we're currently checking the correct hand
                  (left_hand && reset_fingertoggleable.lfinger) ||
                  (!left_hand && reset_fingertoggleable.rfinger)
                ) { //reset state
                    if (reset_fingertoggleable.on) {
                        reset.GetComponent<AudioSource>().Play();
                        for (int i = 0; i < tools.Count; i++)
                            if (tools[i].engaged) DetachTool(tools[i], new Vector3(0.0f, 0.0f, 0.0f));
                        thermo_present.Reset();
                        reset_fingertoggleable.on = false;
                    }
                }
            }
            else reset_fingertoggleable.on = false;
        }

        //centerer
        if (vrcenter_fingertoggleable.finger) //finger hitting vrcenter object
        {
            if ( //we're currently checking the correct hand
              (left_hand && vrcenter_fingertoggleable.lfinger) ||
              (!left_hand && vrcenter_fingertoggleable.rfinger)
            ) { //reset center position
                vrcenter_backing_meshrenderer.material = tab_hisel;
                UnityEngine.XR.InputTracking.Recenter();
                OVRManager.display.RecenterPose();
                Vector3 pos = cam_offset.transform.localPosition - (cam_offset.transform.localPosition + ceye.transform.localPosition);
                pos.y = 0f;
                cam_offset.transform.localPosition = pos;
            }
            else vrcenter_backing_meshrenderer.material = tab_hi;
        }
        else vrcenter_backing_meshrenderer.material = tab_default;

    }

    /*                 
     * Function to update object materials/appearance in response to a "grab" event.
     */
    void UpdateGrabVis() {
        for (int i = 0; i < tools.Count; i++) {
            Tool t = tools[i];
            GameObject g = t.gameObject;

            if (t.storage == null) {
                continue; // tool that does not have a physical mesh
            }

            if (lgrabbed == g || rgrabbed == g) {
                //active
                if (t.active_ghost.tintersect) {
                    t.active_available_meshrenderer.enabled = true;
                    t.active_available_meshrenderer.material = GameDB.Instance.SnapMat;
                }
                else {
                    t.active_available_meshrenderer.enabled = true;
                    t.active_available_meshrenderer.material = GameDB.Instance.AvailableMat;
                }
                //storage
                if (t.storage_ghost.tintersect) {
                    t.storage_meshrenderer.enabled = true;
                    t.storage_meshrenderer.material = GameDB.Instance.SnapMat;
                }
                else {
                    t.storage_meshrenderer.enabled = true;
                    t.storage_meshrenderer.material = GameDB.Instance.AvailableMat;
                }
            }
            else {
                t.active_available_meshrenderer.enabled = false;
                t.storage_meshrenderer.enabled = false;
            }
        }

        Touchable gr;
        bool ltouch = false;
        bool rtouch = false;
        for (int i = 0; i < movables.Count; i++) {
            gr = movables[i];
            if (gr.ltouch) ltouch = true;
            if (gr.rtouch) rtouch = true;
        }
        for (int i = 0; i < tools.Count; i++) {
            gr = tools[i].dial_touchable;
            if (gr.ltouch) ltouch = true;
            if (gr.rtouch) rtouch = true;
        }
        gr = handle_workspace_touchable;
        if (gr.ltouch) ltouch = true;
        if (gr.rtouch) rtouch = true;
        gr = graph_touchable;
        if (gr.ltouch) ltouch = true;
        if (gr.rtouch) rtouch = true;
        Touchable press_gr = halfer_touchable;
        if (press_gr.ltouch) ltouch = true;
        if (press_gr.rtouch) rtouch = true;
        press_gr = reset_touchable;
        if (press_gr.ltouch) ltouch = true;
        if (press_gr.rtouch) rtouch = true;

        if (lgrabbed) lhand_meshrenderer.materials = hand_grabbings;
        else if (ltouch) lhand_meshrenderer.materials = hand_touchings;
        else lhand_meshrenderer.materials = hand_emptys;

        if (rgrabbed) rhand_meshrenderer.materials = hand_grabbings;
        else if (rtouch) rhand_meshrenderer.materials = hand_touchings;
        else rhand_meshrenderer.materials = hand_emptys;
    }

    //give it a list of fingertoggleables, and it manipulates them to act as a singularly-selectable list
    int reconcileDependentSelectables(int known, List<Tab> list) {
        int n_toggled = 0;
        Tab t;
        for (int i = 0; i < list.Count; i++) {
            t = list[i];
            if (t.fingertoggleable.on) n_toggled++;
        }

        if (n_toggled <= 1) {
            known = -1;
            for (int i = 0; i < list.Count; i++) {
                t = list[i];
                if (t.fingertoggleable.on) known = i;
            }
        }
        else //need conflict resolution!
        {
            known = -1;
            for (int i = 0; i < list.Count; i++) {
                t = list[i];
                if (t.fingertoggleable.on) {
                    if (known == -1) known = i;
                    else {
                        //if only t is intersecting, prefer t
                        if (t.fingertoggleable.finger && !list[known].fingertoggleable.finger) {
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

    void updateSelectableVis(int known, List<Tab> list) {
        Tab t;
        for (int i = 0; i < list.Count; i++) {
            t = list[i];
            if (known == i) t.backing_meshrenderer.material = tab_sel;
            else t.backing_meshrenderer.material = tab_default;
        }
    }

    /*
     * Another behemoth, does frame-by-frame updates to state, appearances, transforms, etc.
     * Includes calls to TryHand, UpdateGrabVis, etc. as well as calls to ThermoState functions.
     * Basically, wraps calls to a bunch of other functions, and a hodgepodge of other random tasks,
     * as far as I can tell.
     */
    void FixedUpdate() {
        //hands keep trying to run away- no idea why (this is a silly way to keep them still)
        lactualhand.transform.localPosition = new Vector3(0f, 0f, 0f);
        lactualhand.transform.localEulerAngles = new Vector3(0f, 0f, 90f);
        ractualhand.transform.localPosition = new Vector3(0f, 0f, 0f);
        ractualhand.transform.localEulerAngles = new Vector3(0f, 0f, -90f);

        double delta_time = (double)Time.fixedDeltaTime;

        //apply thermo
        if (!tool_clamp.engaged) {
            ambient_pressure = tool_ambientPressure.dial_dial.map;
            room_temp = tool_roomTemp.dial_dial.map;
            double psi_to_pascal = 6894.76;
            double neutral_pressure = 14.6959; //~1atm in psi
            double weight_pressure = (applied_weight) / thermo_present.get_surfacearea_insqr(); //psi
            weight_pressure += ambient_pressure + neutral_pressure; // TODO: establish logical bounds and units on the ambient pressure dial
            weight_pressure *= psi_to_pascal; //conversion from psi to pascal

            //treat "applied_weight" as target, and iterate toward it, rather than applying it additively
            //(technically, "should" add_pressure(...) with every delta of weight on the piston, but that would result in very jumpy nonsense movements. iterating toward a target smooths it out)
            double delta_pressure = (weight_pressure - thermo_present.get_pressure());

            bool significantChange = Mathf.Abs((float)delta_pressure * Time.fixedDeltaTime) > DELTA_PRESSURE_CUTOFF ? true : false;
            if (significantChange) {
                // we only are applying pressure added/released if we are in gas region.
                // if this ever changes, adapt this "if" check as needed.
                if (thermo_present.get_region() == 2) {
                    arrows.Go(delta_pressure > 0.0f);
                    arrows.SetFlow(delta_pressure);
                }
            }
            else if (arrows.running) {
                arrows.Stop();
            }

            if (System.Math.Abs(delta_pressure) > 0) {
                if (tool_insulator.engaged) thermo_present.add_pressure_insulated_per_delta_time(delta_pressure, delta_time); // Pressure Constrained -> Insulated ->  delta pressure
                else thermo_present.add_pressure_uninsulated_per_delta_time(delta_pressure, delta_time); // Pressure Constrained -> Uninsulated ->  delta pressure
            }
        }

        // stop arrows if there was a jump to another region other than gas
        if (thermo_present.get_region() != 2 && arrows.running) {
            arrows.Stop();
        }

        double insulation_coefficient = tool_insulator.engaged ? 1.0f : CONTAINER_INSULATION_COEFFICIENT;

        /*
            if (!tool_insulator.engaged) //heat leak
            {
              double room_temp = 295.372f;
              double diff = room_temp - thermo.temperature;
              heat_joules += diff * 100 * (double)Time.fixedDeltaTime; //nonsense calculation
              arrows.Go(diff < 0);
              arrows.SetFlow(diff / 100.0f);
            }
            else if (arrows.running) arrows.Stop();
        */
        if (applied_heat != 0) {
            if (tool_clamp.engaged) thermo_present.add_heat_constant_v_per_delta_time(applied_heat, insulation_coefficient, delta_time);    // Volume Constrained -> Insulated && Uninsulated ->  delta energy     // time eqtn 6a
            else thermo_present.add_heat_constant_p_per_delta_time(applied_heat, insulation_coefficient, delta_time);                       // Pressure Constrained -> Insulated && Uninsulated -> delta energy    // time eqtn 6b
        }

        // TODO: apply heat transfer from piston to outer room and vice versa

        //running blended average of hand velocity (transfers this velocity on "release object" for consistent "throwing")
        lhand_vel += (lhand.transform.position - lhand_pos) / Time.fixedDeltaTime;
        lhand_vel *= 0.5f;
        lhand_pos = lhand.transform.position;

        rhand_vel += (rhand.transform.position - rhand_pos) / Time.fixedDeltaTime;
        rhand_vel *= 0.5f;
        rhand_pos = rhand.transform.position;

        //input
        //float lhandt  = OVRInput.Get(OVRInput.Axis1D.PrimaryHandTrigger);
        //float lindext = OVRInput.Get(OVRInput.Axis1D.PrimaryIndexTrigger);
        //float rhandt  = OVRInput.Get(OVRInput.Axis1D.SecondaryHandTrigger);
        //float rindext = OVRInput.Get(OVRInput.Axis1D.SecondaryIndexTrigger);
        float lhandt = OVRInput.Get(OVRInput.RawAxis1D.LHandTrigger);
        float lindext = OVRInput.Get(OVRInput.RawAxis1D.LIndexTrigger);
        float rhandt = OVRInput.Get(OVRInput.RawAxis1D.RHandTrigger);
        float rindext = OVRInput.Get(OVRInput.RawAxis1D.RIndexTrigger);

        /*
            //index compatibility
            if (OVRInput.Get(OVRInput.Button.One, OVRInput.Controller.LTouch))
            {
              lhandt = 1.0f;
            }
            if(OVRInput.Get(OVRInput.Button.One,OVRInput.Controller.RTouch))
            {
              rhandt = 1.0f;
            }
            lhandt += lindext;
            rhandt += rindext;
            if (lindext > 0.0f || rindext > 0.0f)
            {
              ;
            }
        */

        //test effect of hands one at a time ("true" == "left hand", "false" == "right hand")
        TryHand(true, lhandt, lindext, lhand.transform.position, lhand_vel, ref lhtrigger, ref litrigger, ref lhtrigger_delta, ref litrigger_delta, ref lz, ref ly, ref lhand, ref lgrabbed, ref rhand, ref rgrabbed); //left hand
        TryHand(false, rhandt, rindext, rhand.transform.position, rhand_vel, ref rhtrigger, ref ritrigger, ref rhtrigger_delta, ref ritrigger_delta, ref rz, ref ry, ref rhand, ref rgrabbed, ref lhand, ref lgrabbed); //right hand

        UpdateGrabVis();

        //clipboard
        int old_board_mode = board_mode;
        board_mode = reconcileDependentSelectables(board_mode, mode_tabs);
        if (board_mode == -1) board_mode = old_board_mode;
        updateSelectableVis(board_mode, mode_tabs);

        if (!CHALLENGE_ENABLED && board_mode == 1) board_mode = 0;
        if (!QUIZ_ENABLED && board_mode == 2) board_mode = 0;
        switch (board_mode) {
            case 0: //instructions
                if (!instructions_parent.activeSelf) instructions_parent.SetActive(true);
                if (challenge_parent.activeSelf) { challenge_parent.SetActive(false); challenge_dot.SetActive(false); }
                if (quiz_parent.activeSelf) quiz_parent.SetActive(false);
                break;
            case 1: //challenge
                if (instructions_parent.activeSelf) instructions_parent.SetActive(false);
                if (!challenge_parent.activeSelf) { challenge_parent.SetActive(true); challenge_dot.SetActive(true); }
                if (quiz_parent.activeSelf) quiz_parent.SetActive(false);

                if (challenge_ball_collide.win) {
                    challenge_ball_collide.win = false;
                    SetChallengeBall();
                }
                break;
            case 2: //quiz
                if (instructions_parent.activeSelf) instructions_parent.SetActive(false);
                if (challenge_parent.activeSelf) { challenge_parent.SetActive(false); challenge_dot.SetActive(false); }
                if (!quiz_parent.activeSelf) quiz_parent.SetActive(true);
                qselected = reconcileDependentSelectables(qselected, option_tabs);
                updateSelectableVis(qselected, option_tabs);
                if (qselected != -1) //answer selected, can be confirmed
                {
                    if (qconfirm_tab.fingertoggleable.finger) qconfirm_tab.backing_meshrenderer.material = tab_hisel; //touching
                    else qconfirm_tab.backing_meshrenderer.material = tab_sel;   //not touching

                    if (!qconfirmed && qconfirm_tab.fingertoggleable.finger) //newly hit
                    {
                        qconfirmed = true;
                        if (qselected == answers[question])
                            ; //correct
                        else
                            ; //incorrect
                        givens[question] = qselected;
                        question++;
                        SetQuizText();
                    }
                }
                else //no answer selected- can't confirm
                {
                    if (qconfirm_tab.fingertoggleable.finger) qconfirm_tab.backing_meshrenderer.material = tab_hi;      //touching
                    else qconfirm_tab.backing_meshrenderer.material = tab_default; //not touching

                    qconfirmed = false;
                }
                break;
            case 3: //congratulations
                board_mode = 2; //immediately return back to quiz until this section is actually implemented
                break;
        }

        //tooltext
        Tool t;
        for (int i = 0; i < tools.Count; i++) {
            t = tools[i];
            t.dial_dial.examined = false;
            if (t.dial == lgrabbed || t.dial == rgrabbed) t.dial_dial.examined = true;
            if (t.dial_dial.val != t.dial_dial.prev_val) {
                UpdateToolText(t);
                t.dial_dial.examined = true;
            }
            t.dial_dial.prev_val = t.dial_dial.val;
        }

        // Disable the warning on that the tools are not available in this region
        /* 
        switch (thermo.region)
        {
          case 0:
          case 1:
            tool_weight.textn.GetComponent<MeshRenderer>().enabled = true;
            tool_weight.disabled = true;
            tool_balloon.textn.GetComponent<MeshRenderer>().enabled = true;
            tool_balloon.disabled = true;
            break;
          case 2:
            tool_weight.textn.GetComponent<MeshRenderer>().enabled = false;
            tool_weight.disabled = false;
            tool_balloon.textn.GetComponent<MeshRenderer>().enabled = false;
            tool_balloon.disabled = false;
            break;
        }
        */
        for (int i = 0; i < tools.Count; i++) {
            t = tools[i];
            if (t.text_fadable == null) {
                continue;
            }
            if (!t.text_fadable.stale) {
                if (t.text_fadable.alpha == 0f) {
                    t.textv_meshrenderer.enabled = false;
                    t.textl_meshrenderer.enabled = false;
                }
                else {
                    t.textv_meshrenderer.enabled = true;
                    t.textl_meshrenderer.enabled = true;
                    Color32 c = t.disabled ? new Color32(70, 70, 70, (byte)(t.text_fadable.alpha * 255)) : new Color32(0, 0, 0, (byte)(t.text_fadable.alpha * 255));
                    t.textv_tmpro.faceColor = c;
                    t.textl_tmpro.faceColor = c;
                }
            }
        }
        thermo_present.UpdateErrorState();
    }

}

