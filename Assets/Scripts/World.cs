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
using ThermoVR.Controls;
using ThermoVR.Dials;
using BeauUtil;
using ThermoVR.Tools;
using UnityEngine.UIElements;
using ThermoVR;
using static OVRInput;
using ThermoVR.UI.GraphElements;
using ThermoVR.UI;
using ThermoVR.Lab;
using ThermoVR.State;
using System;

public class World : MonoBehaviour
{
    public static World Instance;

    #region Consts

    const float CONTAINER_INSULATION_COEFFICIENT = 0.1f; // 0.1f; // Not really based on a physical material, just a way to roughly simulate imperfect insulation.
    public const double DELTA_PRESSURE_CUTOFF = 100.0;
    const double PSI_TO_PASCAL = 6894.76;
    const double SPECIFIC_HEAT_CAPACITY_LIQ = 4184; // how many J it takes to heat 1 kg of water liquid 1 Kelvin
    const double SPECIFIC_HEAT_CAPACITY_VAP = 1.996; // how many J it takes to heat 1 kg of water vapor 1 Kelvin
    public const float BURNER_MAX = 100000;
    public const float COIL_MAX = -100000;

    #endregion // Consts

    #region Inspector

    public WorldModMgr ModMgr;

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
    [Space(5)]
    [Header("Thermo")]
    [SerializeField] private ThermoPresent thermo_present;
    [SerializeField] private Tablet tablet;
    // TODO: assign these
    [SerializeField] private GraphElement[] graph_elements; // visual aid overlays on graph (region labels, number lines, etc.)

    [Space(5)]
    [Header("Controls")]
    [SerializeField] private GameObject cam_offset;
    [SerializeField] private GameObject ceye;
    [SerializeField] private ControllerAnchor lhand;
    [SerializeField] private ControllerAnchor rhand;
    [SerializeField] private GameObject origin;

    /*
    GameObject vrcenter;
    FingerToggleable vrcenter_fingertoggleable;
    MeshRenderer vrcenter_backing_meshrenderer;
    */

    //[Space(5)]
    //[Header("Moveables")]
    List<Touchable> movables;
    GameObject workspace;
    [SerializeField] private GameObject handle_workspace;
    Touchable handle_workspace_touchable;
    Touchable graph_touchable;
    GameObject lgrabbed = null;
    GameObject rgrabbed = null;
    int lhtrigger_delta = 0;
    int litrigger_delta = 0;
    int rhtrigger_delta = 0;
    int ritrigger_delta = 0;
    Vector3 lpos = new Vector3(0f, 0f, 0f);
    Vector3 rpos = new Vector3(0f, 0f, 0f);

    [Space(5)]
    [Header("Tools")]
    [SerializeField] private Tool tool_insulator;
    //[SerializeField] private Tool tool_clamp;
    [SerializeField] private Tool tool_stop1;
    [SerializeField] private Tool tool_stop2;
    [SerializeField] private Tool tool_burner;
    [SerializeField] private Tool tool_coil;
    [SerializeField] private Tool tool_weight;
    [SerializeField] private Tool tool_balloon;
    [SerializeField] private Tool tool_ambientPressure;
    [SerializeField] private Tool tool_roomTemp;
    // [SerializeField] private Tool tool_percentInsulation;
    // [SerializeField] private Tool tool_clampRange;

    List<Tool> tools;

    [Space(5)]
    [Header("Dials")]
    // [SerializeField] private Dial dial_insulator;
    [SerializeField] private Dial dial_stop1;
    [SerializeField] private Dial dial_stop2;
    [SerializeField] private Dial dial_burner;
    [SerializeField] private Dial dial_coil;
    [SerializeField] private Dial dial_weight;
    [SerializeField] private Dial dial_balloon;
    [SerializeField] private Dial dial_ambientPressure;
    [SerializeField] private Dial dial_roomTemp;
    [SerializeField] private Dial dial_percentInsulation;

    //[SerializeField] private ThermoToggle toggle_heatTransfer;
    [SerializeField] private PhysicalToggle toggle_heatTransfer;

    public List<Dial> dials;
    ParticleSystem flame; //special case

    [Space(5)]
    [Header("Halfables")]
    bool halfed = false;
    List<Halfable> halfables;
    [SerializeField] Pressable reset_button;
    [SerializeField] Pressable halfer_button;

    [Space(5)]
    [Header("Dot Placement")]
    GameObject vessel;
    GameObject graph;
    GameObject state_dot;
    GameObject placement_dot;
    Vector3 placement_thermo;
    bool placement_thermo_reasonable;

    bool lhtrigger = false;
    bool rhtrigger = false;
    bool litrigger = false;
    bool ritrigger = false;

    // sim variables
    double room_temp = 292; // in K
    double applied_heat = 0;
    double applied_weight = 0;
    double ambient_pressure = 0;

    private List<Pressable> m_pressables; // pressables register themselves with this on event

    #endregion // Inspector

    #region Initialization

    private void Awake() {
        if (Instance == null) {
            Instance = this;
        }
        else if (this != Instance) {
            Destroy(this.gameObject);
            return;
        }

        GameMgr.Events?.Register<Pressable>(GameEvents.RegisterPressable, HandleRegisterPressable);
        GameMgr.Events?.Register<Touchable>(GameEvents.RegisterMovable, HandleRegisterMovable);

        GameMgr.Events?.Register(GameEvents.WarpPVT, HandleWarpPVT);

        GameMgr.Events?.Register<Tool>(GameEvents.PressedToolToggle, HandleToolTogglePressed);

        GameMgr.Events?.Register<List<ToolType>>(GameEvents.UpdateAllowedTools, HandleAllowedToolsUpdated);
        GameMgr.Events?.Register(GameEvents.ResetToolRestrictions, HandleResetToolRestrictions);

        movables = new List<Touchable>();
    }


    public void Init() {
        // All this code does, at end of day, is find all the objects to manage,
        // and set initial values and such as needed.
        // thermo_present = GameObject.Find("Oracle").GetComponent<ThermoPresent>();

        hand_emptys = new Material[] { hand_empty };
        hand_touchings = new Material[] { hand_touching };
        hand_grabbings = new Material[] { hand_grabbing };

        lhand.Init(hand_emptys);
        rhand.Init(hand_emptys);

        double kg_corresponding_to_10mpa = thermo_present.get_surfacearea_insqr() * (10 * 1453.8/*MPa->psi*/) * 0.453592/*lb->kg*/;
        double kg_corresponding_to_2mpa = thermo_present.get_surfacearea_insqr() * (2 * 1453.8/*MPa->psi*/) * 0.453592/*lb->kg*/; // 10 MPa seems way too big, sooooo... we'll just do 2 MPa.

        // As we grab them, set ranges on tool dials (sliders).
        tools = new List<Tool> {
            tool_insulator,
            tool_stop1,
            tool_stop2,
            tool_burner,
            tool_coil,
            tool_weight,
            tool_balloon,
            tool_ambientPressure,
            tool_roomTemp,
            // tool_percentInsulation
        };

        tool_insulator.Init(Units.Quality);
        tool_stop1.Init(Units.Volume);
        tool_stop2.Init(Units.Volume);
        tool_burner.Init(Units.Heat, 0.001f);
        tool_coil.Init(Units.Heat, 0.001f);
        tool_weight.Init(Units.Weight);
        tool_balloon.Init(Units.Weight);
        tool_ambientPressure.Init(Units.AmbientPressure, 0.001f); // display in kPa
        tool_roomTemp.Init(Units.TemperatureK);
        // tool_percentInsulation.Init(Units.Percent);

        dials = new List<Dial> {
            // dial_insulator,
            dial_stop1,
            dial_stop2,
            dial_burner,
            dial_coil,
            dial_weight,
            dial_balloon,
            dial_ambientPressure,
            dial_roomTemp,
            dial_percentInsulation
        };

        // dial_insulator.Init(0f, 1f);
        dial_stop1.Init((float)ThermoMath.v_min, (float)ThermoMath.v_max, DigitFormat.Volume);
        dial_stop2.Init((float)ThermoMath.v_min, (float)ThermoMath.v_max, DigitFormat.Volume);
        dial_burner.Init(0f, BURNER_MAX, DigitFormat.Heat);
        dial_coil.Init(0f, COIL_MAX, DigitFormat.Heat);
        dial_weight.Init(0f, (float)kg_corresponding_to_10mpa / 5.0f, DigitFormat.Weight);
        dial_balloon.Init(0f, -(float)kg_corresponding_to_10mpa / 5.0f, DigitFormat.Weight); // 500.0f
        dial_ambientPressure.Init((float)ThermoMath.p_min, (float)ThermoMath.p_max, DigitFormat.AmbientPressure);
        dial_roomTemp.Init(273, 366, DigitFormat.TemperatureK); // -100 to 200 fahrenheit // default val of 0.55 sets to 292 kelvin (72 degrees fahrenheit)
        dial_percentInsulation.Init(0f, 100, DigitFormat.Percent);

        room_temp = tool_roomTemp.get_val(); // in K

        // Gather pressables
        m_pressables = new List<Pressable>();
        GameMgr.Events.Dispatch(GameEvents.GatherPressables);

        // Initialize Buttons
        reset_button.OnPress += HandleResetPressed;
        halfer_button.OnPress += HandleHalferPressed;

        // toggle_heatTransfer.Init();
        // toggle_heatTransfer.Pressable.PressCompleted += HandleHeatTransferToggle;

        // Initialize Tablet (and corresponding buttons)
        tablet.Init();

        flame = GameObject.Find("Flame").GetComponent<ParticleSystem>();

        workspace = GameObject.Find("Workspace");
        handle_workspace_touchable = handle_workspace.GetComponent<Touchable>();

        // set initial states of meshrenderers and transforms for our tools.
        for (int i = 0; i < tools.Count; i++) {
            Tool t = tools[i];
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

            GameMgr.Events?.Dispatch(GameEvents.UpdateToolText, t);
        }

        vessel = GameObject.Find("Vessel");
        graph = GameObject.Find("Graph");
        graph_touchable = graph.GetComponent<Touchable>();
        state_dot = GameObject.Find("gstate");
        placement_dot = GameObject.Find("tstate");
        placement_dot.GetComponent<Renderer>().enabled = false;
        placement_thermo_reasonable = false;

        for (int i = 0; i < tools.Count; i++) movables.Add(tools[i].touchable); //important that tools take priority, so they can be grabbed and removed
        movables.Add(tablet.touchable);

        halfables = new List<Halfable> {
            GameObject.Find("Container").GetComponent<Halfable>(),
            GameObject.Find("Tool_Insulator Variant").GetComponent<Halfable>(),
            GameObject.Find("Tool_Coil Variant").GetComponent<Halfable>()
        };

        SetAllHalfed(true);
    }

    #endregion // Initialization

    #region Callbacks

    /*
     * Another behemoth, does frame-by-frame updates to state, appearances, transforms, etc.
     * Includes calls to TryHand, UpdateGrabVis, etc. as well as calls to ThermoState functions.
     * Basically, wraps calls to a bunch of other functions, and a hodgepodge of other random tasks,
     * as far as I can tell.
     */
    public void ManualFixedUpdate() {
        arrows.ManualFixedUpdate(); // reset arrows for next instruction

        StabilizeHands();

        ApplyTools();

        EnforceLimits();

        ProcessInputs();

        UpdateGrabVis();

        UpdateToolTexts();

        ProcessErrors();
    }

    #endregion // Callbacks

    #region Tools

    // The three functions below are used to manage attach/detach and storage of tools.
    // Generally, they have to set the transforms properly, update state variables,
    // and update text.
    void ActivateTool(Tool t) {
        GameObject o = t.gameObject;
        o.transform.SetParent(t.active.transform);
        t.touchable.grabbed = false;
        t.engaged = true;
        if (t == tool_stop1) {
            thermo_present.add_v_stop(tool_stop1.get_val(), t);
        }
        else if (t == tool_stop2) {
            thermo_present.add_v_stop(tool_stop2.get_val(), t);
        }
        t.stored = false;
        t.boxcollider.isTrigger = true;
        // Not sure the two below are ever used? We don't have dials for these tools in use.
        //     if(t == tool_insulator) t.dial_dial.val = (float)ThermoMath.percent_given_t(thermo.temperature);
        //else if(t == tool_clamp)     t.dial_dial.val = (float)ThermoMath.percent_given_v(thermo.volume);
        GameMgr.Events?.Dispatch(GameEvents.ActivateTool, t);
        GameMgr.Events?.Dispatch(GameEvents.UpdateToolText, t);
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
        if (t == tool_stop1) {
            thermo_present.release_v_stop(t);
        }
        else if (t == tool_stop2) {
            thermo_present.release_v_stop(t);
        }
        t.stored = true;
        o.transform.localPosition = new Vector3(0f, 0f, 0f);
        o.transform.localRotation = Quaternion.identity;
        o.transform.localScale = new Vector3(1f, 1f, 1f);
        float v = t.storage.transform.localScale.x; //can grab any dimension
        Vector3 invscale = new Vector3(1f / v, 1f / v, 1f / v);
        t.text.transform.localScale = invscale;
        Halfable h = o.GetComponent<Halfable>();
        if (h != null) h.setHalf(false); //Un-half when we store a tool.
        GameMgr.Events?.Dispatch(GameEvents.StoreTool, t);
        GameMgr.Events?.Dispatch(GameEvents.UpdateToolText, t);
        UpdateApplyTool(t);
    }
    void DetachTool(Tool t, Vector3 vel) {
        GameObject o = t.gameObject;
        o.transform.SetParent(t.touchable.og_parent);
        t.touchable.grabbed = false;
        if (!t.always_engaged) {
            t.engaged = false;
        }
        if (t == tool_stop1) {
            thermo_present.release_v_stop(t);
        }
        else if (t == tool_stop2) {
            thermo_present.release_v_stop(t);
        }
        t.stored = false;
        o.transform.localScale = new Vector3(1f, 1f, 1f);
        t.text.transform.localScale = new Vector3(1f, 1f, 1f);
        t.rigidbody.isKinematic = false;
        t.rigidbody.velocity = vel;
        GameMgr.Events?.Dispatch(GameEvents.DetachTool, t);
        GameMgr.Events?.Dispatch(GameEvents.UpdateToolText, t);
        UpdateApplyTool(t);
    }

    void AllowTool(Tool t) {
        // TODO: this
        t.gameObject.SetActive(true);
    }

    void DisallowTool(Tool t) {
        // TODO: this
        DetachTool(t, Vector3.zero);
        t.gameObject.SetActive(false);
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
                if (tool_burner.get_val() == 0.0f) {
                    var e = flame.emission;
                    e.enabled = false;
                }
                else {
                    var e = flame.emission;
                    e.enabled = true;
                }
                var vel = flame.velocityOverLifetime;
                vel.speedModifierMultiplier = Mathf.Lerp(0.1f, 0.5f, t.get_val());
            }
            else if (t == tool_coil) {
                //TODO: coil visuals?
            }

            applied_heat = 0;
            if (tool_burner.engaged) applied_heat += tool_burner.get_val();
            if (tool_coil.engaged) applied_heat += tool_coil.get_val();
        }
        else if (t == tool_stop1 || t == tool_stop2) {
            if (t == tool_stop1) {
                thermo_present.update_v_stop(tool_stop1.get_val(), t);
            }
            if (t == tool_stop2) {
                thermo_present.update_v_stop(tool_stop2.get_val(), t);
            }
            applied_weight = 0;
            if (tool_weight.engaged) applied_weight += tool_weight.get_val();
            if (tool_balloon.engaged) applied_weight += tool_balloon.get_val();
        }
        else if (t == tool_weight || t == tool_balloon) {
            if (tool_weight.engaged && tool_balloon.engaged) {
                if (t == tool_weight) DetachTool(tool_balloon, popVector());
                else if (t == tool_balloon) DetachTool(tool_weight, popVector());
            }

            float v = 1f;
            if (t == tool_weight) v += dial_weight.val;
            else if (t == tool_balloon) v += dial_balloon.val;
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
            if (tool_weight.engaged) applied_weight += tool_weight.get_val();
            if (tool_balloon.engaged) applied_weight += tool_balloon.get_val();
        }

    }

    #endregion // Tools

    #region Helpers

    private void StabilizeHands() {
        //hands keep trying to run away- no idea why (this is a silly way to keep them still)
        lhand.actualhand.transform.localPosition = new Vector3(0f, 0f, 0f);
        lhand.actualhand.transform.localEulerAngles = new Vector3(0f, 0f, 90f);
        rhand.actualhand.transform.localPosition = new Vector3(0f, 0f, 0f);
        rhand.actualhand.transform.localEulerAngles = new Vector3(0f, 0f, -90f);
    }

    private void ApplyTools() {
        double delta_time = (double)Time.fixedDeltaTime;

        //apply thermo
        ambient_pressure = tool_ambientPressure.get_val();
        room_temp = tool_roomTemp.get_val();
        double weight_pressure = (applied_weight) / thermo_present.get_surfacearea_insqr(); //psi
        weight_pressure *= PSI_TO_PASCAL; //conversion from psi to pascal
        weight_pressure += ambient_pressure;
        weight_pressure = Math.Clamp(weight_pressure, ThermoMath.p_min, ThermoMath.p_max);

        // get the amount of weight to apply, based on the difference between the total weight to be applied and how much is currently applied
        double delta_weight = (weight_pressure - thermo_present.get_pressure());
        if (System.Math.Abs(delta_weight * delta_time) < World.DELTA_PRESSURE_CUTOFF) {
            // small enough step; finish transition
        }
        else {
            delta_weight *= delta_time;
        }

        double insulation_coefficient;

        if (tool_insulator.engaged) {
            insulation_coefficient = dial_percentInsulation.val;
        }
        else {
            insulation_coefficient = CONTAINER_INSULATION_COEFFICIENT;
        }

        insulation_coefficient = 1 - insulation_coefficient; // invert proportion

        // check if weight was applied
        double temperature_gradient = room_temp - thermo_present.get_temperature();

        if (System.Math.Abs(delta_weight) > 0) {
            if (insulation_coefficient == 0) {
                // perfect insulation
                thermo_present.add_pressure_insulated_per_delta_time(delta_weight, delta_time, weight_pressure, temperature_gradient); // Pressure Constrained -> Insulated ->  delta pressure
            }
            else {
                // insulation is inversely proportional to the rate of weight application
                thermo_present.add_pressure_uninsulated_per_delta_time(delta_weight, delta_time, insulation_coefficient, weight_pressure, temperature_gradient); // Pressure Constrained -> Uninsulated ->  delta pressure
            }
        }


        // heat leak
        if (toggle_heatTransfer.IsOn()) {
            double heat_transfer_delta =
                (room_temp - thermo_present.get_temperature()) // total temperature difference
                * insulation_coefficient // what percentage of that difference is shielded by insulation
                * (SPECIFIC_HEAT_CAPACITY_LIQ) // how much heat is required to raise 1 kg of water 1 Kelvin
                                               // TODO: Replace this specific heat with a function calculating based on quality parameter
                                               // for all processes not constant pressure, use c_v (vs c_p -- to be used in constant pressure)
                / delta_time * 0.5; // halve the immediacy effect so that simulation can handle the change
            // if you have some state, and know r, can calculate heat exchange (based on eqtn 2), 

            if (heat_transfer_delta != 0) {
                // insulation is inversely proportional to the rate of heat transfer (outside insulation)
                thermo_present.add_heat_per_delta_time(heat_transfer_delta, insulation_coefficient, delta_time, weight_pressure, false, temperature_gradient);
            }
        }

        // tool heat
        if (applied_heat != 0) {
            // insulation is inversely proportional to the rate of heat transfer (within insulation)
            thermo_present.add_heat_per_delta_time(applied_heat, (1 - insulation_coefficient), delta_time, weight_pressure, true, temperature_gradient);
        }
    }

    private void ProcessInputs() {
        //running blended average of hand velocity (transfers this velocity on "release object" for consistent "throwing")
        lhand.vel += (lhand.transform.position - lhand.pos) / Time.fixedDeltaTime;
        lhand.vel *= 0.5f;
        lhand.pos = lhand.transform.position;

        rhand.vel += (rhand.transform.position - rhand.pos) / Time.fixedDeltaTime;
        rhand.vel *= 0.5f;
        rhand.pos = rhand.transform.position;

        //input
        //float lhandt  = OVRInput.Get(OVRInput.Axis1D.PrimaryHandTrigger);
        //float lindext = OVRInput.Get(OVRInput.Axis1D.PrimaryIndexTrigger);
        //float rhandt  = OVRInput.Get(OVRInput.Axis1D.SecondaryHandTrigger);
        //float rindext = OVRInput.Get(OVRInput.Axis1D.SecondaryIndexTrigger);
        float lhandt = OVRInput.Get(OVRInput.RawAxis1D.LHandTrigger);
        float lindext = OVRInput.Get(OVRInput.RawAxis1D.LIndexTrigger);
        float rhandt = OVRInput.Get(OVRInput.RawAxis1D.RHandTrigger);
        float rindext = OVRInput.Get(OVRInput.RawAxis1D.RIndexTrigger);

        bool rhandraytoggle = OVRInput.GetDown(OVRInput.Button.One);
        bool lhandraytoggle = OVRInput.GetDown(OVRInput.Button.Three);

        if (rhandraytoggle) {
            rhand.ray.enabled = !rhand.ray.enabled;
        }
        if (lhandraytoggle) {
            lhand.ray.enabled = !lhand.ray.enabled;
        }

        //test effect of hands one at a time ("true" == "left hand", "false" == "right hand")
        TryHand(true, lhandt, lindext, lhand.transform.position, lhand.vel, ref lhtrigger, ref litrigger, ref lhtrigger_delta, ref litrigger_delta, ref lpos, ref lhand.obj, ref lgrabbed, ref rhand.obj, ref rgrabbed); //left hand
        TryHand(false, rhandt, rindext, rhand.transform.position, rhand.vel, ref rhtrigger, ref ritrigger, ref rhtrigger_delta, ref ritrigger_delta, ref rpos, ref rhand.obj, ref rgrabbed, ref lhand.obj, ref lgrabbed); //right hand

    }

    private void UpdateToolTexts() {

        //tooltext
        Dial d;
        for (int i = 0; i < dials.Count; i++) {
            d = dials[i];
            d.set_examined(false);
            if (d.gameObject == lgrabbed || d.gameObject == rgrabbed) d.set_examined(true);
            // TODO: update tool text
            /*
            if (t.get_val() != t.dial_dial.prev_val) {
                UpdateToolText(t);
                t.dial_dial.examined = true;
            }
            d.prev_val = t.get_val();
            */
        }

        Tool t;
        for (int i = 0; i < tools.Count; i++) {
            t = tools[i];
            if (t.text_fadable == null) {
                continue;
            }
            if (!t.text_fadable.stale) {
                if (t.text_fadable.alpha == 0f) {
                    t.textl_meshrenderer.enabled = false;
                }
                else {
                    t.textl_meshrenderer.enabled = true;
                    Color32 c = t.disabled ? new Color32(70, 70, 70, (byte)(t.text_fadable.alpha * 255)) : new Color32(0, 0, 0, (byte)(t.text_fadable.alpha * 255));
                    t.textl_tmpro.faceColor = c;
                }
            }
        }
    }

    private void ProcessErrors() {
        thermo_present.UpdateErrorState();
    }

    //safe to call if not interactable, as it will just do nothing
    bool floatNumeric(float f) {
        if (double.IsNaN(f)) return false;
        if (double.IsInfinity(f)) return false;
        return true;
    }
    void TryInteractable(GameObject actable, Vector3 hand_pos, ref Vector3 r_hand_pos) {
        //grabbing handle
        if (actable == handle_workspace) {
            float dy = (r_hand_pos.y - hand_pos.y);
            workspace.transform.position = new Vector3(workspace.transform.position.x, workspace.transform.position.y - dy, workspace.transform.position.z);
            // origin.transform.Translate(new Vector3(0, dy, 0));
        }
        else if (actable == graph) {
            Vector3 localspace = graph.transform.InverseTransformPoint(hand_pos);
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
            //grabbing dial

            Dial dd = actable.GetComponent<Dial>();

            if (dd != null) {
                dd.update_val(hand_pos, r_hand_pos);

                List<Tool> relevant_tools = dd.get_relevant_tools();
                for (int t = 0; t < relevant_tools.Count; t++) {
                    UpdateApplyTool(relevant_tools[t]);
                }
            }
        }
    }

    /*
     * This function seems to handle all possible interactions between the hand and other objects.
     * Honestly, I haven't quite got a full understanding of this ~200-line behemoth.
     */
    //"left_hand": true -> left, false -> right
    void TryHand(bool left_hand, float htrigger_val, float itrigger_val, Vector3 hand_pos, Vector3 hand_vel, ref bool ref_htrigger, ref bool ref_itrigger, ref int ref_htrigger_delta, ref int ref_itrigger_delta, ref Vector3 ref_hand_pos, ref GameObject ref_hand, ref GameObject ref_grabbed, ref GameObject ref_ohand, ref GameObject ref_ograbbed) {
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
            ref_hand_pos.x = hand_pos.x;
            ref_hand_pos.y = hand_pos.y;
        }
        else if (ref_itrigger && itrigger_val <= itrigger_threshhold) {
            ref_itrigger_delta = -1;
            ref_itrigger = false;
        }

        // fine-selection (used to extract data from UI)
        if (ref_itrigger_delta == 1 && ref_htrigger_delta != 1) {
            // TODO: this
        }

        //find new grabs
        if (ref_grabbed == null && ((ref_htrigger_delta == 1 && ref_itrigger) || (ref_htrigger && ref_itrigger_delta == 1))) {
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
                        if (t == tool_stop1) {
                            thermo_present.release_v_stop(t);
                        }
                        else if (t == tool_stop2) {
                            thermo_present.release_v_stop(t);
                        }
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

                    Cartridge c = ref_grabbed.GetComponent<Cartridge>();
                    if (c) { // newly grabbed object is a cartridge
                        GameMgr.Events.Dispatch(GameEvents.ColliderGrabbed, c.GetComponent<Collider>());
                    }
                }
            }
            //then dials
            if (ref_grabbed == null) {
                for (int i = 0; i < dials.Count; i++) {
                    if ( //dial newly grabbed
                           (left_hand && dials[i].touchable.ltouch) ||
                           (!left_hand && dials[i].touchable.rtouch)
                          ) {
                        ref_grabbed = dials[i].gameObject;
                        dials[i].touchable.grabbed = true;
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
        else if (ref_grabbed && (ref_htrigger_delta == -1 || ref_itrigger_delta == -1)) //something newly released
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

                    if (placement_thermo_reasonable) {
                        for (int i = 0; i < tools.Count; i++) {
                            Tool toDetach = tools[i];
                            if (toDetach.engaged) {
                                DetachTool(toDetach, popVector());
                            }
                        }
                        WarpPVT(placement_thermo.y, placement_thermo.x, placement_thermo.z);
                    }
                    state_dot.GetComponent<Renderer>().enabled = true;
                }

                // newly released is a cartridge
                Cartridge c = ref_grabbed.GetComponent<Cartridge>();

                if (c != null) {
                    GameMgr.Events.Dispatch(GameEvents.ColliderReleased, c.GetComponent<Collider>());
                }

                /* TODO: separate out returning functionality from tools, then add to cartridges
                if (c != null) {
                    c.GetComponent<Rigidbody>().isKinematic = false;
                    c.GetComponent<Rigidbody>().velocity = hand_vel;
                    c.GetComponent<Touchable>().grabbed = false;
                }
                */
            }

            ref_grabbed.GetComponent<Touchable>().grabbed = false;
            ref_grabbed = null;
        }

        if (ref_grabbed) TryInteractable(ref_grabbed, hand_pos, ref ref_hand_pos);

        ref_hand_pos = hand_pos;

        if (ref_grabbed == null) //still not holding anything
        {
            // Check if pressing buttons

            GameMgr.Events.Dispatch(GameEvents.CheckForPress, left_hand);
        }

        /*
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
        */
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

        bool ltouch = false;
        bool rtouch = false;
        update_touches(ref ltouch, ref rtouch);
        update_meshes(ref ltouch, ref rtouch);
    }

    private void update_touches(ref bool ltouch, ref bool rtouch) {
        for (int i = 0; i < m_pressables.Count; i++) {
            Pressable pressable = m_pressables[i];
            pressable.SetFingerTouches(ref ltouch, ref rtouch);
            continue;
        }
    }

    private void update_meshes(ref bool ltouch, ref bool rtouch) {
        if (lgrabbed) lhand.meshrenderer.materials = hand_grabbings;
        else if (ltouch) lhand.meshrenderer.materials = hand_touchings;
        else lhand.meshrenderer.materials = hand_emptys;

        if (rgrabbed) rhand.meshrenderer.materials = hand_grabbings;
        else if (rtouch) rhand.meshrenderer.materials = hand_touchings;
        else rhand.meshrenderer.materials = hand_emptys;
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

    /// <summary>
    /// Given 2 of 3 variables among p, v, t, calculate the third. Then warp.
    /// </summary>
    /// <param name="p"></param>
    /// <param name="v"></param>
    /// <param name="t"></param>
    public void WarpPVTPartial(double p, double v, double t) {
        thermo_present.warp_pv_partial(p, v, t);
    }

    public void WarpPVT(double p, double v, double t) {
        thermo_present.warp_pv(p, v, t);
    }

    /// <summary>
    ///  Used for Reach State lab questions and limits
    /// </summary>
    /// <param name="id"></param>
    /// <returns></returns>
    public double get_state_var(VarID id) {
        switch (id) {
            case VarID.VolumeStop:
                return -1; // should be handled previously
            default:
                return thermo_present.get_state_var(id);
        }
    }

    public Tuple<double, double> get_stop_vals() {
        return new Tuple<double, double>(tool_stop1.get_val(), tool_stop2.get_val());
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
        return new Vector3(UnityEngine.Random.Range(-1f, 1f), 1f, UnityEngine.Random.Range(-1f, 1f));
    }

    #endregion Helpers

    #region Limits

    private void EnforceLimits() {
        if (!ModMgr.LimitsEnabled()) {
            return;
        }

        // Check if limits have been surpassed
        bool crossedLimit = false;
        bool crossedCeiling;

        // Pressure
        if (ModMgr.CrossedLimit(VarID.Pressure, get_state_var(VarID.Pressure), out crossedCeiling)) {
            HandlePressureLimitCrossed(crossedCeiling);
            crossedLimit = true;
        }
        // Temperature
        else if (ModMgr.CrossedLimit(VarID.Temperature, get_state_var(VarID.Temperature), out crossedCeiling)) {
            HandleTemperatureLimitCrossed(crossedCeiling);
            crossedLimit = true;
        }
        // Volume
        else if (ModMgr.CrossedLimit(VarID.Volume, get_state_var(VarID.Volume), out crossedCeiling)) {
            HandleVolumeLimitCrossed(crossedCeiling);
            crossedLimit = true;
        }
        // InternalEnergy
        else if (ModMgr.CrossedLimit(VarID.InternalEnergy, get_state_var(VarID.InternalEnergy), out crossedCeiling)) {
            HandleInternalEnergyLimitCrossed(crossedCeiling);
            crossedLimit = true;
        }
        // Entropy
        else if (ModMgr.CrossedLimit(VarID.Entropy, get_state_var(VarID.Entropy), out crossedCeiling)) {
            HandleEntropyLimitCrossed(crossedCeiling);
            crossedLimit = true;
        }
        // Enthalpy
        else if (ModMgr.CrossedLimit(VarID.Enthalpy, get_state_var(VarID.Enthalpy), out crossedCeiling)) {
            HandleEnthalpyLimitCrossed(crossedCeiling);
            crossedLimit = true;
        }
        // Quality
        else if (ModMgr.CrossedLimit(VarID.Quality, get_state_var(VarID.Quality), out crossedCeiling)) {
            HandleQualityLimitCrossed(crossedCeiling);
            crossedLimit = true;
        }

        if (crossedLimit) {
            Debug.Log("[World] Crossed Limit!");
        }
    }

    private void HandlePressureLimitCrossed(bool crossedCeiling) {
        if (crossedCeiling) {
            // TODO: handle pressure limit crossed
        }
        else {
            // TODO: handle pressure limit crossed
        }
    }

    private void HandleTemperatureLimitCrossed(bool crossedCeiling) {
        if (crossedCeiling) {
            // TODO: handle temperature limit crossed
        }
        else {
            // TODO: handle temperature limit crossed
        }
    }

    private void HandleVolumeLimitCrossed(bool crossedCeiling) {
        if (crossedCeiling) {
            // TODO: handle volume limit crossed
        }
        else {
            // TODO: handle volume limit crossed
        }
    }

    private void HandleInternalEnergyLimitCrossed(bool crossedCeiling) {
        if (crossedCeiling) {
            // TODO: handle internal energy limit crossed
        }
        else {
            // TODO: handle internal energy limit crossed
        }
    }

    private void HandleEntropyLimitCrossed(bool crossedCeiling) {
        if (crossedCeiling) {
            // TODO: handle entropy limit crossed
        }
        else {
            // TODO: handle entropy limit crossed
        }
    }

    private void HandleEnthalpyLimitCrossed(bool crossedCeiling) {
        if (crossedCeiling) {
            // TODO: handle enthalpy limit crossed
        }
        else {
            // TODO: handle enthalpy limit crossed
        }
    }

    private void HandleQualityLimitCrossed(bool crossedCeiling) {
        if (crossedCeiling) {
            // TODO: handle quality limit crossed
        }
        else {
            // TODO: handle quality limit crossed
        }
    }

    #endregion // Limits

    #region Handlers

    private void HandleResetPressed(object sender, System.EventArgs args) {
        for (int i = 0; i < tools.Count; i++) {
            if (tools[i].engaged) DetachTool(tools[i], new Vector3(0.0f, 0.0f, 0.0f));
        }
        thermo_present.Reset();
    }

    private void HandleHalferPressed(object sender, System.EventArgs args) {
        SetAllHalfed(!halfed);
    }

    private void HandleRegisterPressable(Pressable pressable) {
        if (!m_pressables.Contains(pressable)) {
            m_pressables.Add(pressable);
        }
    }

    private void HandleRegisterMovable(Touchable touchable) {
        if (!movables.Contains(touchable)) {
            movables.Add(touchable);
        }
    }

    private void HandleWarpPVT() {
        // set ambient pressure to the pressure picked
        ambient_pressure = thermo_present.get_pressure();
        Debug.Log("[Warp] pressure: " + thermo_present.get_pressure() + " || min " + ThermoMath.p_min);
        dial_ambientPressure.set_val((float)((ambient_pressure - ThermoMath.p_min) / (ThermoMath.p_max - ThermoMath.p_min)));
    }


    private void HandleToolTogglePressed(Tool t) {
        for (int i = 0; i < tools.Count; i++) {
            if (t == tools[i]) {
                if (t.engaged) {
                    DetachTool(t, Vector3.zero);
                    // StoreTool(t);
                    break;
                }
                else {
                    ActivateTool(t);
                    break;
                }
            }
        }
    }

    private void HandleAllowedToolsUpdated(List<ToolType> allowed) {
        for (int i = 0; i < tools.Count; i++) {
            if (allowed.Contains(tools[i].tool_type)) {
                AllowTool(tools[i]);
            }
            else if (tools[i].always_engaged) {
                AllowTool(tools[i]);
            }
            else {
                DisallowTool(tools[i]);
            }
        }
    }

    private void HandleResetToolRestrictions() {
        for (int i = 0; i < tools.Count; i++) {
            AllowTool(tools[i]);
        }
    }

    #endregion // Handlers

}

