/*
DOCUMENTATION- phil, 12/16/19

This class manages all the interaction in the scene.
It relies on ThermoState to keep track of any thermodynamic-centric state, but other than that, this is responsible for everything moving about the scene.
It should be instantiated as a game object "Oracle" at the root of the scene heirarchy.
There are unfortunately somewhat inconsistent patterns of what variables are defined publicly via the editor inspector, and which are set in code, though I tried to err toward the latter where possible.
*/

using System.Collections.Generic;
using UnityEngine;
using ThermoVR.Controls;
using ThermoVR.Dials;
using ThermoVR.Tools;
using ThermoVR;
using ThermoVR.UI.GraphElements;
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

    #endregion // Consts

    #region Inspector

    public WorldModMgr ModMgr;
    [SerializeField] private ToolMgr ToolMgr;

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
    [Header("Dot Placement")]
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

        GameMgr.Events?.Register(GameEvents.ResetPressed, HandleResetPressed);


        movables = new List<Touchable>();
    }


    public void Init() {
        // All this code does, at end of day, is find all the objects to manage,
        // and set initial values and such as needed.
        // thermo_present = GameObject.Find("Oracle").GetComponent<ThermoPresent>();

        hand_emptys = new Material[] { hand_empty };
        hand_touchings = new Material[] { hand_touching };
        hand_grabbings = new Material[] { hand_grabbing };

        if (lhand != null) {
            lhand.Init(hand_emptys);
        }
        if (rhand != null) {
            rhand.Init(hand_emptys);
        }

        ToolMgr.Init();

        room_temp = ToolMgr.GetToolVal(ToolType.SurroundingTemperature); // in K

        // Gather pressables
        m_pressables = new List<Pressable>();
        GameMgr.Events.Dispatch(GameEvents.GatherPressables);

        // toggle_heatTransfer.Init();
        // toggle_heatTransfer.Pressable.PressCompleted += HandleHeatTransferToggle;

        // Initialize Tablet (and corresponding buttons)
        tablet.Init();

        workspace = GameObject.Find("Workspace");
        handle_workspace_touchable = handle_workspace.GetComponent<Touchable>();

        graph = GameObject.Find("Graph");
        graph_touchable = graph.GetComponent<Touchable>();
        state_dot = GameObject.Find("gstate");
        placement_dot = GameObject.Find("tstate");
        placement_dot.GetComponent<Renderer>().enabled = false;
        placement_thermo_reasonable = false;

        movables.Add(tablet.touchable);
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

        if (lhand != null && rhand != null) {
            StabilizeHands();
        }

        ApplyTools();

        EnforceLimits();

        if (lhand != null && rhand != null) {
            ProcessInputs();

            UpdateGrabVis();
        }


        ProcessErrors();
    }

    #endregion // Callbacks

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
        ambient_pressure = ToolMgr.GetToolVal(ToolType.SurroundingPressure);
        room_temp = ToolMgr.GetToolVal(ToolType.SurroundingTemperature);
        double weight_pressure = (ToolMgr.GetAppliedWeight()) / ThermoMath.surfacearea_insqr; //psi
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

        if (ToolMgr.IsToolEngaged(ToolType.Insulator)) {
            insulation_coefficient = ToolMgr.GetDialVal(ToolType.Insulator);
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
        if (ToolMgr.IsHeatToggleOn()) {
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

        double applied_heat = ToolMgr.GetAppliedHeat();

        // tool heat
        if (applied_heat != 0) {
            // insulation is inversely proportional to the rate of heat transfer (within insulation)
            thermo_present.add_heat_per_delta_time(applied_heat, (1 - insulation_coefficient), delta_time, weight_pressure, true, temperature_gradient);
        }

        // Debug.Log("[warp] current temp: " + thermo_present.get_temperature()); // useful for determining exact temp needed for set values in labs
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

    private void ProcessErrors() {
        thermo_present.UpdateErrorState();
    }

    /// <summary>
    /// Tries to grab things
    /// </summary>
    /// <param name="actable"></param>
    /// <param name="hand_pos">curr hand position</param>
    /// <param name="r_hand_pos">ref to prev hand position</param>
    public void TryInteractable(GameObject actable, Vector3 hand_pos, ref Vector3 r_hand_pos) {
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

            if (MathUtility.floatNumeric(localguess.x) && MathUtility.floatNumeric(localguess.y) && MathUtility.floatNumeric(localguess.z)) {
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
                    ToolMgr.UpdateApplyTool(relevant_tools[t]);
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
                for (int i = 0; i < ToolMgr.Dials.Count; i++) {
                    if ( //dial newly grabbed
                           (left_hand && ToolMgr.Dials[i].touchable.ltouch) ||
                           (!left_hand && ToolMgr.Dials[i].touchable.rtouch)
                          ) {
                        ref_grabbed = ToolMgr.Dials[i].gameObject;
                        ToolMgr.Dials[i].touchable.grabbed = true;
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
                    ToolMgr.DeactivateAllTools();
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
        return new Tuple<double, double>(ToolMgr.GetToolVal(ToolType.Stops, 1), ToolMgr.GetToolVal(ToolType.Stops, 2));
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

    private void HandleResetPressed() {
        thermo_present.Reset();
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

    #endregion // Handlers

}

