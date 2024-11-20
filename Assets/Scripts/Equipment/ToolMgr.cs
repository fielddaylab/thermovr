using BeauUtil;
using BeauUtil.Extensions;
using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR.Analytics;
using ThermoVR.Dials;
using ThermoVR.State;
using UnityEngine;

namespace ThermoVR.Tools
{
    public class ToolMgr : MonoBehaviour
    {
        public static ToolMgr Instance;

        public const float BURNER_MAX = 100000 * 100;
        public const float COIL_MAX = -100000 * 100;
        private const float DEFAULT_CHAMBER_PRESSURE = 101325;
        private const float DEFAULT_CHAMBER_TEMPERATURE = 300;

        [Space(5)]
        [Header("Tools")]
        [SerializeField] private Tool tool_insulator;
        [SerializeField] private Tool tool_stop1;
        [SerializeField] private Tool tool_stop2;
        [SerializeField] private Tool tool_burner;
        [SerializeField] private Tool tool_coil;
        [SerializeField] private Tool tool_weight;
        [SerializeField] private Tool tool_negativeWeight;
        [SerializeField] private Tool tool_surroundingPressure;
        [SerializeField] private Tool tool_surroundingTemp;

        List<Tool> tools;
        ParticleSystem flame; //special case

        [Space(5)]
        [Header("Dials")]
        [SerializeField] private Dial dial_stop1;
        [SerializeField] private Dial dial_stop2;
        [SerializeField] private Dial dial_burner;
        [SerializeField] private Dial dial_coil;
        [SerializeField] private Dial dial_weight;
        [SerializeField] private Dial dial_negativeWeight;
        [SerializeField] private Dial dial_surroundingPressure;
        [SerializeField] private Dial dial_surroundingTemp;
        [SerializeField] private Dial dial_percentInsulation;
        [SerializeField] private Dial dial_selector;

        [SerializeField] private PhysicalToggle toggle_heatTransfer;

        public float DetentStep = 0.05f;

        [HideInInspector] public List<Dial> Dials;
        public List<VolumeStop> VStops { get; private set; }


        [Space(5)]
        [Header("Halfables")]
        bool halfed = false;
        List<Halfable> halfables;
        [SerializeField] Pressable reset_button;
        // [SerializeField] Pressable halfer_button;

        [Space(5)]
        [Header("Other")]
        [SerializeField] private Pressable tare_heat_button;
        [SerializeField] private Pressable tare_weight_button;

        private double m_accumulatedHeatEnergy;
        private double m_accumulatedWeightEnergy;

        private AnalyticsService.SliderPanelLogData m_panelLogState;

        private void Awake() {
            if (Instance == null) {
                Instance = this;
            }
            else if (this != Instance) {
                Destroy(this);
            }
        }

        public void Init() {
            VStops = new List<VolumeStop>();

            // As we grab them, set ranges on tool dials (sliders).
            tools = new List<Tool> {
                tool_insulator,
                tool_stop1,
                tool_stop2,
                tool_burner,
                tool_coil,
                tool_weight,
                tool_negativeWeight,
                tool_surroundingPressure,
                tool_surroundingTemp,
            };

            tool_insulator.Init(Units.Quality);
            tool_stop1.Init(Units.Volume);
            tool_stop2.Init(Units.Volume);
            tool_burner.Init(Units.Heat, 0.001f);
            tool_coil.Init(Units.Heat, 0.001f);
            tool_weight.Init(Units.Weight);
            tool_negativeWeight.Init(Units.Weight);
            tool_surroundingPressure.Init(Units.AmbientPressure, 0.001f); // display in kPa
            tool_surroundingTemp.Init(Units.TemperatureK);

            // TODO: make explicit assignment
            flame = GameObject.Find("Flame").GetComponent<ParticleSystem>();

            // keep in same order as dial dashboard
            Dials = new List<Dial> {
                dial_percentInsulation,
                dial_stop2,
                dial_stop1,
                dial_weight,
                dial_negativeWeight,
                dial_burner,
                dial_coil,
                dial_surroundingTemp,
                dial_surroundingPressure,
                // dial_selector
            };

            double kg_corresponding_to_10mpa = ThermoState.surfacearea_insqr * (10 * 1453.8/*MPa->psi*/) * 0.453592/*lb->kg*/;
            // double kg_corresponding_to_2mpa = ThermoMath.surfacearea_insqr * (2 * 1453.8/*MPa->psi*/) * 0.453592/*lb->kg*/; // 10 MPa seems way too big, sooooo... we'll just do 2 MPa.

            dial_stop1.Init((float)ThermoMath.v_min, (float)ThermoMath.v_max, DigitFormat.Volume);
            dial_stop2.Init((float)ThermoMath.v_min, (float)ThermoMath.v_max, DigitFormat.Volume);
            dial_burner.Init(0f, BURNER_MAX, DigitFormat.Heat);
            dial_coil.Init(0f, COIL_MAX, DigitFormat.Heat);
            dial_weight.Init(0f, (float)kg_corresponding_to_10mpa / 5.0f, DigitFormat.Weight);
            dial_negativeWeight.Init(0f, -(float)kg_corresponding_to_10mpa / 5.0f, DigitFormat.Weight); // 500.0f
            dial_surroundingPressure.Init((float)ThermoMath.p_min, (float)ThermoMath.p_max, DigitFormat.AmbientPressure);
            dial_surroundingTemp.Init(273, (float)ThermoMath.t_max, DigitFormat.TemperatureK);
            dial_percentInsulation.Init(0f, 100, DigitFormat.Percent, preserveDuringLockedOverride: true);
            // dial_selector.Init(0, 1, DigitFormat.Percent);

            toggle_heatTransfer.IsActiveImpl = () => { return tool_surroundingTemp.allowed; };

            ResetDefaults();

            // Initialize Buttons
            reset_button.OnPress += HandleResetPressed;
            tare_heat_button.OnPress += HandleTareHeatPressed;
            tare_weight_button.OnPress += HandleTareWeightPressed;


            EventMgr.Events?.Register<STuple<double, double, double>>(GameEvents.WarpPVT, HandleWarpPVT);

            EventMgr.Events?.Register<Tool>(GameEvents.PressedToolToggle, HandleToolTogglePressed);

            EventMgr.Events?.Register<List<ToolType>>(GameEvents.UpdateAllowedTools, HandleAllowedToolsUpdated);
            EventMgr.Events?.Register(GameEvents.ResetToolRestrictions, HandleResetToolRestrictions);

            EventMgr.Events?.Register<GameObject>(GameEvents.ObjectGrabbed, HandleObjectGrabbed);
            EventMgr.Events?.Register<GameObject>(GameEvents.ObjectReleased, HandleObjectReleased);

            EventMgr.Events?.Register<int>(GameEvents.NudgeUpClicked, HandleNudgeUpClicked);
            EventMgr.Events?.Register<int>(GameEvents.NudgeDownClicked, HandleNudgeDownClicked);

            m_panelLogState = new AnalyticsService.SliderPanelLogData();
            m_panelLogState.Insulation = new AnalyticsService.SliderSettings();
            m_panelLogState.LowerStop = new AnalyticsService.SliderSettings();
            m_panelLogState.UpperStop = new AnalyticsService.SliderSettings();
            m_panelLogState.Weight = new AnalyticsService.SliderSettings();
            m_panelLogState.NegativeWeight = new AnalyticsService.SliderSettings();
            m_panelLogState.Heat = new AnalyticsService.SliderSettings();
            m_panelLogState.Cooling = new AnalyticsService.SliderSettings();
            m_panelLogState.ChamberPressure = new AnalyticsService.SliderSettings();
            m_panelLogState.ChamberTemperature = new AnalyticsService.SliderSettings();
        }

        #region Accessors

        public double GetToolVal(ToolType type, int uniqueID = 0) {
            switch (type) {
                case ToolType.Burner:
                    return tool_burner.GetVal();
                case ToolType.Coil:
                    return tool_coil.GetVal();
                case ToolType.Weight:
                    return tool_weight.GetVal();
                case ToolType.NegativeWeight:
                    return tool_negativeWeight.GetVal();
                case ToolType.Insulator:
                    return tool_insulator.GetVal();
                case ToolType.SurroundingTemperature:
                    return tool_surroundingTemp.GetVal();
                case ToolType.SurroundingPressure:
                    return tool_surroundingPressure.GetVal();
                case ToolType.Stops:
                    if (uniqueID == 1) { return tool_stop1.GetVal(); }
                    else if (uniqueID == 2) { return tool_stop2.GetVal(); }
                    else { break; }
                default:
                    break;
            }

            return -1;
        }

        public double GetDialVal(ToolType type, int uniqueID = 0) {
            switch (type) {
                case ToolType.Burner:
                    return dial_burner.get_val();
                case ToolType.Coil:
                    return dial_coil.get_val();
                case ToolType.Weight:
                    return dial_weight.get_val();
                case ToolType.NegativeWeight:
                    return dial_negativeWeight.get_val();
                case ToolType.Insulator:
                    return dial_percentInsulation.get_val();
                case ToolType.SurroundingTemperature:
                    return dial_surroundingTemp.get_val();
                case ToolType.SurroundingPressure:
                    return dial_surroundingPressure.get_val();
                case ToolType.Stops:
                    if (uniqueID == 1) { return dial_stop1.get_val(); }
                    else if (uniqueID == 2) { return dial_stop2.get_val(); }
                    else { break; }
                default:
                    break;
            }

            return -1;
        }

        public bool IsToolEngaged(ToolType type, int uniqueID = 0) {
            switch (type) {
                case ToolType.Burner:
                    return tool_burner.engaged;
                case ToolType.Coil:
                    return tool_coil.engaged;
                case ToolType.Weight:
                    return tool_weight.engaged;
                case ToolType.NegativeWeight:
                    return tool_negativeWeight.engaged;
                case ToolType.Insulator:
                    return tool_insulator.engaged;
                case ToolType.SurroundingTemperature:
                    return tool_surroundingTemp.engaged;
                case ToolType.SurroundingPressure:
                    return tool_surroundingPressure.engaged;
                case ToolType.Stops:
                    if (uniqueID == 1) { return tool_stop1.engaged; }
                    else if (uniqueID == 2) { return tool_stop2.engaged; }
                    else { break; }
                default:
                    break;
            }

            return false;
        }

        public bool IsHeatToggleOn() {
            return toggle_heatTransfer.IsOn();
        }

        public int IdentifyStop(Tool tool)
        {
            if (tool == tool_stop1)
            {
                return 1;
            }
            else if (tool == tool_stop2)
            {
                return 2;
            }

            return 0;
        }

        #endregion // Accessors

        #region Volume Stops


        /// <summary>
        /// Remove clamp volume bounds
        /// </summary>
        public void ReleaseVStop(Tool source) {
            if (!StopExists(source)) {
                return;
            }

            RemoveStop(source);
        }

        private void AddVStop(double v_stop, Tool source, ConstrainType constraint) {
            if (StopExists(source)) {
                return;
            }

            // constrain stop's values to global bounds
            v_stop = MathUtility.Clampd(v_stop, ThermoMath.v_min, ThermoMath.v_max);

            VolumeStop new_stop = new VolumeStop(v_stop, source, constraint);
            VStops.Add(new_stop);
        }

        private bool StopExists(Tool source) {
            if (VStops == null) { return false; }

            for (int i = 0; i < VStops.Count; i++) {
                if (VStops[i].Source == source) {
                    return true;
                }
            }
            return false;
        }

        private void RemoveStop(Tool source) {

            for (int i = 0; i < VStops.Count; i++) {
                if (VStops[i].Source == source) {
                    VStops.RemoveAt(i);
                    return;
                }
            }
        }

        public void UpdateVStop(double v_stop_val, Tool source) {
            for (int i = 0; i < VStops.Count; i++) {
                if (VStops[i].Source == source) {
                    VolumeStop temp_stop = VStops[i];
                    temp_stop.Volume = v_stop_val;
                    VStops[i] = temp_stop;
                    return;
                }
            }
        }

        public bool CheckForStopsViolation(double currV)
        {
            for (int i = 0; i < VStops.Count; i++)
            {
                if (VStops[i].StoppedDir == ConstrainType.Min && currV < VStops[i].Volume)
                {
                    return true;
                }
                if (VStops[i].StoppedDir == ConstrainType.Max && currV > VStops[i].Volume)
                {
                    return true;
                }
            }

            return false;
        }

        #endregion // Volume Stops

        public void ActivateTool(Tool t, bool auto = true) {
            // trigger tool's entry animations
            t.TriggerActivation();

            GameObject o = t.gameObject;
            t.was_engaged = t.engaged;
            t.engaged = true;
            int uniqueStopID = 0;
            if (t == tool_stop1) {
                AddVStop(tool_stop1.GetVal(), t, ConstrainType.Max);
                uniqueStopID = 1;
            }
            else if (t == tool_stop2) {
                AddVStop(tool_stop2.GetVal(), t, ConstrainType.Min);
                uniqueStopID = 2;
            }

            if (!auto)
            {
                EventMgr.Events?.Dispatch(GameEvents.ToolTogglePressed, EvtArgs.Create(new STuple<ToolType, bool, bool, int>(t.tool_type, true, false, uniqueStopID)));
            }

            EventMgr.Events?.Dispatch(GameEvents.ActivateTool, EvtArgs.Ref(t));
            UpdateApplyTool(t);

            Halfable h = o.GetComponent<Halfable>();
            if (h != null) h.setHalf(halfed); //conform to half-ness while engaged

            RecordPanelUpdate();
        }

        public void DeactivateTool(Tool t, bool auto = true) {
            bool toolReset = false;
            if (!t.always_engaged) {
                if (!(auto && t.preserve_during_locked))
                {
                    // Trigger tool's deactivation animations
                    t.TriggerDeactivation();

                    t.engaged = false;
                    toolReset = true;
                }
            }
            int uniqueStopID = 0;
            if (t == tool_stop1) {
                ReleaseVStop(t);
                uniqueStopID = 1;
            }
            else if (t == tool_stop2) {
                ReleaseVStop(t);
                uniqueStopID = 2;
            }

            if (!auto)
            {
                EventMgr.Events?.Dispatch(GameEvents.ToolTogglePressed, EvtArgs.Create(new STuple<ToolType, bool, bool, int>(t.tool_type, false, toolReset, uniqueStopID)));
            }

            EventMgr.Events?.Dispatch(GameEvents.DeactivateTool, EvtArgs.Ref(t));
            UpdateApplyTool(t);

            RecordPanelUpdate();
        }

        public void EngageTool(Tool t) {
            t.TriggerEngage();
        }

        public void DisengageTool(Tool t) {
            t.TriggerDisengage();
        }

        public void BeginAdjustTool(Tool t) {
            t.TriggerBeginAdjust();
        }

        public void EndAdjustTool(Tool t) {
            t.TriggerEndAdjust();

            if (t == tool_stop1)
            {
                ReleaseVStop(t);
            }
            else if (t == tool_stop2)
            {
                ReleaseVStop(t);
            }

            if (t == tool_stop1)
            {
                AddVStop(tool_stop1.GetVal(), t, ConstrainType.Max);
            }
            else if (t == tool_stop2)
            {
                AddVStop(tool_stop2.GetVal(), t, ConstrainType.Min);
            }
        }


        public void AllowTool(Tool t) {
            // TODO: show on buttons (dispatch event)
            t.allowed = true;

            EventMgr.Events?.Dispatch(GameEvents.AllowTool, EvtArgs.Ref(t));

            RecordPanelUpdate();
        }

        public void DisallowTool(Tool t) {
            // TODO: show on buttons (dispatch event)
            t.allowed = false;

            DeactivateTool(t);

            EventMgr.Events?.Dispatch(GameEvents.DisallowTool, EvtArgs.Ref(t));

            RecordPanelUpdate();
        }

        public void UpdateApplyTool(Tool t) //alters "applied_x"
        {
            if (t == tool_insulator) {
                //do nothing
                return;
            }
            else if (t == tool_burner || t == tool_coil) {
                if (tool_burner.engaged && tool_coil.engaged) {
                    if (t == tool_burner) DeactivateTool(tool_coil);
                    if (t == tool_coil) DeactivateTool(tool_burner);
                }

                if (t == tool_burner) {
                    if (tool_burner.GetVal() == 0.0f) {
                        var e = flame.emission;
                        e.enabled = false;
                    }
                    else {
                        var e = flame.emission;
                        e.enabled = true;
                    }
                    var vel = flame.velocityOverLifetime;
                    vel.speedModifierMultiplier = Mathf.Lerp(0.1f, 0.5f, t.GetVal());
                }
                else if (t == tool_coil) {
                    //TODO: coil visuals?
                }
            }
            else if (t == tool_stop1 || t == tool_stop2) {
                if (t == tool_stop1) {
                    UpdateVStop(tool_stop1.GetVal(), t);
                }
                if (t == tool_stop2) {
                    UpdateVStop(tool_stop2.GetVal(), t);
                }
            }
            else if (t == tool_weight || t == tool_negativeWeight) {
                if (tool_weight.engaged && tool_negativeWeight.engaged) {
                    if (t == tool_weight) DeactivateTool(tool_negativeWeight);
                    else if (t == tool_negativeWeight) DeactivateTool(tool_weight);
                }

                float v = 1f;
                if (t == tool_weight) v += dial_weight.val;
                else if (t == tool_negativeWeight) v += dial_negativeWeight.val;
            }

            RecordPanelUpdate();
        }

        private void SetAllHalfed(bool h) {
            halfed = h;
            for (int i = 0; i < halfables.Count; i++) {
                halfables[i].setHalf(halfed);
            }
            //special case, only halfed when engaged
            if (!tool_coil.engaged) tool_coil.gameObject.GetComponent<Halfable>().setHalf(false);
            if (!tool_insulator.engaged) tool_insulator.gameObject.GetComponent<Halfable>().setHalf(false);
        }

        /// <summary>
        /// Deactivates all tools that aren't omnipresent
        /// </summary>
        /// <param name="forceDetach">Force tools to detach on initial pass</param>
        public void DeactivateAllTools(bool forceDetach = false) {
            for (int i = 0; i < tools.Count; i++) {
                Tool toDetach = tools[i];
                if (toDetach.engaged || forceDetach) {
                    DeactivateTool(toDetach);
                }
            }
        }

        public void ActivateConstantTools() {
            for (int i = 0; i < tools.Count; i++) {
                Tool toActivate = tools[i];
                if (toActivate.always_engaged) {
                    ActivateTool(toActivate);
                }
            }
        }

        public double GetAppliedHeat() {
            double applied_heat = 0;
            if (tool_burner.engaged) applied_heat += tool_burner.GetVal();
            if (tool_coil.engaged) applied_heat += tool_coil.GetVal();

            return applied_heat;
        }

        public void RecordToAccumulatedHeatEnergy(double deltaHeat)
        {
            m_accumulatedHeatEnergy += deltaHeat;

            EventMgr.Events.Dispatch(GameEvents.AccumHeatEnergyUpdated);
        }

        public void ResetAccumulatedHeatEnergy()
        {
            m_accumulatedHeatEnergy = 0;

            EventMgr.Events.Dispatch(GameEvents.AccumHeatEnergyUpdated);

            EventMgr.Events.Dispatch(GameEvents.ClickClearHeatMeter);
        }

        public void RecordToAccumulatedWeightEnergy(double deltaWeight)
        {
            m_accumulatedWeightEnergy += Math.Abs(deltaWeight);

            EventMgr.Events.Dispatch(GameEvents.AccumWeightEnergyUpdated);

            EventMgr.Events.Dispatch(GameEvents.ClickClearWorkMeter);
        }

        public void ResetAccumulatedWeightEnergy()
        {
            m_accumulatedWeightEnergy = 0;

            EventMgr.Events.Dispatch(GameEvents.AccumWeightEnergyUpdated);
        }

        public double GetAccumulatedHeatEnergy()
        {
            return m_accumulatedHeatEnergy;
        }

        public double GetAccumulatedWeightEnergy()
        {
            return m_accumulatedWeightEnergy;
        }

        public double GetAppliedWeight() {
            double applied_weight = 0;
            if (tool_weight.engaged) {
                applied_weight += tool_weight.GetVal();
            }
            if (tool_negativeWeight.engaged) {
                applied_weight += tool_negativeWeight.GetVal();
            }

            return applied_weight;
        }

        private void RecordPanelUpdate()
        {
            m_panelLogState.Insulation.Enabled = tool_insulator.enabled;
            m_panelLogState.Insulation.SliderVal = tool_insulator.GetVal();
            m_panelLogState.LowerStop.Enabled = tool_stop2.enabled;
            m_panelLogState.LowerStop.SliderVal = tool_stop2.GetVal();
            m_panelLogState.UpperStop.Enabled = tool_stop1.enabled;
            m_panelLogState.UpperStop.SliderVal = tool_stop1.GetVal();
            m_panelLogState.Weight.Enabled = tool_weight.enabled;
            m_panelLogState.Weight.SliderVal = tool_weight.GetVal();
            m_panelLogState.NegativeWeight.Enabled = tool_negativeWeight.enabled;
            m_panelLogState.NegativeWeight.SliderVal = tool_negativeWeight.GetVal();
            m_panelLogState.Heat.Enabled = tool_burner.enabled;
            m_panelLogState.Heat.SliderVal = tool_burner.GetVal();
            m_panelLogState.Cooling.Enabled = tool_coil.enabled;
            m_panelLogState.Cooling.SliderVal = tool_coil.GetVal();
            m_panelLogState.ChamberPressure.Enabled = tool_surroundingPressure.enabled;
            m_panelLogState.ChamberPressure.SliderVal = tool_surroundingPressure.GetVal();
            m_panelLogState.ChamberTemperature.Enabled = tool_surroundingTemp.enabled;
            m_panelLogState.ChamberTemperature.SliderVal = tool_surroundingTemp.GetVal();

            EventMgr.Events.Dispatch(GameEvents.SliderPanelUpdated, EvtArgs.Box(m_panelLogState));
        }

        #region Handlers


        private void HandleResetPressed(object sender, System.EventArgs args) {
            if (GameMgr.I.AudioEnabled) { reset_button.ClickAudio.Play(); }

            EventMgr.Events.Dispatch(GameEvents.ResetPressed);

            ResetDefaults();
        }

        private void HandleTareHeatPressed(object sender, System.EventArgs args)
        {
            Debug.Log("[ToolMgr] Heat tare pressed!");

            ResetAccumulatedHeatEnergy();
        }

        private void HandleTareWeightPressed(object sender, System.EventArgs args)
        {
            Debug.Log("[ToolMgr] Weight tare pressed!");

            ResetAccumulatedWeightEnergy();
        }

        private void ResetDefaults() {
            DeactivateAllTools(true);
            ActivateConstantTools();

            float targetMap = (float)((DEFAULT_CHAMBER_PRESSURE - ThermoMath.p_min) / (ThermoMath.p_max - ThermoMath.p_min));
            dial_surroundingPressure.set_mapped_val(targetMap);
            SetChamberTemp(DEFAULT_CHAMBER_TEMPERATURE);

            // Insulator starts engaged
            ActivateTool(tool_insulator);

            toggle_heatTransfer.ResetToggle();
        }

        public void SetChamberTemp(float targetVal)
        {
            float targetMap = (float)((targetVal - 273) / (ThermoMath.t_max - 273));
            dial_surroundingTemp.set_mapped_val(targetMap);
        }

        public void SetInsulation(float targetVal)
        {
            if (targetVal > 1)
            {
                // convert to decimal if given as percent 
                targetVal /= 100;
            }

            dial_percentInsulation.set_mapped_val(targetVal);
        }

        private void HandleHalferPressed(object sender, System.EventArgs args) {
            SetAllHalfed(!halfed);
        }

        private void HandleWarpPVT(STuple<double, double, double> pvt) {
            // Pop off applied tools after warp
            DeactivateAllTools(true);

            // set ambient pressure to the pressure picked
            float targetMap = (float)((pvt.Item1 - ThermoMath.p_min) / (ThermoMath.p_max - ThermoMath.p_min));
            dial_surroundingPressure.set_mapped_val(targetMap);

            // Insulator starts engaged
            ActivateTool(tool_insulator);

            // TODO: activate toggle on a task by task basis
            // toggle_heatTransfer.ResetToggle();
        }

        private void HandleToolTogglePressed(Tool t) {
            for (int i = 0; i < tools.Count; i++) {
                if (t == tools[i]) {
                    if (t.engaged) {
                        DeactivateTool(t, false);
                        break;
                    }
                    else {
                        ActivateTool(t, false);
                        break;
                    }
                }
            }
        }


        private void HandleAllowedToolsUpdated(List<ToolType> allowed) {
            // Between lab tasks

            for (int i = 0; i < tools.Count; i++) {
                if (allowed.Contains(tools[i].tool_type)) {
                    AllowTool(tools[i]);
                }
                /*
                else if (tools[i].always_engaged) {
                    AllowTool(tools[i]);
                }
                */
                else {
                    DisallowTool(tools[i]);
                }
            }
        }

        private void HandleResetToolRestrictions() {
            // When lab is removed
            ResetDefaults();

            for (int i = 0; i < tools.Count; i++) {
                AllowTool(tools[i]);
            }
        }


        private void HandleObjectGrabbed(GameObject obj) {
            // Check if it was a tool dial knob
            List<Tool> relevantTools = null;

            // Check if it was a tool dial knob
            foreach (var dial in Dials) {
                if (dial.gameObject == obj) {
                    relevantTools = dial.get_relevant_tools();
                    dial.OnGrabbed();
                    break;
                }
            }

            if (relevantTools != null) {
                foreach (var tool in relevantTools) {
                    BeginAdjustTool(tool);
                }
            }
        }

        private void HandleObjectReleased(GameObject obj) {
            List<Tool> relevantTools = null;

            // Check if it was a tool dial knob
            foreach(var dial in Dials) {
                if (dial.gameObject == obj) {
                    relevantTools = dial.get_relevant_tools();
                    dial.DeactivateNudge();
                    dial.OnReleased();
                    break;
                }
            }

            if (relevantTools != null) {
                foreach(var tool in relevantTools) {
                    EndAdjustTool(tool);
                }
            }
        }

        private void HandleNudgeUpClicked(int dialIndex)
        {
            Dials[dialIndex].NudgeUp(false, Vector3.zero);
        }

        private void HandleNudgeDownClicked(int dialIndex)
        {
            Dials[dialIndex].NudgeDown(false, Vector3.zero);
        }

        #endregion // Handlers
    }
}