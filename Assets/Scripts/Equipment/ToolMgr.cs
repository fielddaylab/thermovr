using BeauUtil;
using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR.Dials;
using UnityEngine;

namespace ThermoVR.Tools
{
    public class ToolMgr : MonoBehaviour
    {
        public static ToolMgr Instance;

        public const float BURNER_MAX = 100000;
        public const float COIL_MAX = -100000;

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

        [SerializeField] private PhysicalToggle toggle_heatTransfer;

        [HideInInspector] public List<Dial> Dials;
        public List<VolumeStop> VStops { get; private set; }


        [Space(5)]
        [Header("Halfables")]
        bool halfed = false;
        List<Halfable> halfables;
        [SerializeField] Pressable reset_button;
        [SerializeField] Pressable halfer_button;

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

            DeactivateAllTools(true);

            // Insulator starts engaged
            ActivateTool(tool_insulator);

            Dials = new List<Dial> {
                dial_stop1,
                dial_stop2,
                dial_burner,
                dial_coil,
                dial_weight,
                dial_negativeWeight,
                dial_surroundingPressure,
                dial_surroundingTemp,
                dial_percentInsulation
            };

            double kg_corresponding_to_10mpa = ThermoMath.surfacearea_insqr * (10 * 1453.8/*MPa->psi*/) * 0.453592/*lb->kg*/;
            // double kg_corresponding_to_2mpa = ThermoMath.surfacearea_insqr * (2 * 1453.8/*MPa->psi*/) * 0.453592/*lb->kg*/; // 10 MPa seems way too big, sooooo... we'll just do 2 MPa.

            dial_stop1.Init((float)ThermoMath.v_min, (float)ThermoMath.v_max, DigitFormat.Volume);
            dial_stop2.Init((float)ThermoMath.v_min, (float)ThermoMath.v_max, DigitFormat.Volume);
            dial_burner.Init(0f, BURNER_MAX, DigitFormat.Heat);
            dial_coil.Init(0f, COIL_MAX, DigitFormat.Heat);
            dial_weight.Init(0f, (float)kg_corresponding_to_10mpa / 5.0f, DigitFormat.Weight);
            dial_negativeWeight.Init(0f, -(float)kg_corresponding_to_10mpa / 5.0f, DigitFormat.Weight); // 500.0f
            dial_surroundingPressure.Init((float)ThermoMath.p_min, (float)ThermoMath.p_max, DigitFormat.AmbientPressure);
            dial_surroundingTemp.Init(273, 366, DigitFormat.TemperatureK); // -100 to 200 fahrenheit // default val of 0.55 sets to 292 kelvin (72 degrees fahrenheit)
            dial_percentInsulation.Init(0f, 100, DigitFormat.Percent);

            // Initialize Buttons
            reset_button.OnPress += HandleResetPressed;
            halfer_button.OnPress += HandleHalferPressed;

            /*
            halfables = new List<Halfable> {
                // GameObject.Find("Container").GetComponent<Halfable>(),
                GameObject.Find("Tool_Insulator Variant").GetComponent<Halfable>(),
                GameObject.Find("Tool_Coil Variant").GetComponent<Halfable>()
            };

            SetAllHalfed(true);
            */

            GameMgr.Events?.Register<Tuple<double, double, double>>(GameEvents.WarpPVT, HandleWarpPVT);

            GameMgr.Events?.Register<Tool>(GameEvents.PressedToolToggle, HandleToolTogglePressed);

            GameMgr.Events?.Register<List<ToolType>>(GameEvents.UpdateAllowedTools, HandleAllowedToolsUpdated);
            GameMgr.Events?.Register(GameEvents.ResetToolRestrictions, HandleResetToolRestrictions);
        }

        #region Accessors

        public double GetToolVal(ToolType type, int uniqueID = 0) {
            switch (type) {
                case ToolType.Burner:
                    return tool_burner.get_val();
                case ToolType.Coil:
                    return tool_coil.get_val();
                case ToolType.Weight:
                    return tool_weight.get_val();
                case ToolType.NegativeWeight:
                    return tool_negativeWeight.get_val();
                case ToolType.Insulator:
                    return tool_insulator.get_val();
                case ToolType.SurroundingTemperature:
                    return tool_surroundingTemp.get_val();
                case ToolType.SurroundingPressure:
                    return tool_surroundingPressure.get_val();
                case ToolType.Stops:
                    if (uniqueID == 1) { return tool_stop1.get_val(); }
                    else if (uniqueID == 2) { return tool_stop2.get_val(); }
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

        private void AddVStop(double v_stop, Tool source) {
            if (StopExists(source)) {
                return;
            }

            // constrain stop's values to global bounds
            v_stop = MathUtility.Clampd(v_stop, ThermoMath.v_min, ThermoMath.v_max);

            VolumeStop new_stop = new VolumeStop(v_stop, source);
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

        #endregion // Volume Stops

        public void ActivateTool(Tool t) {
            // TODO: trigger tool's entry animations
            // tool.ActivateRoutine
            t.gameObject.SetActive(true);

            GameObject o = t.gameObject;
            t.engaged = true;
            if (t == tool_stop1) {
                AddVStop(tool_stop1.get_val(), t);
            }
            else if (t == tool_stop2) {
                AddVStop(tool_stop2.get_val(), t);
            }
            GameMgr.Events?.Dispatch(GameEvents.ActivateTool, t);
            UpdateApplyTool(t);

            Halfable h = o.GetComponent<Halfable>();
            if (h != null) h.setHalf(halfed); //conform to half-ness while engaged
        }

        public void DeactivateTool(Tool t) {
            if (!t.always_engaged) {
                t.engaged = false;
                t.gameObject.SetActive(false);
                // TODO: trigger tool's entry animations
                // tool.DeactivateRoutine
            }
            if (t == tool_stop1) {
                ReleaseVStop(t);
            }
            else if (t == tool_stop2) {
                ReleaseVStop(t);
            }
            GameMgr.Events?.Dispatch(GameEvents.DeactivateTool, t);
            UpdateApplyTool(t);
        }

        public void AllowTool(Tool t) {
            // TODO: show buttons
            t.gameObject.SetActive(true);
        }

        public void DisallowTool(Tool t) {
            // TODO: disable buttons
            DeactivateTool(t);
            t.gameObject.SetActive(false);
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
            }
            else if (t == tool_stop1 || t == tool_stop2) {
                if (t == tool_stop1) {
                    UpdateVStop(tool_stop1.get_val(), t);
                }
                if (t == tool_stop2) {
                    UpdateVStop(tool_stop2.get_val(), t);
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

        public double GetAppliedHeat() {
            double applied_heat = 0;
            if (tool_burner.engaged) applied_heat += tool_burner.get_val();
            if (tool_coil.engaged) applied_heat += tool_coil.get_val();

            return applied_heat;
        }

        public double GetAppliedWeight() {
            double applied_weight = 0;
            if (tool_weight.engaged) applied_weight += tool_weight.get_val();
            if (tool_negativeWeight.engaged) applied_weight += tool_negativeWeight.get_val();

            return applied_weight;
        }

        #region Handlers


        private void HandleResetPressed(object sender, System.EventArgs args) {
            GameMgr.Events.Dispatch(GameEvents.ResetPressed);

            for (int i = 0; i < tools.Count; i++) {
                if (tools[i].engaged) DeactivateTool(tools[i]);
            }

            VStops.Clear();
        }

        private void HandleHalferPressed(object sender, System.EventArgs args) {
            SetAllHalfed(!halfed);
        }

        private void HandleWarpPVT(Tuple<double, double, double> pvt) {
            // set ambient pressure to the pressure picked
            dial_surroundingPressure.set_val((float)((pvt.Item1 - ThermoMath.p_min) / (ThermoMath.p_max - ThermoMath.p_min)));
        }

        private void HandleToolTogglePressed(Tool t) {
            for (int i = 0; i < tools.Count; i++) {
                if (t == tools[i]) {
                    if (t.engaged) {
                        DeactivateTool(t);
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
}