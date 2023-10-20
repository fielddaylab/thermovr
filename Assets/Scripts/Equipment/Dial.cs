using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR.Tools;
using TMPro;
using UnityEngine;

namespace ThermoVR.Dials
{
    public enum Effect {
        None,
        Value,      // change the value of a sim value (e.g. burner heat, or stop volume)
        Movement    // move the tool in world space
    }

    // TODO: incorporate these triggers (stops should only get enacted on release, for example)
    public enum EffectTrigger {
        OnMove,     // takes effect immediately
        OnRelease   // takes effect once dial is released
    }

    [Serializable]
    public struct EffectResponse {
        public Effect Effect;
        public float Modifier;
        public EffectTrigger Trigger;

        public EffectResponse(Effect effect, float modifier, EffectTrigger trigger) {
            Effect = effect;
            Modifier = modifier;
            Trigger = trigger;
        }
    }

    [Serializable]
    public struct EffectGroup {
        public List<Tool> ToAffect;
        public List<EffectResponse> Effects;

        public EffectGroup(List<Tool> toAffect, List<EffectResponse> effects) {
            ToAffect = toAffect;
            Effects = effects;
        }
    }

    [RequireComponent(typeof(Touchable))]
    public class Dial : MonoBehaviour
    {
        [SerializeField] private Transform max_pos;
        [SerializeField] private Transform min_pos;
        public TextMeshPro textv;
        [SerializeField] private GameObject meter; // (Knob)

        [System.NonSerialized]
        public TextMeshPro textv_tmpro;

        public int response_power; //sometimes, we want log mapping
        [System.NonSerialized]
        public float val = 0.0f; //abstract 0-1 representing knob position
        [System.NonSerialized]
        public float prev_val = 0.0f;
        public float default_val = 0;
        [System.NonSerialized]
        public Vector3 orientation_dir;
        [System.NonSerialized]
        public float map = 0.0f; //meaningful value mapped from val to [min_map,max_map]
        [System.NonSerialized]
        public float min_map = 0.0f;
        [System.NonSerialized]
        public float max_map = 1.0f;
        [System.NonSerialized]
        public string unit = "";
        // [System.NonSerialized]
        // public string display_unit = "";
        [System.NonSerialized]
        public float display_mul = 1.0f; //multiplied with map before displaying with display_unit

        public List<EffectGroup> m_effect_map; // maps all effects to each tool affected by this dial's value

        private List<Tool> relevant_tools;

        [SerializeField] private ToolActivator activator_button;
        [HideInInspector] public Touchable touchable;

        private string valFormat;

        private float total_dist;
        private Vector3 initial_offset;

        public void Init(float min_map, float max_map, string valFormat) {
            this.min_map = min_map;
            this.max_map = max_map;
            this.valFormat = valFormat;
            if (min_pos == null) {
                orientation_dir = new Vector3(1, 0, 0);
            }
            else {
                orientation_dir = (max_pos.position - min_pos.position).normalized;
            }

            relevant_tools = new List<Tool>();
            for (int g = 0; g < m_effect_map.Count; g++) {
                EffectGroup group = m_effect_map[g];

                // for each tool in that group
                for (int t = 0; t < group.ToAffect.Count; t++) {
                    Tool tool = group.ToAffect[t];
                    relevant_tools.Add(tool);
                }
            }

            if (activator_button != null) {
                activator_button.SetTools(relevant_tools);
            }

            total_dist = Vector3.Distance(max_pos.position, min_pos.position);
            initial_offset = meter.transform.localPosition;

            touchable = this.GetComponent<Touchable>();
            textv_tmpro = textv.GetComponent<TextMeshPro>();

            GameMgr.Events?.Register<Tool>(GameEvents.ActivateTool, HandleActivateTool, this)
                .Register<Tool>(GameEvents.DeactivateTool, HandleDeactivateTool, this);

            Reset();
        }

        // Update is called once per frame, and after val updated
        void Update() {
            if (!AnyToolsActive()) {
                return;
            }

            RecalibratePos();
            /*
            float magnitude = 0.05f - val * 0.1f;
            meter.transform.localPosition = Quaternion.Euler(0, 2.0f, 0) * meter.transform.forward * magnitude;
            forceMap();
            */
        }

        private void RecalibratePos() {
            Vector3 lp = meter.transform.localPosition;
            lp.x = total_dist / 2 - val * total_dist - initial_offset.x * 2;
            meter.transform.localPosition = lp;
            forceMap();
        }

        public void SetValText(float value) {
            string updateText = string.Format(this.valFormat, value);
            textv_tmpro.SetText(updateText);
        }

        private bool AnyToolsActive() {
            if (relevant_tools == null) {
                // vacuously true, I guess?
                return true;
            }

            bool anyActive = false;
            for (int i = 0; i < relevant_tools.Count; i++) {
                if (relevant_tools[i].engaged) {
                    anyActive = true;
                }
            }

            return anyActive;
        }

        public void forceMap() {
            //map = min_map+(max_map-min_map)*val;
            map = response_power > 1 ? mapSharp() : mapLinear();
        }

        public void Reset() {
            float prev_val = val;

            val = default_val;

            forceMap();

            apply_change(map, val, prev_val);

            RecalibratePos();

            // reset active button materials
            if (activator_button != null) {
                activator_button.UpdateActiveMaterial();
            }
        }

        public void set_val(float new_val) {
            val = new_val;

            forceMap();

            apply_change(map, val, prev_val);

            RecalibratePos();
        }

        public float get_val() {
            return val;
        }

        /*
         * Standard way to map from 0-1 slder "val" range to min-max "tool" range.
         */
        private float mapLinear() { return min_map + (max_map - min_map) * val; }

        /*
         * Non-linear mapping, which raises the 0-1 "val" to given power.
         * This gives us super-linear behavior, so for a 0-1 range we have "slow" growth of mapped value initially,
         * accelerating as we approach 1.
         * Useful for sliders whose tool's influence tends to apply too rapidly at small values.
         */
        private float mapSharp() { return min_map + (max_map - min_map) * Mathf.Pow(val, response_power); }

        public void update_val(Vector3 hand_pos, Vector3 r_hand_pos) {
            if (!AnyToolsActive()) {
                return;
            }

            float dx = r_hand_pos.x - hand_pos.x;
            float dy = r_hand_pos.y - hand_pos.y;
            float dz = r_hand_pos.z - hand_pos.z;

            Vector3 movement_vector = new Vector3(
                (dx) * orientation_dir.x,
                (dy) * orientation_dir.y,
                (dz) * orientation_dir.z
                );

            movement_vector *= -10f;
            // float dx = (r_hand_pos.x - hand_pos.x) * -10f;
            // constrain vector to relative orientation
            float magnitude = movement_vector.magnitude;

            bool orient_positive = (orientation_dir.x + orientation_dir.y + orientation_dir.z) >= 0; // whether orientation vector is overall positively directioned
            bool delta_positive = (dx + dy + dz) >= 0; // whether movement vector is overall positively directioned

            if (orient_positive != delta_positive) {
                // if not heading in same direction, decrease value
                magnitude *= -1;
            }

            float prev_val = val;
            // float prev_map = map;

            float new_val = Mathf.Clamp(prev_val - magnitude, 0f, 1f);
            //if this close to either end, assume user wants min/max
            if (new_val < 0.005) new_val = 0f;
            if (new_val > 0.995) new_val = 1f;

            // TODO: need a check here for if tool is enabled?
            val = new_val;
            forceMap();

            apply_change(map, new_val, prev_val);
        }

        public List<Tool> get_relevant_tools() {
            return relevant_tools;
        }

        private void apply_change(float map, float new_val, float prev_val) {
            // for each group
            for (int g = 0; g < m_effect_map.Count; g++) {
                EffectGroup group = m_effect_map[g];

                // for each tool in that group
                for (int t = 0; t < group.ToAffect.Count; t++) {
                    Tool tool = group.ToAffect[t];

                    // apply each effect listed
                    for (int e = 0; e < group.Effects.Count; e++) {
                        EffectResponse effectResponse = group.Effects[e];

                        switch (effectResponse.Effect) {
                            case Effect.Value:
                                tool.UpdateVal(map, this);
                                break;
                            case Effect.Movement:
                                tool.Move((new_val - prev_val) * effectResponse.Modifier);
                                break;
                            default:
                                break;
                        }
                    }
                }
            }
        }

        #region Handlers

        private void HandleActivateTool(Tool tool) {
            if (relevant_tools.Contains(tool)) {
                Reset();
            }
        }

        private void HandleDeactivateTool(Tool tool) {
            if (relevant_tools.Contains(tool)) {
                Reset();
            }
        }

        #endregion // Handlers
    }
}
