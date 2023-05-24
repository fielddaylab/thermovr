using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using TMPro;

/**
 * Runs a "dial" object, which is really a slider.
 * This controls the amount of heat I/O or weight from a tool
 **/
public class Dial : MonoBehaviour
{
    [SerializeField] private Transform max_pos;
    [SerializeField] private Transform min_pos;

    public int response_power; //sometimes, we want log mapping
    [System.NonSerialized]
    public float val = 0.0f; //abstract 0-1 representing knob position
    [System.NonSerialized]
    public float prev_val = 0.0f;
    [System.NonSerialized]
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
    [System.NonSerialized]
    public string display_unit = "";
    [System.NonSerialized]
    public float display_mul = 1.0f; //multiplied with map before displaying with display_unit
    public GameObject tool;

    [System.NonSerialized]
    public bool examined = false;
    GameObject meter;

    public void Init(float min_map, float max_map, string unit, string display_unit, float display_mul) {
        meter = gameObject.transform.GetChild(0).gameObject;

        this.min_map = min_map;
        this.max_map = max_map;
        this.unit = unit;
        this.display_unit = display_unit;
        this.display_mul = display_mul;
        if (min_pos == null) {
            orientation_dir = new Vector3(1, 0, 0);
        }
        else {
            orientation_dir = (max_pos.position - min_pos.position).normalized;
        }
    }

    // Update is called once per frame, and after val updated
    void Update() {
        Vector3 lp = meter.transform.localPosition;
        lp.x = 0.05f - val * 0.1f;
        meter.transform.localPosition = lp;
        forceMap();

        /*
        float magnitude = 0.05f - val * 0.1f;
        meter.transform.localPosition = Quaternion.Euler(0, 2.0f, 0) * meter.transform.forward * magnitude;
        forceMap();
        */
    }

    public void forceMap() {
        //map = min_map+(max_map-min_map)*val;
        map = response_power > 1 ? mapSharp() : mapLinear();
    }

    public void Reset() {
        val = default_val;
        Update();
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
}
