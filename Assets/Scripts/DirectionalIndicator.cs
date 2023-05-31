using System.Collections;
using System.Collections.Generic;
using ThermoVR;
using UnityEngine;

/**
 * Simple class to handle the arrow indicators for heat in/out of the system.
 * Not sure it will be useful for whatever we use in the future, but good enough for now.
 * -Luke
 **/
public class DirectionalIndicator : MonoBehaviour
{
    private enum FLOW_DIRECTIONS { OUT, IN };
    const float BASE_FRAME_RATE = 1f / 120f; // base increment will take 120 frames if coefficient = 1.
    const float BASE_RATE_SCALING = 1f / 1000f; // further, we drop based on observations that the flow coefficients are often in the thousands.

    public Material indicator_material;
    public Transform T_left;
    public Transform T_right;
    [System.NonSerialized]
    public bool running;
    public float flow_coefficient;
    private float current_offset;
    private FLOW_DIRECTIONS flow_direction;

    // Start is called before the first frame update
    void Start() {
        current_offset = 0.0f;
        flow_direction = FLOW_DIRECTIONS.OUT;
        flow_coefficient = 1.0f;
        running = true;
        this.Stop(); // we start the game stopped.
    }

    /*
     * Function to stop the indicator, if heat is no longer actually flowing.
     */
    public void Stop() {
        if (running) {
            running = false;
            MeshRenderer[] mr_list = gameObject.GetComponentsInChildren<MeshRenderer>();
            foreach (MeshRenderer mr in mr_list) {
                mr.enabled = false;
            }
        }
    }

    /*
     * Function to show the indicator, when heat starts to flow.
     */
    public void Go(bool flow_outward) {
        // If we weren't running already, then start.
        if (!running) {
            running = true;
            MeshRenderer[] mr_list = gameObject.GetComponentsInChildren<MeshRenderer>();
            foreach (MeshRenderer mr in mr_list) {
                mr.enabled = true;
            }
        }
        // If we aren't pointing in the same direction as flow, fix it.
        FLOW_DIRECTIONS correct_direction = flow_outward ? FLOW_DIRECTIONS.OUT : FLOW_DIRECTIONS.IN;
        if (flow_direction != correct_direction) {
            flow_direction = correct_direction;
            float y_rot = flow_outward ? 90f : -90f;
            T_left.rotation = Quaternion.Euler(90f, y_rot, 0f);
            T_right.rotation = Quaternion.Euler(90f, y_rot, 0f);
        }
    }

    // Update is called once per frame
    void Update() {
        if (running) {
            float increment = flow_coefficient * BASE_FRAME_RATE * BASE_RATE_SCALING; // / BASE_FRAME_PER_CYCLE;
            current_offset += increment;
            // Cycle back to 0 when we reach offset of 10, just so we avoid getting huge values over time.
            if (Mathf.Abs(current_offset) >= 10.0f) {
                current_offset = 0.0f;
            }
            indicator_material.SetTextureOffset("_MainTex", new Vector2(current_offset, 0.0f));
        }
    }

    public void SetFlow(double coefficient) {
        flow_coefficient = (float)coefficient;
    }

}
