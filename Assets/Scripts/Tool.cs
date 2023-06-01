/*
DOCUMENTATION- phil, 12/16/19 [intended to be a description as of a point in time, NOT nec a prescription for how it should be evolved- feel free to uproot]

The various tools which can be variably engaged to the container, thrown, or stored.
*/

using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using TMPro;
using ThermoVR.Dials;

namespace ThermoVR.Tools {
    public struct VolumeStop
    {
        public double Volume;
        public Tool Source;

        public VolumeStop(double volume, Tool source) {
            Volume = volume;
            Source = source;
        }
    }

    public class Tool : MonoBehaviour
    {
        [System.NonSerialized]
        public bool engaged = false;
        [System.NonSerialized]
        public bool stored = false;
        [System.NonSerialized]
        public Touchable touchable;
        [System.NonSerialized]
        public BoxCollider boxcollider;
        [System.NonSerialized]
        public Rigidbody rigidbody;
        [System.NonSerialized]
        public float t_free = 0.0f;

        public GameObject mesh;

        public GameObject storage;
        [System.NonSerialized]
        public float default_storage_scale;
        [System.NonSerialized]
        public Ghost storage_ghost;
        [System.NonSerialized]
        public Touchable storage_touchable;
        [System.NonSerialized]
        public GameObject storage_obj;
        [System.NonSerialized]
        public MeshRenderer storage_meshrenderer;
        [System.NonSerialized]
        public GameObject storage_snap;

        public GameObject active;
        [System.NonSerialized]
        public Ghost active_ghost;
        [System.NonSerialized]
        public Touchable active_touchable;
        [System.NonSerialized]
        public GameObject active_available;
        [System.NonSerialized]
        public MeshRenderer active_available_meshrenderer;
        [System.NonSerialized]
        public GameObject active_snap;

        private float val = 0.0f; // value of the tool

        public GameObject text;
        //[HideInInspector]
        //public GameObject textv;
        

        public GameObject textl;
        [System.NonSerialized]
        public TextMeshPro textl_tmpro;
        [System.NonSerialized]
        public MeshRenderer textl_meshrenderer;

        public GameObject textd;
        [System.NonSerialized]
        public TextMeshPro textd_tmpro;
        [System.NonSerialized]
        public MeshRenderer textd_meshrenderer;

        public GameObject textn;

        [System.NonSerialized]
        public Fadable text_fadable;
        [System.NonSerialized]
        public bool disabled;

        [System.NonSerialized]
        public string display_unit = "";
        [System.NonSerialized]
        public float display_mul = 1.0f; //multiplied with map before displaying with display_unit

        [System.NonSerialized]
        public AudioSource audioS;

        private bool examined;

        public void Init(string unit, float mul = 1) {
            this.display_unit = unit;
            this.display_mul = mul;

            Setup();
        }

        private void Setup() {
            touchable = gameObject.GetComponent<Touchable>();
            boxcollider = gameObject.GetComponent<BoxCollider>();
            rigidbody = gameObject.GetComponent<Rigidbody>();

            if (storage != null) {
                default_storage_scale = storage.transform.localScale.x; //could grab any dimension
                storage_ghost = storage.GetComponent<Ghost>();
                storage_ghost.set_tool(this);
                storage_touchable = storage.GetComponent<Touchable>();
                storage_obj = storage_ghost.obj;
                storage_meshrenderer = storage_obj.GetComponent<MeshRenderer>();
            }

            if (active != null) {
                active_ghost = active.GetComponent<Ghost>();
                active_ghost.set_tool(this);
                active_touchable = active.GetComponent<Touchable>();
                active_available = active_ghost.obj;
                active_available_meshrenderer = active_available.GetComponent<MeshRenderer>();
            }

            // TODO: separate textv from dial, move to tool
            // textv = dial_dial.textv.gameObject;
            //textv_tmpro = textv.GetComponent<TextMeshPro>();
            //textv_meshrenderer = textv.GetComponent<MeshRenderer>();
            textl_tmpro = textl.GetComponent<TextMeshPro>();
            textl_meshrenderer = textl.GetComponent<MeshRenderer>();
            text_fadable = GetComponent<Fadable>();
            if (textd) {
                textd_tmpro = textd.GetComponent<TextMeshPro>();
                textd_meshrenderer = textd.GetComponent<MeshRenderer>();
            }

            audioS = gameObject.GetComponent<AudioSource>();

            disabled = false;
            enabled = true;
            examined = false;
        }

        // Update is called once per frame
        void Update() {
            t_free += Time.deltaTime;

            if (text_fadable) text_fadable.set_factive(touchable.touch || examined);

            // If tool has been sitting unused and unmoved for long enough, start it moving back to storage.
            if (!engaged && !stored && !touchable.grabbed && t_free > 1.0f && storage != null) {
                rigidbody.isKinematic = true;
                gameObject.transform.position = Vector3.Lerp(gameObject.transform.position, storage.transform.position, 0.1f);
                gameObject.transform.localRotation = Quaternion.Lerp(gameObject.transform.localRotation, storage.transform.localRotation, 0.1f);
                if (Vector3.Distance(gameObject.transform.position, storage.transform.position) < 0.01f) {
                    gameObject.transform.SetParent(storage.transform);
                    stored = true;
                    gameObject.transform.localPosition = new Vector3(0f, 0f, 0f);
                    gameObject.transform.localScale = new Vector3(1f, 1f, 1f);
                    gameObject.transform.localRotation = Quaternion.identity;
                    float v = storage.transform.localScale.x; //can grab any dimension
                    Vector3 invscale = new Vector3(1f / v, 1f / v, 1f / v);
                    text.transform.localScale = invscale;
                }
            }
            else if (engaged || stored || touchable.grabbed)
                t_free = 0.0f;
        }

        public void update_val(float new_val, Dial dial) {
            val = new_val;

            update_tool_text(dial);
        }

        public void update_text(Dial dial) {
            update_tool_text(dial);
        }

        public float get_val() {
            return val;
        }

        /// <summary>
        /// Moves the engaged tool position
        /// </summary>
        /// <param name="dv">the dial's value difference</param>
        public void move(float dv) {
            active.transform.position += new Vector3(0, dv, 0);
        }

        public void set_examined(bool examined) {
            this.examined = true;
        }

        private void update_tool_text(Dial dial) {
            // find dial
            // if (!dial_dial) return;
            // do not apply to insulator
            // if (t == tool_insulator) return;

            string txt = string.Format("{0:3} " + display_unit, (float)(dial.map * display_mul));
            //t.textv_tmpro.SetText(txt);
            //if(t.textd_tmpro) t.textd_tmpro.SetText(txt);
            dial.textv_tmpro.SetText("{0:3} " + display_unit, (float)(dial.map * display_mul));
            if (textd_tmpro) textd_tmpro.SetText("{0:3} " + display_unit, (float)(dial.map * display_mul));
        }

    }
}

