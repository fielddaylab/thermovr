using System.Collections;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using UnityEngine;
using UnityEngine.Scripting;

/**
 * Class to detect when an object gets touched and/or grabbed by the hand controller,
 * and to mark itself as touched, for later processing in "World."
 **/
public class Touchable : MonoBehaviour
{
    protected GameObject lhand;
    protected GameObject rhand;
    [System.NonSerialized]
    public bool grabbed = false;
    protected Collider lhand_c;
    protected Collider rhand_c;
    [System.NonSerialized]
    public Transform og_parent;
    private Lightable[] lightables;
    public bool has_lightables = true;

    void Awake() {
        lhand = GameObject.Find("LeftControllerAnchor");
        rhand = GameObject.Find("RightControllerAnchor");
        //lhand_c = lhand.GetComponentsInChildren<Collider>();
        //rhand_c = rhand.GetComponentsInChildren<Collider>();
        if (lhand != null ) {
            lhand_c = lhand.GetComponent<SphereCollider>();
        }
        if (rhand != null) {
            rhand_c = rhand.GetComponent<SphereCollider>();
        }
        og_parent = gameObject.transform.parent;

        if (has_lightables) {
            lightables = gameObject.GetComponentsInChildren<Lightable>();
        }
        else {
            lightables = null;
        }
    }

    private void OnDisable()
    {
        rtouch = ltouch = false;
    }

    /*
     * Check if the given collider is in the list of colliders found on
     * the left/right hand.
     */
    private bool cInHandList(Collider c, bool left) {
        Collider hc = left ? lhand_c : rhand_c;
        if (c == hc) return true; //big volume grab collider
        return false;
    }

    //[System.NonSerialized]
    [HideInInspector] public bool ltouch = false;
    //[System.NonSerialized]
    [HideInInspector] public bool rtouch = false;
    //[System.NonSerialized]
    [HideInInspector] public bool touch = false;

    public void deltaTouch(Collider c, bool delta) {
        //bool prev_l_touch = ltouch;
        //bool prev_r_touch = rtouch;
        //bool in_either_hand = cInHandList(c, true) || cInHandList(c, false);
        if (cInHandList(c, true)) ltouch = delta;
        if (cInHandList(c, false)) rtouch = delta;

        touch = (ltouch || rtouch);

        /*
        string name = this.gameObject.name;
        if (touch && in_either_hand)
        {
            Debug.Log("[Touch] " + name + " started touching " + c.name + ": l--" + ltouch + "|| r--" + rtouch);
            Debug.Log("[Touch] prev was : l--" + prev_l_touch + "|| r--" + prev_r_touch);
        }
        else if (in_either_hand)
        {
            Debug.Log("[Touch] " + name + " stopped touching " + c.name + ": l--" + ltouch + "|| r--" + rtouch);
            Debug.Log("[Touch] prev was : l--" + prev_l_touch + "|| r--" + prev_r_touch);
        }
        */

        if (lightables != null) {
            foreach (Lightable light in lightables) {
                light.SetLit(touch);
            }
        }
    }

    /// <summary>
    /// 
    /// </summary>
    /// <param name="queryingLeft">true if querying the left finger, false if the right finger</param>
    /// <returns></returns>
    public void SetFingerTouches(ref bool ltouch, ref bool rtouch) {
        if (this.ltouch) { ltouch = true; }
        if (this.rtouch) { rtouch = true; }
    }

    protected void OnTriggerEnter(Collider c) {
        deltaTouch(c, true);
    }

    protected void OnTriggerExit(Collider c) {
        deltaTouch(c, false);
    }

}

