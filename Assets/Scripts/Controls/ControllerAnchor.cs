using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.XR.Interaction.Toolkit;

namespace ThermoVR.Controls
{
    public class ControllerAnchor : MonoBehaviour
    {
        public GameObject obj;
        public GameObject actualhand;
        public XRRayInteractor ray;
        [HideInInspector] public Vector3 pos;
        [HideInInspector] public Vector3 vel;
        public SkinnedMeshRenderer meshrenderer;

        public void Init(Material[] hand_emptys) {
            pos = this.transform.position;
            vel = new Vector3(0, 0, 0);
            meshrenderer.materials = hand_emptys;
        }
    }
}