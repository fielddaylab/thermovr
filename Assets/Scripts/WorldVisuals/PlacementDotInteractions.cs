using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using ThermoVR.Controls;
using ThermoVR.Tools;

namespace ThermoVR
{
    public class PlacementDotInteractions : MonoBehaviour
    {
        public GameObject state_dot;
        public GameObject placement_dot;
        [HideInInspector] public Vector3 placement_thermo;
        [HideInInspector] public bool placement_thermo_reasonable;

        [SerializeField] private new MeshCollider collider;

        public GameObject graph;
        private Touchable graph_touchable;

        public void Init()
        {
            state_dot = GameObject.Find("gstate");
            placement_dot = GameObject.Find("tstate");
            placement_dot.GetComponent<Renderer>().enabled = false;
            placement_thermo_reasonable = false;

            graph = GameObject.Find("Graph");
            graph_touchable = graph.GetComponent<Touchable>();
        }

        public void AssignMesh(Mesh mesh)
        {
            collider.sharedMesh = mesh;
        }

        public void BeginInteract(Hand handType)
        {
            graph_touchable.SetGrabbed(true, handType);

            state_dot.GetComponent<Renderer>().enabled = false;
            placement_dot.GetComponent<Renderer>().enabled = true;
        }

        public void CancelInteract()
        {
            placement_thermo_reasonable = false;
        }

        public void FinishInteract(Hand handType)
        {
            graph_touchable.SetGrabbed(false, handType);

            placement_dot.GetComponent<Renderer>().enabled = false;
            state_dot.GetComponent<Renderer>().enabled = true;

            if (placement_thermo_reasonable)
            {
                ToolMgr.Instance.DeactivateAllTools();
                World.Instance.WarpPVT(
                    placement_thermo.y,
                    placement_thermo.x,
                    placement_thermo.z
                    );
            }
        }

        public void ContinueInteract(Vector3 interactPos)
        {
            Vector3 localspace = graph.transform.InverseTransformPoint(interactPos);
            Vector3 correctedspace = new Vector3(localspace.z, localspace.y, -localspace.x) * 4.0f; //rotate 90, mul by 4 (inverse transform of gmodel)

            //Vector3 thermoguess = thermo.guessPlot(ThermoMath.t_neutral, correctedspace.y, correctedspace.x);
            Vector3 thermoguess = ThermoPresent.Instance.guessMeshPlot(correctedspace.x, correctedspace.y, correctedspace.z);
            Vector3 localguess = ThermoPresent.Instance.plot(thermoguess.y, thermoguess.x, thermoguess.z); //note swizzle!

            if (MathUtility.floatNumeric(localguess.x) && MathUtility.floatNumeric(localguess.y) && MathUtility.floatNumeric(localguess.z))
            {
                placement_dot.transform.localPosition = localguess;
                placement_thermo = thermoguess;
                placement_thermo_reasonable = true;
            }
            else
            {
                placement_thermo_reasonable = false;
            }
        }
    }
}