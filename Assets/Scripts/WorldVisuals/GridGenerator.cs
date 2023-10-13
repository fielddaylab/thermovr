using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR
{
    public class GridGenerator : MonoBehaviour
    {
        [SerializeField] private GameObject gridMeshPrefab;
        [SerializeField] private int numXSegments, numZSegments = 8;
        [SerializeField] private float step = 1;

        private List<GameObject> previousGenerations;

        [ContextMenu("Generate 3D Grid")]
        private void Generate3DGrid() {
            if (previousGenerations == null) {
                previousGenerations = new List<GameObject>();
            }

            for (int i = 0; i < previousGenerations.Count; i++) {
                DestroyImmediate(previousGenerations[i]);
            }

            previousGenerations.Clear();

            float step = 1;
            int numRows = numXSegments;
            int numCols = numZSegments;

            for (int r = 0; r < numRows; r++) {
                GameObject gridPlane = Instantiate(gridMeshPrefab, this.transform);
                gridPlane.transform.Translate(new Vector3(-step * r * this.transform.localScale.x, 0, 0), Space.World);
                gridPlane.GetComponent<GridMesh>().GenerateAll();
                previousGenerations.Add(gridPlane);
            }

            for (int c = 0; c < numCols; c++) {
                GameObject gridPlane = Instantiate(gridMeshPrefab, this.transform);
                gridPlane.transform.Rotate(new Vector3(-90f, 0, 0));
                gridPlane.transform.Translate(new Vector3(-step * numCols / 2 * this.transform.localScale.x, 0, -step * numCols/2 * this.transform.localScale.z + step * c * this.transform.localScale.x), Space.World);
                gridPlane.GetComponent<GridMesh>().GenerateAll();
                previousGenerations.Add(gridPlane);
            }
        }
    }
}