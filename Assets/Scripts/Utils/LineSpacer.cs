using System.Collections;
using System.Collections.Generic;
using UnityEditor;
using UnityEngine;

public class LineSpacer : MonoBehaviour
{
    [SerializeField] private GameObject[] m_ToSpace;
    [SerializeField] private Transform m_StartPos;
    [SerializeField] private float m_SpacingPerChunk = 75;


    [ContextMenu("ApplySpacing")]
    public void ApplySpacing() {
        float currY = m_StartPos.position.y;

        for (int i = 0; i < m_ToSpace.Length; i++) {
            m_ToSpace[i].transform.position = new Vector3(m_StartPos.position.x, currY, m_StartPos.position.z);

            currY -= m_SpacingPerChunk;
        }
    }
}