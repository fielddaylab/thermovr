using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Tools
{
    [RequireComponent(typeof(Collider))]
    [RequireComponent(typeof(MeshRenderer))]
    public class CritBarrierMgr : MonoBehaviour
    {
        public static CritBarrierMgr Instance;

        [SerializeField] private Collider m_graphBall;
        [SerializeField] private Material m_neutralMat;
        [SerializeField] private Material m_activeMat;

        private Collider m_barrier;
        private MeshRenderer m_meshRenderer;

        private void Awake()
        {
            if (Instance == null)
            {
                Instance = this;
            }
            else if (Instance != this)
            {
                Destroy(this.gameObject);
                return;
            }

            m_barrier = this.GetComponent<Collider>();
            m_meshRenderer = this.GetComponent<MeshRenderer>();
        }

        private void OnTriggerEnter(Collider other)
        {
            if (other == m_graphBall)
            {
                Debug.Log("[CritBarrierMgr] graph ball is in barrier zone!");
                EventMgr.Events.Dispatch(GameEvents.RestoreLastSimState);
                SetMat(m_activeMat);
            }
        }

        private void OnTriggerStay(Collider other)
        {
            if (other == m_graphBall)
            {
                Debug.Log("[CritBarrierMgr] graph ball is in barrier zone!");
                EventMgr.Events.Dispatch(GameEvents.RestoreLastSimState);
            }
        }

        private void OnTriggerExit(Collider other)
        {
            if (other == m_graphBall)
            {
                Debug.Log("[CritBarrierMgr] graph ball left the barrier zone!");
                SetMat(m_neutralMat);
            }
        }

        private void SetMat(Material newMat)
        {
            var mats = m_meshRenderer.sharedMaterials;
            mats[0] = newMat;
            m_meshRenderer.sharedMaterials = mats;
        }
    }
}