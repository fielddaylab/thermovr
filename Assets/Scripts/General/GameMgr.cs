using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR
{
    public class GameMgr : MonoBehaviour
    {
        [SerializeField] private World m_world;

        private void Start() {
            m_world.Init();
        }

        private void FixedUpdate() {
            m_world.ManualFixedUpdate();
        }
    }
}