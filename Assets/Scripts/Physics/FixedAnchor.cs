using BeauUtil;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Physics
{
    public class FixedAnchor : MonoBehaviour
    {
        [Required] [SerializeField] private GameObject m_AnchorTo;

        private void Update() {
            this.transform.position = m_AnchorTo.transform.position;
        }
    }
}
