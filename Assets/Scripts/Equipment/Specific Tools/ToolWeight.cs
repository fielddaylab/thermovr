using System.Collections;
using System.Collections.Generic;
using ThermoVR.Physics;
using UnityEngine;

namespace ThermoVR.Tools
{
    public class ToolWeight : Tool
    {
        #region Inspector

        [SerializeField] private FixedAnchor m_HydraulicPressAnchor;

        #endregion // Inspector

        #region Tool

        protected override IEnumerator ActivationRoutine() {
            // disable anchor until activation completed
            m_HydraulicPressAnchor.enabled = false;

            gameObject.SetActive(true);

            m_HydraulicPressAnchor.enabled = true;
            yield return null;
        }

        protected override IEnumerator DeactivationRoutine() {
            // disable anchor until deactivation completed
            m_HydraulicPressAnchor.enabled = false;

            gameObject.SetActive(false);

            m_HydraulicPressAnchor.enabled = true;
            yield return null;
        }

        #endregion // Tool
    }
}