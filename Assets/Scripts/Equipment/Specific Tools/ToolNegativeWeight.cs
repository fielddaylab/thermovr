using System.Collections;
using System.Collections.Generic;
using ThermoVR.Physics;
using UnityEngine;

namespace ThermoVR.Tools
{
    public class ToolNegativeWeight : Tool
    {
        #region Inspector

        [SerializeField] private FixedAnchor m_ElasticAnchor;

        #endregion // Inspector

        #region Tool

        protected override IEnumerator ActivationRoutine() {
            // disable anchor until activation completed
            m_ElasticAnchor.enabled = false;

            gameObject.SetActive(true);

            m_ElasticAnchor.enabled = true;
            yield return null;
        }

        protected override IEnumerator DeactivationRoutine() {
            // disable anchor for deactivation
            m_ElasticAnchor.enabled = false;

            gameObject.SetActive(false);
            m_ElasticAnchor.enabled = true;
            yield return null;
        }

        #endregion // Tool
    }
}