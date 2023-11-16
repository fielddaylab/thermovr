using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Dials
{
    [RequireComponent(typeof(Dial))]
    public class DialConstrainer : MonoBehaviour
    {
        [SerializeField] private ConstrainType m_ConstrainType;
        [SerializeField] private Dial m_ConstrainingDial; // If the constrainer is a different dial value

        private Dial m_Dial; // The dial being constrained

        private void OnEnable() {
            m_Dial = this.GetComponent<Dial>();

            if (m_ConstrainingDial) {
                m_ConstrainingDial.DialMoved.AddListener(HandleDialMoved);
            }
        }

        private void OnDisable() {
            if (m_ConstrainingDial) {
                m_ConstrainingDial.DialMoved.RemoveListener(HandleDialMoved);
            }
        }

        #region Handlers

        private void HandleDialMoved() {
            if (m_Dial) {
                m_Dial.SetConstraint(m_ConstrainingDial.get_val(), m_ConstrainType);
            }
        }

        #endregion // Handlers
    }

}