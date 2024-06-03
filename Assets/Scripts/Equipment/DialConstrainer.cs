using System.Collections;
using System.Collections.Generic;
using ThermoVR.State;
using UnityEngine;

namespace ThermoVR.Dials
{
    [RequireComponent(typeof(Dial))]
    public class DialConstrainer : MonoBehaviour
    {
        [SerializeField] private ConstrainType m_ConstrainType;

        [Header("Dial Constraint")]
        [SerializeField] private Dial m_ConstrainingDial; // If the constrainer is a different dial value

        [Header("Sim Constraint")]
        [SerializeField] private VarID m_ConstrainingVar; // If the constraint is a simulation variable

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

        private void Update()
        {
            if (!m_ConstrainingDial)
            {
                var currVal = (float)ThermoPresent.Instance.get_state_var(m_ConstrainingVar);
                float targetMap = (float)((currVal - ThermoMath.v_min) / (ThermoMath.v_max - ThermoMath.v_min));
                float dialVal = m_Dial.MapToDialVal(targetMap);
                Debug.Log("[Constraint] currVal: " + currVal + " || dialVal: " + dialVal);
                m_Dial.SetConstraint(dialVal, m_ConstrainType, 0.01f);
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