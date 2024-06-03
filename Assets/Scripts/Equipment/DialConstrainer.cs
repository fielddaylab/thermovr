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
                float margin = 0.05f;
                if (m_ConstrainType == ConstrainType.Max) { margin *= -1; }

                var currVal = (float)ThermoPresent.Instance.get_state_var(m_ConstrainingVar);
                float targetMap = (float)((currVal - ThermoMath.v_min) / (ThermoMath.v_max - ThermoMath.v_min));
                float marginMap = (float)((currVal + margin - ThermoMath.v_min) / (ThermoMath.v_max - ThermoMath.v_min));
                float targetVal = m_Dial.MapToDialVal(targetMap);
                float marginVal = m_Dial.MapToDialVal(marginMap);
                m_Dial.SetConstraint(targetVal, m_ConstrainType, marginVal - targetVal);
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