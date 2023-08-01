using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Lab
{

    public abstract class Evaluable : MonoBehaviour, IEvaluable
    {
        protected bool m_evaluated;

        #region IEvaluable

        public abstract bool IsCorrect();

        public abstract bool AnswerSelected();

        public virtual void ResetState() {
            m_evaluated = false;
        }

        public abstract void HandleEvaluation(bool correct);

        #endregion // IEvaluable
    }
}
