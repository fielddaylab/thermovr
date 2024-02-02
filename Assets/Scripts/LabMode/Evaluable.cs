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

        public virtual bool HasBeenEvaluated() {
            return m_evaluated;
        }

        public virtual void ResetState() {
            m_evaluated = false;
        }

        public abstract void HandleEvaluation(bool correct);

        public void LoadCompleted(bool completed)
        {
            if (completed)
            {
                m_evaluated = true;
            }
        }

        #endregion // IEvaluable
    }
}
