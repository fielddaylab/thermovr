using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR {
    public interface IEvaluable
    {
        public abstract bool IsCorrect();

        public abstract void ResetState();

        public abstract void HandleEvaluation(bool correct);
    }
}

