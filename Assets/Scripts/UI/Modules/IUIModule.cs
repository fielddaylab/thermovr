using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.UI.Interfaces
{
    public interface IUIModule
    {
        UIID ID {
            get;
        }

        public abstract void Init();

        public abstract void Open();

        public abstract void Close();
    }

}