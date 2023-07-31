using System.Collections;
using System.Collections.Generic;
using ThermoVR.UI.Interfaces;
using UnityEngine;

namespace ThermoVR.UI
{
    // [RequireComponent(typeof(Canvas))]
    public class UIModule : MonoBehaviour, IUIModule
    {
        [SerializeField] private UIID m_id;

        public bool AlwaysAvailable = false;

        protected UIHub m_hub;

        public UIID ID {
            get { return m_id; }
        }

        public virtual void Init() {
            m_hub = this.GetComponent<UIHub>();
        }

        public virtual void Close() {
            this.gameObject.SetActive(false);
        }

        public virtual void Open() { 
            this.gameObject.SetActive(true);
        }
    }
}