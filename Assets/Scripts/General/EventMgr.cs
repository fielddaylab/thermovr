using BeauUtil;
using BeauUtil.Extensions;
using System.Collections;
using System.Collections.Generic;
using ThermoVR;
using UnityEngine;

namespace ThermoVR {
    [DefaultExecutionOrder(-5000)]
    public class EventMgr : MonoBehaviour
    {
        public static EventMgr I;

        private readonly EventDispatcher<EvtArgs> m_EventDispatcher = new EventDispatcher<EvtArgs>();

        protected void Awake()
        {
            if (I == null)
            {
                I = this;
                DontDestroyOnLoad(this.gameObject);
            }
            else if (I != this)
            {
                Destroy(this.gameObject);
                return;
            }
        }

        private void LateUpdate()
        {
            m_EventDispatcher.Flush();
        }


        /// <summary>
        /// Global game event dispatcher.
        /// </summary>
        static public EventDispatcher<EvtArgs> Events
        {
            get { return EventMgr.I?.m_EventDispatcher; }
        }
    }
}
