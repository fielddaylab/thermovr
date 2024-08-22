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

        private readonly EventDispatcher<object> m_EventDispatcher = new EventDispatcher<object>();

        protected void Awake()
        {
            if (I == null)
            {
                I = this;
            }
            else if (I != this)
            {
                Destroy(this.gameObject);
                return;
            }
        }

        private void LateUpdate()
        {
            m_EventDispatcher.FlushQueue();
        }


        /// <summary>
        /// Global game event dispatcher.
        /// </summary>
        static public EventDispatcher<object> Events
        {
            get { return EventMgr.I?.m_EventDispatcher; }
        }
    }
}
