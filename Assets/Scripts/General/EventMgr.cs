using BeauUtil;
using BeauUtil.Extensions;
using System.Collections;
using System.Collections.Generic;
using ThermoVR;
using UnityEngine;

public class EventMgr : Singleton<EventMgr>
{
    private readonly EventDispatcher<object> m_EventDispatcher = new EventDispatcher<object>();

    protected override void Awake()
    {
        base.Awake();
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
