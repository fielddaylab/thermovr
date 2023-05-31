using BeauUtil;
using BeauUtil.Extensions;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR
{
    public class GameMgr : Singleton<GameMgr>
    {
        private readonly EventDispatcher<object> m_EventDispatcher = new EventDispatcher<object>();

        [SerializeField] private World m_world;

        protected override void Awake() {
            base.Awake();
        }

        private void Start() {
            m_world.Init();
        }

        private void FixedUpdate() {
            m_world.ManualFixedUpdate();
        }

        private void LateUpdate() {
            m_EventDispatcher.FlushQueue();
        }

        /// <summary>
        /// Gloabl game event dispatcher.
        /// </summary>
        static public EventDispatcher<object> Events {
            get { return GameMgr.I?.m_EventDispatcher; }
        }
    }
}