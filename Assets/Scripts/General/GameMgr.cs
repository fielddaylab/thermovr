using BeauRoutine;
using BeauUtil;
using BeauUtil.Debugger;
using BeauUtil.Extensions;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR
{
    public class GameMgr : Singleton<GameMgr>
    {
        public bool AudioEnabled = false;

        public bool IsAlphaRelease = true; // temp solution to managing alpha release channel
        public bool IsDesktop = false;

        private readonly EventDispatcher<object> m_EventDispatcher = new EventDispatcher<object>();

        [SerializeField] private World m_world;
        [SerializeField] private ThermoPresent m_thermo_present;

        [SerializeField] private SaveSystem m_SaveSystem = null;
        private string m_ProfileName;

        protected override void Awake() {
            base.Awake();
        }

        private void Start() {
            AudioEnabled = false;

            m_thermo_present.Init();
            m_world.Init();
            m_thermo_present.Reset();

            AudioEnabled = true;

            m_ProfileName = string.Empty;

            Events.Register(GameEvents.TryNewName, OnTryNewName, this);

            Events.Dispatch(GameEvents.InitialLoadComplete);
        }

        private void FixedUpdate() {
            m_world.ManualFixedUpdate();
        }

        private void LateUpdate() {
            m_EventDispatcher.FlushQueue();
        }

        /// <summary>
        /// Global game event dispatcher.
        /// </summary>
        static public EventDispatcher<object> Events {
            get { return GameMgr.I?.m_EventDispatcher; }
        }


        #region New Game

        private void OnTryNewName()
        {
            if (m_ProfileName.Equals(string.Empty))
            {
                Debug.Log("[Analytics] New name try...");
                OGD.Player.NewId(OnNewNameSuccess, OnNewNameFail);
            }
        }

        private void OnNewNameSuccess(string inName)
        {
            Debug.Log("[Analytics] New name success! " + inName);

            GameMgr.Events.Dispatch(GameEvents.NewNameGenerated, inName);
            m_ProfileName = inName;

            GameMgr.Events.Dispatch(GameEvents.StartSession);
            GameMgr.Events.Dispatch(GameEvents.StartGame);
        }

        private void OnNewNameFail(OGD.Core.Error error)
        {
            Debug.Log("[Analytics] New failed.");

            Log.Error("[Game] Generating new player id failed: {0}", error.Msg);

            GameMgr.Events.Dispatch(GameEvents.StartSession);
            GameMgr.Events.Dispatch(GameEvents.StartGame);
        }

        #endregion // New Game
    }
}