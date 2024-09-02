using BeauRoutine;
using BeauUtil;
using BeauUtil.Debugger;
using BeauUtil.Extensions;
using System.Collections;
using System.Collections.Generic;
using ThermoVR.UI;
using UnityEngine;

namespace ThermoVR
{
    public class GameMgr : Singleton<GameMgr>
    {
        public bool AudioEnabled = false;

        public bool IsAlphaRelease = true; // temp solution to managing alpha release channel
        public bool IsDesktop = false;

        [SerializeField] private World m_world;
        [SerializeField] private ThermoPresent m_thermo_present;

        [SerializeField] private SaveSystem m_SaveSystem = null;
        [SerializeField] private UILoading m_UILoading;
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

            EventMgr.Events.Register(GameEvents.TryNewName, OnTryNewName, this);

            EventMgr.Events.Dispatch(GameEvents.StartGame);

            EventMgr.Events.Dispatch(GameEvents.TryNewName);

            EventMgr.Events.Dispatch(GameEvents.InitialLoadComplete);
        }

        private void FixedUpdate() {
            m_world.ManualFixedUpdate();
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

            EventMgr.Events.Dispatch(GameEvents.NewNameGenerated, inName);
            m_ProfileName = inName;

            EventMgr.Events.Dispatch(GameEvents.StartSession);
        }

        private void OnNewNameFail(OGD.Core.Error error)
        {
            Debug.Log("[Analytics] New failed.");

            Log.Error("[Game] Generating new player id failed: {0}", error.Msg);

            EventMgr.Events.Dispatch(GameEvents.StartSession);
        }

        #endregion // New Game
    }
}