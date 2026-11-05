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

        // False until InitRoutine() has finished; gates FixedUpdate.
        private bool m_initComplete = false;

        protected override void Awake() {
            base.Awake();
        }

        /// <summary>
        /// 0..1 fraction of the initial load that has completed. Monotonic non-decreasing,
        /// computed from completed work units (never from a timer), exactly 1f when done.
        /// </summary>
        public static float LoadProgress {
            get { return LoadProgressTracker.Progress; }
        }

        private void Start() {
            AudioEnabled = false;

            // The init chain used to run fully synchronously inside this single frame, which
            // meant no Update() ran and no loading feedback could animate. It is now a
            // coroutine that yields between phases.
            StartCoroutine(InitRoutine());
        }

        /// <summary>
        /// Initial load, broken into phases separated by `yield return null` so frames
        /// render during the load. Progress is published through LoadProgressTracker.
        /// </summary>
        private IEnumerator InitRoutine() {
            System.Diagnostics.Stopwatch sw_total = System.Diagnostics.Stopwatch.StartNew();
            System.Diagnostics.Stopwatch sw_phase = new System.Diagnostics.Stopwatch();

            // Fixed denominator, established before any unit is completed, so the reported
            // fraction can never move backwards: ThermoPresent.Init phases + chunks,
            // plus world.init, thermo.reset and the final dispatch phase.
            const int GAMEMGR_UNITS = 3; // world.init, thermo_present.reset, dispatch
            LoadProgressTracker.Begin(m_thermo_present.CountInitWorkUnits() + GAMEMGR_UNITS);

            IEnumerator thermo_init = m_thermo_present.InitRoutine();
            while (thermo_init.MoveNext()) {
                yield return thermo_init.Current;
            }

            yield return null;

            sw_phase.Restart();
            m_world.Init();
            Debug.Log("[LoadTiming] phase=world.init ms=" + sw_phase.ElapsedMilliseconds);
            LoadProgressTracker.Complete();
            sw_phase.Stop();
            yield return null;

            sw_phase.Restart();
            m_thermo_present.Reset();
            Debug.Log("[LoadTiming] phase=thermo_present.reset ms=" + sw_phase.ElapsedMilliseconds);
            LoadProgressTracker.Complete();
            sw_phase.Stop();
            yield return null;

            sw_phase.Restart();

            m_initComplete = true;

            AudioEnabled = true;

            m_ProfileName = string.Empty;

            EventMgr.Events.Register(GameEvents.TryNewName, OnTryNewName, this);

            EventMgr.Events.Dispatch(GameEvents.StartGame);

            EventMgr.Events.Dispatch(GameEvents.TryNewName);

            EventMgr.Events.Dispatch(GameEvents.InitialLoadComplete);

            Debug.Log("[LoadTiming] phase=dispatch ms=" + sw_phase.ElapsedMilliseconds);
            sw_phase.Stop();

            // Every reserved work unit is now done; pin the fraction to exactly 1.
            LoadProgressTracker.Finish();

            Debug.Log("[LoadTiming] total ms=" + sw_total.ElapsedMilliseconds);
            sw_total.Stop();
        }

        private void FixedUpdate() {
            // Init used to finish inside Start(), i.e. before the first FixedUpdate. Now
            // that it yields across frames, FixedUpdate can fire before World.Init() has
            // run, so it must be held off until the chain reports completion.
            if (!m_initComplete) {
                return;
            }

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