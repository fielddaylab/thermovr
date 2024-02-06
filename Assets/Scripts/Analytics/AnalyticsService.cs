#if UNITY_EDITOR || DEVELOPMENT_BUILD
    #define DEVELOPMENT
#endif // UNITY_EDITOR || DEVELOPMENT_BUILD


using BeauUtil;
using System;
using System.Collections.Generic;
using UnityEngine;
using FieldDay;
using BeauUtil.Tags;
using BeauPools;

namespace ThermoVR.Analytics
{

    public partial class AnalyticsService : MonoBehaviour
    {
        private enum GamePlatform { 
            VR,
            WEB
        }


        #region Inspector

        [SerializeField, Required] private string m_AppId = "THERMOVR";
        [SerializeField, Required] private string m_AppVersion = "1.0";
        [SerializeField] private FirebaseConsts m_Firebase = default(FirebaseConsts);

        #endregion // Inspector

        #region Logging Enums & Structs



        #endregion // Logging Structs

        #region Logging Variables

        private OGDLog m_Log;

        private GamePlatform m_Platform;

        [NonSerialized] private bool m_Debug;


        #endregion // Logging Variables

        private void Awake() {
            Initialize();
        }

        #region IService

        protected void Initialize()
        {
            // General Events
            //Game.Events.Register(GameEvents.StoryEvalBegin, OnFeedbackBegin, this)
            //    .Register<string>(GameEvents.ProfileStarting, SetUserCode, this)

            
            // Analytics Events
            // text click
            //Game.Events.Register(GameEvents.TextClicked, LogTextClick, this)
            // display text dialog
                
            // SceneHelper.OnSceneLoaded += LogSceneChanged;

            // CrashHandler.OnCrash += OnCrash;

            // NetworkStats.OnError.Register(OnNetworkError);

            m_Log = new OGDLog(new OGDLogConsts() {
                AppId = m_AppId,
                AppVersion = m_AppVersion,
                ClientLogVersion = 1
            });
            m_Log.UseFirebase(m_Firebase);

            #if DEVELOPMENT
                m_Debug = true;
            #endif // DEVELOPMENT

            m_Log.SetDebug(m_Debug);

            m_Platform = GamePlatform.VR;
#if UNITY_WEBGL
            m_Platform = GamePlatform.WEB;
#elif UNITY_ANDROID
            m_Platform = GamePlatform.VR;
#endif

            LogStartGame();
        }

        private void SetUserCode(string userCode)
        {
            Debug.Log("[Analytics] Setting user code: " + userCode);
            m_Log.Initialize(new OGDLogConsts() {
                AppId = m_AppId,
                AppVersion = m_AppVersion,
                ClientLogVersion = 1,
                AppBranch = BuildInfo.Branch()
            });
            m_Log.SetUserId(userCode);
        }

        protected void Shutdown()
        {
            // Game.Events?.DeregisterAll(this);
        }
        #endregion // IService

        #region Log Events

        private void LogStartGame()
        {
            Debug.Log("[Analytics] event: start_game" + "\n" + "Platform: " + m_Platform);

            using (var e = m_Log.NewEvent("start_game"))
            {
                e.Param("mode", m_Platform.ToString());
            }
        }

        #endregion // Log Events

        #region Other Events



        #endregion // Other Events

        #region Helpers



        #endregion // Helpers
    }
}
