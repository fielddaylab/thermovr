using System.Collections;
using System.Collections.Generic;
using System.Threading;
using ThermoVR.UI;
using UnityEditor;
using UnityEngine;
using UnityEngine.SceneManagement;

namespace ThermoVR {
    public class Boot : MonoBehaviour
    {
        [SerializeField] private string m_firstScene;
        [SerializeField] private UILoading m_loadingCanvas;

        private AsyncOperation m_asyncLoad;

        // Start is called before the first frame update
        void Start()
        {
            EventMgr.Events.Register(GameEvents.StartGameClicked, HandleStartGameClicked);

            StartCoroutine(LoadAsync());
        }

        private IEnumerator LoadAsync()
        {
            m_asyncLoad = SceneManager.LoadSceneAsync(m_firstScene);
            m_asyncLoad.allowSceneActivation = false;

            while (m_asyncLoad.progress < 0.9f)
            {
                yield return null;
            }

            while (!m_loadingCanvas.MinLoadTimeCompleted())
            {
                yield return null;
            }

            if (ModeMgr.Instance.IsDesktop)
            {
                // Wait for button
                EventMgr.Events.Dispatch(GameEvents.AsyncLoadComplete);
            }
            else
            {
                // Load immediately
                m_asyncLoad.allowSceneActivation = true;
            }
        }

        #region Handlers

        private void HandleStartGameClicked()
        {
            m_asyncLoad.allowSceneActivation = true;

            if (m_loadingCanvas.GetFullscreen())
            {
                // trigger fullscreen
                Screen.fullScreen = true;
                #if UNITY_WEBGL && !UNITY_EDITOR
                // NativeFullscreen_SetFullscreen(true);
                #endif // UNITY_WEBGL && !UNITY_EDITOR
            }

            if (PersistentState.Instance.Bools.ContainsKey(PersistentVars.MusicOnStart))
            {
                PersistentState.Instance.Bools.Add(PersistentVars.MusicOnStart, m_loadingCanvas.GetMusic());
            }
            else
            {
                PersistentState.Instance.Bools[PersistentVars.MusicOnStart] = m_loadingCanvas.GetMusic();
            }
        }

        #endregion // Handlers
    }
}
