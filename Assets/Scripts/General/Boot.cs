using System.Collections;
using System.Collections.Generic;
using ThermoVR.UI;
using UnityEditor;
using UnityEngine;
using UnityEngine.SceneManagement;

namespace ThermoVR {
    public class Boot : MonoBehaviour
    {
        [SerializeField] private string m_firstScene;
        [SerializeField] private UILoading m_loadingCanvas;

        // Start is called before the first frame update
        void Start()
        {
            StartCoroutine(LoadAsync());
        }

        private IEnumerator LoadAsync()
        {
            AsyncOperation asyncLoad = SceneManager.LoadSceneAsync(m_firstScene);
            asyncLoad.allowSceneActivation = false;

            while (asyncLoad.progress < 0.9f)
            {
                yield return null;
            }

            while (!m_loadingCanvas.MinLoadTimeCompleted())
            {
                yield return null;
            }

            asyncLoad.allowSceneActivation = true;
        }
    }
}
