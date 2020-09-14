using UnityEngine;
using System.Collections;
using UnityEngine.VR;
using UnityEngine.SceneManagement;

public class Splash : MonoBehaviour
{

  public float minimumTimeToShowLogo = 5f;
  public string levelToLoad = "MainScene";

  public string introName = "SplashIn";
  public string outroName = "SplashOut";

  IEnumerator Start()
  {
    float minimumTimeEnd = Time.realtimeSinceStartup + minimumTimeToShowLogo;

    // intro

    // background load the new scene (but don't activate it yet)

    AsyncOperation o = SceneManager.LoadSceneAsync(levelToLoad);
    o.allowSceneActivation = true;
    while(o.isDone)
    {
      yield return new WaitForEndOfFrame();
    }

    // delay until minimum time is reached

    if (Time.realtimeSinceStartup < minimumTimeEnd)
    {
      yield return new WaitForSeconds(minimumTimeEnd - Time.realtimeSinceStartup);
    }
  }
}