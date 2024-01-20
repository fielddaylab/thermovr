using System;
using UnityEngine;

public class ListenerFinder : MonoBehaviour
{

    AudioSource[] sources;
    void Start()
    {

        //Get every single audio sources in the scene.
        sources = GameObject.FindObjectsOfType(typeof(AudioSource)) as AudioSource[];
    }

    void Update()
    {

        // When a key is pressed list all the gameobjects that are playing an audio


            foreach (AudioSource audioSource in sources)
            {
                if (audioSource.isPlaying)
                {
                    Debug.Log(audioSource.name + " is playing " + audioSource.clip.name);
                    Debug.Log("---------------------------"); //to avoid confusion next time
                    Debug.Break(); //pause the editor
              }
            }
        // }
    }
}