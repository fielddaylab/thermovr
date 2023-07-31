using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class GameDB : MonoBehaviour
{
    public static GameDB Instance;

    public Material AvailableMat;
    public Material SnapMat;

    public Sprite MCFill, MCCorrect, MCIncorrect;

    private void Awake() {
        Instance = this;
    }
}
