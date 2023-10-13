using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class GameDB : MonoBehaviour
{
    public static GameDB Instance;

    public Material AvailableMat;
    public Material SnapMat;

    public Sprite MCFill;
    public Sprite SocketEmpty;
    public Sprite Correct, Incorrect, Missed;

    public Color TabSelectedColor, TabDefaultColor;

    private void Awake() {
        Instance = this;
    }
}
