using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class GameDB : MonoBehaviour
{
    public static GameDB Instance;

    public Material AvailableMat;
    public Material SnapMat;

    public Material InactiveButtonMaterial;
    public Material ActiveButtonMaterial;

    public Sprite MCFill;
    public Sprite SocketEmpty;
    public Sprite Correct, Incorrect, Missed;

    public Sprite ReachStateIncomplete, ReachStateComplete;
    public Sprite LabTaskTabInactive, LabTaskTabActive;
    public Sprite LabTopicTabInactive, LabTopicTabActive;

    public Color TabSelectedColor, TabDefaultColor;
    public Color MCSelectedBG, MCUnselectedBG, MCIncorrectBG;

    private void Awake() {
        Instance = this;
    }
}
