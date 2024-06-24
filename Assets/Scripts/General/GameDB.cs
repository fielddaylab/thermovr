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

    [Space(5)]
    [Header("Tool Activation")]

    public Material KnobActive;
    public Material KnobNudge;
    public Material KnobGrabbed;
    public Material KnobInactive;
    public Material KnobLockedActive;
    public Material KnobLockedInactive;

    [Space(5)]
    [Header("Audio")]
    public AudioClip TaskCompleteClip;
    public AudioClip SimResetClip;

    private void Awake() {
        Instance = this;
    }
}
