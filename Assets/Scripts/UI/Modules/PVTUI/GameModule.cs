using BeauRoutine;
using BeauUtil;
using System;
using System.Collections;
using System.Collections.Generic;
using ThermoVR;
using ThermoVR.Lab;
using ThermoVR.State;
using ThermoVR.UI;
using TMPro;
using UnityEngine;

public class GameModule : UIModule
{
    #region Consts

    private static float P_RANGE = 100000;
    private static float V_RANGE = 0.01f;
    private static float T_RANGE = 5;

    private static float P_MARGIN = 10000;
    private static float V_MARGIN = 0f;
    private static float T_MARGIN = 5;

    private static float FAILURES_BEFORE_CANCEL = 5;

    #endregion // Consts

    #region Inspector

    [Header("Home")]
    [SerializeField] private CanvasGroup m_homeGroup;
    [SerializeField] private ThermoButton m_beginButton;

    [Header("Game")]
    [SerializeField] private CanvasGroup m_gameGroup;
    [SerializeField] private ThermoButton m_homeButton;
    [SerializeField] private ReachStateHub m_reachStateHub;
    [SerializeField] private TMP_Text m_scoreNumText;
    [SerializeField] private GameObject m_graph;
    [SerializeField] private PlacementDotInteractions m_pdInteractions;

    #endregion //  Inspector

    private int m_currScore;
    private Routine m_reachStateRoutine;
    private int m_failedTargetCount;

    private ReachStateDefinition m_currTargetDef;

    #region Unity Callbacks

    private void Update()
    {
        if (m_gameGroup.alpha == 1 && m_reachStateHub.IsCorrect())
        {
            OnStateReached();

            // TODO: dispatch logging event
        }
    }

    #endregion // Unity Callbacks

    #region Helpers

    private void MoveToHomeWindow()
    {
        if (m_gameGroup.alpha == 1)
        {
            // came from game mode
            EventMgr.Events.Dispatch(GameEvents.GameModeExited);
        }

        SetHomePanelVisible(true);
        SetGamePanelVisible(false);

        // disable reach state hub
        m_reachStateHub.gameObject.SetActive(false);
    }

    private void MoveToGameWindow()
    {
        // entering game mode
        EventMgr.Events.Dispatch(GameEvents.GameModeStarted);

        SetHomePanelVisible(false);
        SetGamePanelVisible(true);

        // generate a new target on open
        m_failedTargetCount = 0;
        GenerateTarget();

        // reset score to 0 on open
        SetScore(0);

        // enable reach state hub
        m_reachStateHub.gameObject.SetActive(true);
    }

    private void AddListeners()
    {
        m_beginButton.OnButtonPressed += HandleBeginButtonPressed;
        m_homeButton.OnButtonPressed += HandleHomeButtonPressed;
    }

    private void RemoveListeners()
    {
        m_beginButton.OnButtonPressed += HandleBeginButtonPressed;
        m_homeButton.OnButtonPressed += HandleHomeButtonPressed;
    }

    private void SetHomePanelVisible(bool isVisible)
    {
        m_homeGroup.alpha = isVisible ? 1 : 0;
        m_homeGroup.interactable = isVisible;
        m_homeGroup.blocksRaycasts = isVisible;
    }

    private void SetGamePanelVisible(bool isVisible)
    {
        m_gameGroup.alpha = isVisible ? 1 : 0;
        m_gameGroup.interactable = isVisible;
        m_gameGroup.blocksRaycasts = isVisible;
    }

    private void SetScore(int newScore)
    {
        m_currScore = newScore;

        m_scoreNumText.SetText(newScore.ToStringLookup());
    }

    private void GenerateTarget()
    {
        // randomly pick position
        double v = UnityEngine.Random.Range(2.0f, 3.0f);
        //double v = UnityEngine.Random.Range((float)ThermoMath.v_min + V_RANGE, (float)ThermoMath.v_max - V_RANGE); ;

        // graph bounds are 0.25f
        // x = t, y = p, z = v
        float xPos = UnityEngine.Random.Range(m_graph.transform.position.x + 0.5f, m_graph.transform.position.x - 0.1f);
        float yPos = UnityEngine.Random.Range(m_graph.transform.position.y + 0.5f, m_graph.transform.position.y + 0.1f);
        float zPos = UnityEngine.Random.Range(m_graph.transform.position.z + 0.5f, m_graph.transform.position.z + 0.1f);

        Vector3 interactPos = new Vector3(xPos, yPos, zPos);

        Vector3 localspace = m_graph.transform.InverseTransformPoint(interactPos);
        Vector3 correctedspace = new Vector3(localspace.z, localspace.y, localspace.x); // * 4.0f; //rotate 90, mul by 4 (inverse transform of gmodel)

        //Vector3 thermoguess = thermo.guessPlot(ThermoMath.t_neutral, correctedspace.y, correctedspace.x);
        Vector3 thermoguess = ThermoPresent.Instance.guessMeshPlot(correctedspace.x, correctedspace.y, correctedspace.z);
        Vector3 localguess = ThermoPresent.Instance.plot(thermoguess.y, thermoguess.x, thermoguess.z); //note swizzle!

        bool isValid = true;

        if (MathUtility.floatNumeric(localguess.x) && MathUtility.floatNumeric(localguess.y) && MathUtility.floatNumeric(localguess.z))
        {
            // interactPos = thermoguess;
            isValid = true;
        }
        else
        {
            isValid = false;
        }

        // x = v, y = p, z = t
        var pvt = thermoguess; // ThermoPresent.Instance.invplot(interactPos.y, interactPos.z, interactPos.x);

        if (pvt.y < ThermoMath.p_min || pvt.x < ThermoMath.v_min || pvt.z < ThermoMath.t_min
            || pvt.y > ThermoMath.p_max || pvt.x > ThermoMath.v_max || pvt.z > ThermoMath.t_max)
        {
            isValid = false;
        }

        if (!isValid)
        {
            if (m_failedTargetCount > FAILURES_BEFORE_CANCEL)
            {
                // TODO: handle failed target generation visually
                return;
            }

            m_failedTargetCount++;
            GenerateTarget();
            return;
        }

        // Create a new reach state def
        SimStateTarget pTarget = new SimStateTarget();
        pTarget.TargetID = VarID.Pressure;
        pTarget.TargetVal = (float)pvt.y / 1000f;
        pTarget.TargetRange = P_RANGE / 1000f;

        SimStateTarget vTarget = new SimStateTarget();
        vTarget.TargetID = VarID.Volume;
        vTarget.TargetVal = (float)pvt.x;
        vTarget.TargetRange = V_RANGE;

        SimStateTarget tTarget = new SimStateTarget();
        tTarget.TargetID = VarID.Temperature;
        tTarget.TargetVal = (float)pvt.z;
        tTarget.TargetRange = T_RANGE;

        List<SimStateTarget> target = new List<SimStateTarget>() {
            pTarget,
            vTarget,
            tTarget
        };

        m_currTargetDef = new ReachStateDefinition(
            String.Empty,
            null,
            target
            );

        m_reachStateHub.SetDefinition(m_currTargetDef);
    }

    private void OnStateReached()
    {
        m_reachStateRoutine.Replace(StateReachedRoutine());
    }

    #endregion // Helpers

    #region Routines

    private IEnumerator StateReachedRoutine() {
        // play animation

        // play sound

        // increment score
        SetScore(m_currScore + 1);

        // pick a new target
        GenerateTarget();

        yield return null;
    }


    #endregion // Routines

    #region IUIModule

    public override void Init()
    {
        base.Init();
    }

    public override void Open() {
        this.gameObject.SetActive(true);

        AddListeners();
        MoveToHomeWindow();
    }

    public override void Close() {
        MoveToHomeWindow();

        this.gameObject.SetActive(false);

        RemoveListeners();
        m_reachStateRoutine.Stop();
    }

    #endregion // IUIModule

    #region Handlers

    private void HandleBeginButtonPressed(object sender, EventArgs args)
    {
        MoveToGameWindow();
    }

    private void HandleHomeButtonPressed(object sender, EventArgs args)
    {
        MoveToHomeWindow();
    }

    #endregion // Handlers
}
