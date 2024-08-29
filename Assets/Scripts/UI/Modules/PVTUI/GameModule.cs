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
using UnityEngine.UI;

public class GameModule : UIModule
{
    #region Consts

    private static float P_RANGE = 0;
    private static float V_RANGE = 0;
    private static float T_RANGE = 0;

    private static float P_MARGIN = 10000;
    private static float V_MARGIN = 0f;
    private static float T_MARGIN = 5;

    private static float FAILURES_BEFORE_CONTINUE = 5;

    #endregion // Consts

    #region Inspector

    [Header("Home")]
    [SerializeField] private CanvasGroup m_homeGroup;
    [SerializeField] private ThermoButton m_beginButton;

    [Header("Game")]
    [SerializeField] private CanvasGroup m_gameGroup;
    [SerializeField] private ThermoButton m_homeButton;
    [SerializeField] private ReachStateHub m_reachStateHub;
    [SerializeField] private Image m_stateIndicatorImg;
    [SerializeField] private TMP_Text m_scoreNumText;
    [SerializeField] private GameObject m_graph;
    [SerializeField] private PlacementDotInteractions m_pdInteractions;
    [SerializeField] private ThermoButton m_resetScoreButton;
    [SerializeField] private GameObject m_targetZone;

    #endregion //  Inspector

    private int m_currScore;
    private Routine m_reachStateRoutine;
    private Routine m_generateTargetRoutine;
    private int m_failedTargetCount;

    private ReachStateDefinition m_currTargetDef;

    private float m_debugTimer = 0;
    private bool m_generatingNewTarget = false;

    #region Unity Callbacks

    private void Update()
    {
        if (m_gameGroup.alpha == 1 && m_reachStateHub.IsCorrect() && !m_generatingNewTarget)
        {
            OnStateReached();

            // TODO: dispatch logging event
        }

        m_debugTimer -= Time.deltaTime;
        if (m_debugTimer <= 0)
        { 
            m_debugTimer = 0.5f;
            // GenerateTarget();
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

        // allow dragging
        World.Instance.ModMgr.EnableGraphBallInteractions();

        // hide target zone
        m_targetZone.gameObject.SetActive(false);

        // disable reach state hub
        m_reachStateHub.gameObject.SetActive(false);
    }

    private void MoveToGameWindow()
    {
        // entering game mode
        EventMgr.Events.Dispatch(GameEvents.GameModeStarted);

        SetHomePanelVisible(false);
        SetGamePanelVisible(true);

        // Disallow dragging
        World.Instance.ModMgr.DisableGraphBallInteractions();

        // generate a new target on open
        m_failedTargetCount = 0;
        m_generateTargetRoutine.Replace(GenerateTargetRoutine());

        // enable reach state hub
        m_reachStateHub.gameObject.SetActive(true);
    }

    private void AddListeners()
    {
        m_beginButton.OnButtonPressed += HandleBeginButtonPressed;
        m_homeButton.OnButtonPressed += HandleHomeButtonPressed;
        m_resetScoreButton.OnButtonPressed += HandleResetScorePressed;
    }

    private void RemoveListeners()
    {
        m_beginButton.OnButtonPressed -= HandleBeginButtonPressed;
        m_homeButton.OnButtonPressed -= HandleHomeButtonPressed;
        m_resetScoreButton.OnButtonPressed -= HandleResetScorePressed;
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

    private IEnumerator GenerateTargetRoutine()
    {
        m_stateIndicatorImg.sprite = GameDB.Instance.Circle;

        bool isValid = false;
        int numTriesThisFrame = 0;
        Vector3 finalPVT = Vector3.zero;
        while (!isValid)
        {
            if (numTriesThisFrame > FAILURES_BEFORE_CONTINUE)
            {
                numTriesThisFrame = 0;
                yield return null;
            }

            // randomly pick position
            // x = t, y = p, z = v
            float xPos = UnityEngine.Random.Range(m_graph.transform.position.x + 0.3f, m_graph.transform.position.x - 0.3f);
            float yPos = UnityEngine.Random.Range(m_graph.transform.position.y + 0.3f, m_graph.transform.position.y - 0.3f);
            float zPos = UnityEngine.Random.Range(m_graph.transform.position.z + 0.3f, m_graph.transform.position.z - 0.3f);

            Vector3 interactPos = new Vector3(xPos, yPos, zPos);

            Vector3 localspace = m_graph.transform.InverseTransformPoint(interactPos);
            Vector3 correctedspace = new Vector3(localspace.z, localspace.y, localspace.x); // * 4.0f; //rotate 90, mul by 4 (inverse transform of gmodel)

            //Vector3 thermoguess = thermo.guessPlot(ThermoMath.t_neutral, correctedspace.y, correctedspace.x);
            Vector3 thermoguess = ThermoPresent.Instance.guessMeshPlot(correctedspace.x, correctedspace.y, correctedspace.z);
            Vector3 localguess = ThermoPresent.Instance.plot(thermoguess.y, thermoguess.x, thermoguess.z); //note swizzle!

            if (MathUtility.floatNumeric(localguess.x) && MathUtility.floatNumeric(localguess.y) && MathUtility.floatNumeric(localguess.z))
            {
                // interactPos = thermoguess;
                isValid = true;
            }

            // x = v, y = p, z = t

            finalPVT = thermoguess; // ThermoPresent.Instance.invplot(interactPos.y, interactPos.z, interactPos.x);

            if (finalPVT.y < ThermoMath.p_min || finalPVT.x < ThermoMath.v_min || finalPVT.z < ThermoMath.t_min
                || finalPVT.y > ThermoMath.p_max || finalPVT.x > ThermoMath.v_max || finalPVT.z > ThermoMath.t_max)
            {
                isValid = false;
            }

            // keep bounds off edge cases
            var newPos = ThermoPresent.Instance.plot(thermoguess.y, thermoguess.x, thermoguess.z);
            if (newPos.x < 0.06f || newPos.y > 0.8f || newPos.z > 0.9f)
            {
                isValid = false;
            }

            numTriesThisFrame++;
        }

        // Create a new reach state def
        SimStateTarget pTarget = new SimStateTarget();
        pTarget.TargetID = VarID.Pressure;
        pTarget.TargetVal = (float)finalPVT.y / 1000f;
        pTarget.TargetRange = P_RANGE / 1000f;

        SimStateTarget vTarget = new SimStateTarget();
        vTarget.TargetID = VarID.Volume;
        vTarget.TargetVal = (float)finalPVT.x;
        vTarget.TargetRange = V_RANGE;

        SimStateTarget tTarget = new SimStateTarget();
        tTarget.TargetID = VarID.Temperature;
        tTarget.TargetVal = (float)finalPVT.z;
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

        m_generatingNewTarget = false;
    }

    private void OnStateReached()
    {
        m_generatingNewTarget = true;
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
        m_generateTargetRoutine.Replace(GenerateTargetRoutine());

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
        m_generateTargetRoutine.Stop();
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

    private void HandleResetScorePressed(object sender, EventArgs args)
    {
        SetScore(0);
    }

    #endregion // Handlers
}
