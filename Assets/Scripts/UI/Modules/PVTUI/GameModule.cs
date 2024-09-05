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

    private static float FAILURES_BEFORE_CONTINUE = 1;

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
    [SerializeField] private Image m_loadingOverlayImg;
    [SerializeField] private TMP_Text m_scoreNumText;
    [SerializeField] private GameObject m_graph;
    [SerializeField] private PlacementDotInteractions m_pdInteractions;
    [SerializeField] private ThermoButton m_resetScoreButton;
    [SerializeField] private GameObject m_targetZone;
    [SerializeField] private SphereCollider m_targetCollider;
    [SerializeField] private Collider m_pvtOverlayCollider;

    [SerializeField] private GameObject m_xExtentAnchor;
    [SerializeField] private GameObject m_yExtentAnchor;
    [SerializeField] private GameObject m_zExtentAnchor;


    #endregion //  Inspector

    private int m_currScore = 0;
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

        m_generateTargetRoutine.Stop();

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

        m_loadingOverlayImg.gameObject.SetActive(false);

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

        EventMgr.Events.Dispatch(GameEvents.GameModeScoreUpdated, newScore);
    }

    private IEnumerator GenerateTargetRoutine()
    {
        m_stateIndicatorImg.gameObject.SetActive(false);
        m_loadingOverlayImg.gameObject.SetActive(true);
        m_targetZone.gameObject.SetActive(false);
        EventMgr.Events.Dispatch(GameEvents.GameModeBeginGenerateTarget);

        var defaultRange = 2f;

        var lowerX = m_graph.transform.position.x;
        var lowerY = m_graph.transform.position.y;
        var lowerZ = m_graph.transform.position.z;

        // Calculate the distance in world (not local units)
        // from the bottom-left corner of the graph (by p and v) to the outer bounds.
        // Will change if graph changes in size
        Vector3 extentDist = new Vector3(
            m_xExtentAnchor.transform.position.x - m_graph.transform.position.x,
            m_yExtentAnchor.transform.position.y - m_graph.transform.position.y,
            m_zExtentAnchor.transform.position.z - m_graph.transform.position.z
            );
            // 0.468314f * 0.95f, 0.468314f /*0.1243223f*/);

        var upperX = m_graph.transform.position.x + extentDist.x;
        var upperY = m_graph.transform.position.y + extentDist.y;
        var upperZ = m_graph.transform.position.z + extentDist.z;

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
            //float xPos = m_graph.transform.position.x + 0.5f; // UnityEngine.Random.Range(lowerX, upperX);
            //float yPos = m_graph.transform.position.y + 0.5f;  // UnityEngine.Random.Range(lowerY, upperY);
            //float zPos = m_graph.transform.position.z + 0.5f; // UnityEngine.Random.Range(lowerZ, upperZ);

            float xPos = UnityEngine.Random.Range(lowerX, upperX);
            float yPos = UnityEngine.Random.Range(lowerY, upperY);
            float zPos = UnityEngine.Random.Range(lowerZ, upperZ);

            Vector3 interactPos = new Vector3(xPos, yPos, zPos);

            Vector3 localspace = m_graph.transform.InverseTransformPoint(interactPos);
            Vector3 correctedspace = new Vector3(localspace.z, localspace.y, localspace.x); // * 4.0f; //rotate 90, mul by 4 (inverse transform of gmodel)

            var debugOffGraph = new Vector3(0.51972f, 2314042f, 493.0317f);

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
            var newPos = localguess;

            if (newPos.x < 0.1f
                || newPos.x > 0.9f
                || newPos.y < 0.1f
                || newPos.y > 0.9f
                || newPos.z < 0.1f
                || newPos.z > 0.9f
                )
            {
                isValid = false;
            }

            /*
            // ensure pvt and ball overlap
            if (!Physics.ComputePenetration(m_targetCollider, m_targetZone.transform.position, m_targetZone.transform.rotation,
                m_pvtOverlayCollider, m_pvtOverlayCollider.gameObject.transform.position, m_pvtOverlayCollider.gameObject.transform.rotation, out Vector3 dir, out float dist))
            {
                isValid = false;
            }
            */

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

        m_stateIndicatorImg.gameObject.SetActive(true);
        m_loadingOverlayImg.gameObject.SetActive(false);
        m_targetZone.gameObject.SetActive(true);
        m_generatingNewTarget = false;

        EventMgr.Events.Dispatch(GameEvents.GameModeCompleteGenerateTarget,
            new Tuple<float, float, float>(pTarget.TargetVal, vTarget.TargetVal, tTarget.TargetVal));
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
        EventMgr.Events.Dispatch(GameEvents.ClickGameStart);
        MoveToGameWindow();
    }

    private void HandleHomeButtonPressed(object sender, EventArgs args)
    {
        EventMgr.Events.Dispatch(GameEvents.ClickGameStop);

        MoveToHomeWindow();
    }

    private void HandleResetScorePressed(object sender, EventArgs args)
    {
        EventMgr.Events.Dispatch(GameEvents.ClickGameScoreReset);

        SetScore(0);
    }

    #endregion // Handlers
}
