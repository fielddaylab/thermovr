using BeauRoutine;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using BeauUtil;
using System;
using BeauUtil.Debugger;


#if UNITY_EDITOR
using UnityEditor;
#endif // UNITY_EDITOR

namespace ThermoVR.UI
{
    public struct CreditsBlockData
    {
        public string Header;
        public string Names;
    }

    public class UICredits : MonoBehaviour
    {
        private static string CREDITS_DELIM = "::";

        private static float TRANSITION_TIME = 0.5f;

        #region Inspector

        [SerializeField] private Button m_ReturnButton;
        [SerializeField] private RectTransform m_Rect;
        // [SerializeField] private Canvas m_Canvas;

        [Space(5)]
        [Header("Text")]
        [SerializeField] private TextAsset m_CreditsTextAsset;
        [SerializeField] private RectTransform m_TextContainer;
        [SerializeField] private LayoutGroup m_TextLayout;

        [SerializeField] private GameObject m_CreditsBlockPrefab;

        [Space(5)]
        [Header("Settings")]
        [SerializeField] private float m_ScrollTime;
        [SerializeField] private float m_ScrollTimeBuffer;

        [Space(5)]
        [Header("Spacing")]
        [SerializeField] private Vector2 m_StartPos;
        [SerializeField] private float m_SpacingPerChunk = 75;
        [SerializeField] private float m_SpacingPerGroupLine = 15;
        [SerializeField] private float m_SpacingPerNamesLine = 5;

        private bool m_Scrolling;
        private float m_ScrollSpeed;
        private float m_ScrollTimer;

        private bool m_Initialized;

        private Routine m_ConceptRoutine;

        private List<CreditsBlockData> m_CreditsBlockDatas;
        private List<CreditsBlock> m_CreditsBlocks;

        #endregion // Inspector

        private Routine m_TransitionRoutine;
        private bool m_Transitioning;

        #region Unity Callbacks

        private void OnEnable()
        {
            if (m_Initialized)
            {
                SetupScroll();
            }
        }

        private void Update()
        {
            if (m_Scrolling && m_ScrollTimer > -m_ScrollTimeBuffer)
            {
                m_TextContainer.transform.position = new Vector3(
                    m_TextContainer.transform.position.x,
                    m_TextContainer.transform.position.y + m_ScrollSpeed * Time.deltaTime,
                    m_TextContainer.transform.position.z
                    );

                m_ScrollTimer -= Time.deltaTime;
            }
        }

        #endregion // Unity Callbacks

        public void Init(bool isLaunch)
        {
            m_ReturnButton.onClick.AddListener(HandleReturnClicked);

            ParseCredits();
            CreditsBlocksToText();

            EventMgr.Events.Register(GameEvents.TitleCreditsOpened, OpenPanelImmediate);

            m_Initialized = true;

            SetupScroll();
            /*
            if (isLaunch)
            {
                m_Canvas.renderMode = RenderMode.ScreenSpaceOverlay;
            }
            else
            {
                m_Canvas.renderMode = RenderMode.ScreenSpaceCamera;
            }
            */
        }

        #region Button Handlers

        private void HandleReturnClicked()
        {
            EventMgr.Events.Dispatch(GameEvents.TitleCreditsClosed);
            if (m_Transitioning) { return; }
            m_TransitionRoutine.Replace(this, ExitRoutine()).ExecuteWhileDisabled();
        }

        #endregion // Button Handlers

        #region External

        public void OpenPanel()
        {
            if (m_Transitioning) { return; }
            m_TransitionRoutine.Replace(this, EnterRoutine(TRANSITION_TIME)).ExecuteWhileDisabled();
        }

        public void OpenPanelImmediate()
        {
            if (m_Transitioning) { return; }
            m_TransitionRoutine.Replace(this, EnterRoutine(0)).ExecuteWhileDisabled();
        }

        #endregion // External

        #region Routines

        private IEnumerator EnterRoutine(float enterTime)
        {
            m_Transitioning = true;
            m_ScrollTimer = m_ScrollTime;
            this.gameObject.SetActive(true);
            // yield return m_Rect.AnchorPosTo(0, TRANSITION_TIME, Axis.X).Ease(Curve.CubeOut);
            yield return null;
            m_Scrolling = true;
            m_Transitioning = false;
        }

        private IEnumerator ExitRoutine()
        {
            m_Transitioning = true;
            this.gameObject.SetActive(false);
            // m_Rect.SetAnchorsPreservePosition(new Vector2(1, 0), new Vector2(2, 1));
            // yield return m_Rect.AnchorPosTo(0, TRANSITION_TIME, Axis.X).Ease(Curve.CubeOut);
            yield return null;

            ResetLayout();

            m_Transitioning = false;
            // ZavalaGame.Events.Dispatch(GameEvents.CreditsExited);
        }

        #endregion // Routines

        private void ResetLayout()
        {
            // Text position
            m_TextContainer.transform.localPosition = new Vector3(
                m_TextContainer.transform.localPosition.x,
                -700 - m_TextContainer.rect.size.y,
                m_TextContainer.transform.localPosition.z
                );

            m_Scrolling = false;
        }


        private void SetupScroll()
        {
            var viewHeight = 1400;

            // Calculate scroll speed based on scroll time and distance to travel
            Assert.True(m_ScrollTime > 0);
            if (m_CreditsBlocks.Count > 0)
            {
                m_ScrollSpeed = (m_CreditsBlocks[0].transform.position.y - m_CreditsBlocks[m_CreditsBlocks.Count - 1].transform.position.y + m_TextContainer.rect.size.y + viewHeight) / (m_ScrollTime);
            }

            ResetLayout();
        }

        #region Credits Parsing

        private void ParseCredits()
        {
            List<string> blocksToParse = TextIO.TextAssetToList(m_CreditsTextAsset, CREDITS_DELIM);
            m_CreditsBlockDatas = new List<CreditsBlockData>();

            for (int i = 0; i < blocksToParse.Count; i++)
            {
                // skip first space
                if (i == 0) { continue; }

                m_CreditsBlockDatas.Add(ParseCreditsBlock(blocksToParse[i]));
            }
        }

        private CreditsBlockData ParseCreditsBlock(string blocktoParse)
        {
            CreditsBlockData newBlock = new CreditsBlockData();

            // Trim leading delim
            int startIdx = blocktoParse.IndexOf(CREDITS_DELIM) + CREDITS_DELIM.Length;

            // Split at first newline
            int newlineIdx;

            int nIndex = blocktoParse.IndexOf("\n");
            int rIndex = blocktoParse.IndexOf("\r");

            newlineIdx = Math.Max(nIndex, rIndex);

            newBlock.Header = blocktoParse.Substring(startIdx, newlineIdx - startIdx).Trim();
            newBlock.Names = blocktoParse.Substring(newlineIdx + 1, blocktoParse.Length - newlineIdx - 1);

            return newBlock;
        }

        private void CreditsBlocksToText()
        {
            m_CreditsBlocks = new List<CreditsBlock>();

            for (int i = 0; i < m_CreditsBlockDatas.Count; i++)
            {
                var newTextBlock = Instantiate(m_CreditsBlockPrefab, m_TextContainer.transform).GetComponent<CreditsBlock>();
                newTextBlock.Header.SetText(m_CreditsBlockDatas[i].Header);
                newTextBlock.Header.ForceMeshUpdate();
                newTextBlock.Names.SetText(m_CreditsBlockDatas[i].Names);
                newTextBlock.Names.ForceMeshUpdate();

                m_CreditsBlocks.Add(newTextBlock);
            }

            // m_TextLayout.ForceRebuild();

            m_TextLayout.enabled = false;

            ApplySpacing();
        }

        public void ApplySpacing()
        {
            float currX = m_CreditsBlocks[0].transform.position.x;
            float currY = m_CreditsBlocks[0].transform.position.y;

            for (int i = 0; i < m_CreditsBlocks.Count; i++)
            {
                m_CreditsBlocks[i].transform.position = new Vector2(currX, currY);

                int numGroupLines = m_CreditsBlocks[i].Header.textInfo.lineCount;
                int numNamesLines = m_CreditsBlocks[i].Names.textInfo.lineCount;

                currY -= (m_SpacingPerChunk + (m_CreditsBlocks[i].Header.fontSize + m_SpacingPerGroupLine) * numGroupLines + (m_CreditsBlocks[i].Names.fontSize + m_SpacingPerNamesLine) * numNamesLines);

                m_CreditsBlocks[i].Layout.ForceRebuild();
            }
        }

        #endregion // Credits Parsing
    }
}