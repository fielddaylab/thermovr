using BeauRoutine;
using System;
using System.Collections;
using System.Collections.Generic;
using TMPro;
using UnityEngine;
using UnityEngine.UI;

namespace ThermoVR
{
    public class ButtonPressMovement : MonoBehaviour
    {
        [SerializeField] private Pressable m_target;
        [SerializeField] private float m_indentAmt;
        [SerializeField] private Axis m_axis;

        [SerializeField] private bool m_explicitCalls;

        [SerializeField] private bool m_hasTimer;
        [SerializeField] private float m_resetTime;
        private float m_resetTimer;

        [Header("Materials")]
        [SerializeField] private bool m_changesMat;
        [SerializeField] private MeshRenderer m_meshRenderer;
        [SerializeField] private Material m_changeToMat;
        [SerializeField] private int m_matIndex;
        private Material m_defaultMat;

        [SerializeField] private TMP_Text m_textRenderer;
        [SerializeField] private SpriteRenderer m_imageRenderer;
        [SerializeField] private Color m_changeToColor;

        private Color m_defaultColor;

        private bool m_indented;

        private Vector3 m_defaultLocalPos;
        private Vector3 m_indentedLocalPos;

        private void OnEnable()
        {
            m_indented = false;
            m_target.OnPress += HandlePress;

            m_defaultLocalPos = transform.localPosition;
            Vector3 newPos = this.transform.localPosition;
            if (m_axis == Axis.X)
            {
                newPos.x -= m_indentAmt;
            }
            else if (m_axis == Axis.Y)
            {
                newPos.y -= m_indentAmt;
            }
            else if (m_axis == Axis.Z)
            {
                newPos.z -= m_indentAmt;
            }

            m_indentedLocalPos = newPos;

            if (m_changesMat)
            {
                m_defaultMat = m_meshRenderer.sharedMaterials[m_matIndex];
                if (m_imageRenderer)
                {
                    m_defaultColor = m_imageRenderer.color;
                }
                if (m_textRenderer)
                {
                    m_defaultColor = m_textRenderer.color;
                }
            }
        }

        private void OnDisable()
        {
            m_target.OnPress -= HandlePress;
        }

        private void Update()
        {
            if (!m_hasTimer) { return; }

            if (m_resetTimer > 0)
            {
                m_resetTimer -= Time.deltaTime;

                if (m_resetTimer <= 0)
                {
                    ResetPosition();
                }
            }
        }

        public void ResetPosition()
        {
            this.transform.localPosition = m_defaultLocalPos;

            m_indented = false;

            if (m_changesMat)
            {
                var mats = m_meshRenderer.sharedMaterials;
                mats[m_matIndex] = m_defaultMat;
                m_meshRenderer.sharedMaterials = mats;

                // set image / text color
                if (m_imageRenderer)
                {
                    m_imageRenderer.color = m_defaultColor;
                }
                if (m_textRenderer)
                {
                    m_textRenderer.color = m_defaultColor;
                }
            }
        }

        public void Indent()
        {
            this.transform.localPosition = m_indentedLocalPos;

            if (m_changesMat)
            {
                var mats = m_meshRenderer.sharedMaterials;
                mats[m_matIndex] = m_changeToMat;
                m_meshRenderer.sharedMaterials = mats;

                // set image / text color
                if (m_imageRenderer)
                {
                    m_imageRenderer.color = m_changeToColor;
                }
                if (m_textRenderer)
                {
                    m_textRenderer.color = m_changeToColor;
                }
            }
        }

        private void HandlePress(object sender, EventArgs args)
        {
            if (m_indented || m_explicitCalls) { return; }

            // move button in/out
            Indent();

            if (m_hasTimer)
            {
                m_resetTimer = m_resetTime;
            }

            m_indented = true;
        }
    }
}