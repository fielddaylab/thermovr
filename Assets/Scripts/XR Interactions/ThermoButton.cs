using System;
using System.Collections;
using System.Collections.Generic;
using TMPro;
using UnityEngine;
using UnityEngine.EventSystems;
using UnityEngine.UI;

namespace ThermoVR.UI
{
    [RequireComponent(typeof(Pressable))]
    [RequireComponent (typeof(BoxCollider))]
    public class ThermoButton : MonoBehaviour
    {
        public event EventHandler OnButtonPressed; // Wrapper for the Pressable event
        private Pressable Pressable;
        [SerializeField] private Button m_button;
        [SerializeField] private TMP_Text m_text;

        public void OnEnable() {
            // Register physical touch
            Pressable = GetComponent<Pressable>();

            Pressable.OnPress += HandlePress;

            // set the collider to be the size of the button
            RectTransform rect = GetComponent<RectTransform>();
            BoxCollider collider = GetComponent<BoxCollider>();
            Vector3 currSize = collider.size;
            // GetComponent<BoxCollider>().size = new Vector3(rect.sizeDelta.x, rect.sizeDelta.y, currSize.z);

            // Register ray touch
            m_button.onClick.AddListener(HandleClick);
        }

        public void OnDisable() {
            Pressable.OnPress -= HandlePress;
            m_button.onClick.RemoveListener(HandleClick);
        }

        public void SetText(string text) {
            m_text.text = text;
        }

        public void SetColor(Color color) {
            m_button.image.color = color;
        }

        public void ClearListeners()
        {
            OnButtonPressed = null;
        }

        public void SetInteractable(bool interactable)
        {
            m_button.interactable = interactable;
        }

        #region Handlers

        private void HandlePress(object sender, EventArgs args) {
            // TODO: UI click audio?

            // Reroutes Press through the button's onClick event
            ExecuteEvents.Execute(m_button.gameObject, new BaseEventData(EventSystem.current), ExecuteEvents.submitHandler);
        }

        private void HandleClick() {
            OnButtonPressed?.Invoke(this, EventArgs.Empty);
        }

        #endregion // Handlers
    }
}