using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Events;
using UnityEngine.EventSystems;
using UnityEngine.UI;

namespace ThermoVR.UI
{
    public class ClickToggle : Toggle
    {
        public UnityEvent PointerClicked;

        public override void OnPointerClick(PointerEventData eventData)
        {
            // base.OnPointerClick(eventData);
            PointerClicked?.Invoke();
        }
    }
}