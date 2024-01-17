using System;
using ThermoVR;
using ThermoVR.UI;
using UnityEngine;

public class QuizLoadedModule : UIModule
{
    [SerializeField] private ThermoButton m_beginButton;

    #region IUIModule

    public override void Open() {
        this.gameObject.SetActive(true);

        m_beginButton.OnButtonPressed += HandleBeginPressed;
    }

    public override void Close() {
        this.gameObject.SetActive(false);

        m_beginButton.OnButtonPressed -= HandleBeginPressed;
    }

    #endregion // IUIModule

    #region Handlers

    private void HandleBeginPressed(object sender, EventArgs args) {
        GameMgr.Events.Dispatch(GameEvents.BeginLab);
    }

    #endregion // Handlers
}
