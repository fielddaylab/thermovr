using System;
using ThermoVR;
using ThermoVR.UI;
using UnityEngine;

public class QuizLoadedModule : UIModule
{
    [SerializeField] private ThermoButton m_beginButton;

    #region IUIModule

    public override void Open() {
        base.Open();

        m_beginButton.OnButtonPressed += HandleBeginPressed;
    }

    public override void Close() {
        base.Close();

        m_beginButton.OnButtonPressed -= HandleBeginPressed;
    }

    #endregion // IUIModule

    #region Handlers

    private void HandleBeginPressed(object sender, EventArgs args) {
        GameMgr.Events.Dispatch(GameEvents.BeginLab);
    }

    #endregion // Handlers
}
