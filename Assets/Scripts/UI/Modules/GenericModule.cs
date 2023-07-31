using ThermoVR.UI;

/// <summary>
/// Useful for screens that only open and close
/// </summary>
public class GenericModule : UIModule
{
    #region IUIModule

    public override void Open() {
        base.Open();
    }

    public override void Close() {
        base.Close();
    }

    #endregion // IUIModule
}
