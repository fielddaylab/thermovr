using ThermoVR.UI;

/// <summary>
/// Useful for screens that only open and close
/// </summary>
public class GenericModule : UIModule
{
    #region IUIModule

    public override void Open() {
        this.gameObject.SetActive(true);
    }

    public override void Close() {
        this.gameObject.SetActive(false);
    }

    #endregion // IUIModule
}
