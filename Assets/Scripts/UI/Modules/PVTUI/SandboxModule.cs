using ThermoVR.UI;

public class SandboxModule : UIModule
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
