#if UNITY_2019_1_OR_NEWER
#define USE_SRP
#endif // UNITY_2019_1_OR_NEWER

using BeauUtil;
using BeauUtil.Debugger;
using ThermoVR;
using UnityEngine;


public sealed class RenderMgr
{
    private bool m_LastKnownFullscreen;
    private Resolution m_LastKnownResolution;

    private Camera m_PrimaryCamera;

    private RingBuffer<CameraClampToVirtualViewport> m_ClampedViewportCameras = new RingBuffer<CameraClampToVirtualViewport>(2, RingBufferMode.Expand);
    private Rect m_VirtualViewport = new Rect(0, 0, 1, 1);

    private float m_MinAspect;
    private float m_MaxAspect;
    private bool m_HasLetterboxing;

    #region Callbacks

    public readonly CastableEvent<bool> OnFullscreenChanged = new CastableEvent<bool>(2);
    public readonly CastableEvent<Resolution> OnResolutionChanged = new CastableEvent<Resolution>(2);

    #endregion // Callbacks

    #region Events

    internal void Initialize()
    {
        GameMgr.Events.Register(GameEvents.CanvasPreUpdate, OnCanvasPreUpdate);
        GameMgr.Events.Register(GameEvents.ApplicationPreRender, OnApplicationPreRender);
    }

    internal void PollScreenSettings()
    {
        bool fullscreen = ScreenUtility.GetFullscreen();
        if (m_LastKnownFullscreen != fullscreen)
        {
            m_LastKnownFullscreen = fullscreen;
            OnFullscreenChanged.Invoke(fullscreen);
        }

        Resolution resolution = ScreenUtility.GetResolution();
        if (resolution.width != m_LastKnownResolution.width || resolution.height != m_LastKnownResolution.height || resolution.refreshRate != m_LastKnownResolution.refreshRate)
        {
            m_LastKnownResolution = resolution;
            OnResolutionChanged.Invoke(resolution);
        }
    }

    internal void Shutdown()
    {
        GameMgr.Events?.Deregister(GameEvents.CanvasPreUpdate, OnCanvasPreUpdate);
        GameMgr.Events?.Deregister(GameEvents.ApplicationPreRender, OnApplicationPreRender);

        OnResolutionChanged.Clear();
        OnFullscreenChanged.Clear();
    }

    #endregion // Events

    #region World Camera

    public void SetPrimaryCamera(Camera camera)
    {
        if (m_PrimaryCamera != null)
        {
            Log.Warn("[RenderMgr] Primary world camera already set to '{0}' - make sure to deregister it first", m_PrimaryCamera);
        }
        m_PrimaryCamera = camera;
        Log.Msg("[RenderMgr] Assigned primary world camera as '{0}'", camera);
    }

    public void RemovePrimaryCamera(Camera camera)
    {
        if (camera == null || m_PrimaryCamera != camera)
        {
            return;
        }

        m_PrimaryCamera = null;
        Log.Msg("[RenderMgr] Removed primary world camera");
    }

    #endregion // World Camera

    #region Clamped Viewport

    public void EnableAspectClamping(int width, int height)
    {
        m_MinAspect = (float)width / height;
        m_MaxAspect = m_MinAspect;
    }

    public void EnableMinimumAspectClamping(int width, int height)
    {
        m_MinAspect = (float)width / height;
        m_MaxAspect = float.MaxValue;
    }

    public void EnableAspectClamping(Vector2Int min, Vector2Int max)
    {
        m_MinAspect = (float)min.x / min.y;
        m_MaxAspect = (float)max.x / max.y;
    }

    public void DisableAspectClamping()
    {
        m_MinAspect = m_MaxAspect = 0;
        m_VirtualViewport = new Rect(0, 0, 1, 1);
    }

    public Rect VirtualViewport
    {
        get { return m_VirtualViewport; }
    }

    public void AddClampedViewportCamera(CameraClampToVirtualViewport camera)
    {
        Assert.NotNull(camera);
        m_ClampedViewportCameras.PushBack(camera);
    }

    public void RemoveClampedViewportCamera(CameraClampToVirtualViewport camera)
    {
        Assert.NotNull(camera);
        m_ClampedViewportCameras.FastRemove(camera);
    }

    #endregion // Clamped Viewport

    #region Handlers

    private void OnCanvasPreUpdate()
    {
        if (m_MinAspect <= 0 || m_MaxAspect <= 0)
        {
            m_HasLetterboxing = false;
            return;
        }

        m_VirtualViewport = UpdateAspectRatioClamping();

        for (int i = 0; i < m_ClampedViewportCameras.Count; i++)
        {
            ref var c = ref m_ClampedViewportCameras[i];
            Rect r = c.Viewport;
            r.x = m_VirtualViewport.x + r.x * m_VirtualViewport.width;
            r.y = m_VirtualViewport.y + r.y * m_VirtualViewport.height;
            r.width = r.width * m_VirtualViewport.width;
            r.height = r.height * m_VirtualViewport.height;
            c.Camera.rect = r;
        }

        m_HasLetterboxing = true;
    }

    private void OnApplicationPreRender()
    {
        if (m_HasLetterboxing && m_ClampedViewportCameras.Count > 0)
        {
            CameraHelper.RenderLetterboxing(m_VirtualViewport, Color.black);
        }
    }

    private Rect UpdateAspectRatioClamping()
    {
        float currentAspect = (float)m_LastKnownResolution.width / m_LastKnownResolution.height;
        float finalAspect = Mathf.Clamp(currentAspect, m_MinAspect, m_MaxAspect);

        float aspectW = finalAspect;
        float aspectH = 1;

        if (aspectW > currentAspect)
        {
            aspectH = currentAspect / finalAspect;
            aspectW = aspectH * finalAspect;
        }

        float diffX = 1 - (aspectW / currentAspect),
            diffY = 1 - (aspectH / 1);

        Rect r = default;
        r.x = diffX / 2;
        r.y = diffY / 2;
        r.width = 1 - diffX;
        r.height = 1 - diffY;

        return r;
    }

    #endregion // Handlers
}
