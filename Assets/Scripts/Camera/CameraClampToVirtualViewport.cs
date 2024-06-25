using System;
using BeauUtil;
using ThermoVR;
using UnityEngine;

/// <summary>
/// Sets the given camera to clamp to the virtual viewport.
/// </summary>
[DisallowMultipleComponent, RequireComponent(typeof(Camera))]
[DefaultExecutionOrder(-10000)]
public sealed class CameraClampToVirtualViewport : MonoBehaviour
{
    public Rect Viewport = new Rect(0, 0, 1, 1);
    [NonSerialized] private Camera m_Camera;

    public Camera Camera
    {
        get { return this.CacheComponent(ref m_Camera); }
    }

    private void Start()
    {
        this.CacheComponent(ref m_Camera);
        GameMgr.RenderMgr.AddClampedViewportCamera(this);
    }

    private void OnDisable()
    {
        GameMgr.RenderMgr.RemoveClampedViewportCamera(this);
    }
}