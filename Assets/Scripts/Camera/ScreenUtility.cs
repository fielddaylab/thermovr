#if UNITY_WEBGL && !UNITY_EDITOR
#define USE_JSLIB
#endif // UNITY_WEBGL && !UNITY_EDITOR

using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using UnityEngine;

static public class ScreenUtility
{
#if UNITY_WEBGL && !UNITY_EDITOR

        [DllImport("__Internal")]
        static private extern void NativeFullscreen_SetFullscreen(bool fullscreen);

#endif // UNITY_WEBGL && !UNITY_EDITOR

    /// <summary>
    /// Sets the fullscreen mode.
    /// </summary>
    static public void SetFullscreen(bool fullscreen)
    {
#if USE_JSLIB
            NativeFullscreen_SetFullscreen(fullscreen);
#elif UNITY_EDITOR
        // TODO: fullscreen within editor?
        Screen.fullScreen = fullscreen;
#else
            Screen.fullScreen = fullscreen;
#endif // UNITY_WEBGL && !UNITY_EDITOR
    }

    /// <summary>
    /// Returns the fullscreen mode.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    static public bool GetFullscreen()
    {
        return Screen.fullScreen;
    }

    /// <summary>
    /// Returns the current resolution.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    static public Resolution GetResolution()
    {
#if UNITY_EDITOR || UNITY_WEBGL
        return new Resolution()
        {
            width = Screen.width,
            height = Screen.height,
            refreshRate = 60
        };
#else
            return Screen.currentResolution;
#endif // UNITY_EDITOR || UNITY_WEBGL
    }
}
