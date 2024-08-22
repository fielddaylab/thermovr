using BeauRoutine;
using BeauUtil;
using BeauUtil.Debugger;
using BeauUtil.Extensions;
using System.Collections;
using System.Collections.Generic;
using ThermoVR.UI;
using UnityEngine;

namespace ThermoVR
{
    public class ViewportMgr : Singleton<ViewportMgr>
    {
        static public RenderMgr RenderMgr { get; internal set; }

        protected override void Awake()
        {
            base.Awake();

            RenderMgr = new RenderMgr();
            RenderMgr.Initialize();

            RenderMgr.EnableAspectClamping(1920, 1080);
        }

        private void FixedUpdate()
        {
            EventMgr.Events.Dispatch(GameEvents.CanvasPreUpdate);
            EventMgr.Events.Dispatch(GameEvents.ApplicationPreRender);
        }

        private void LateUpdate()
        {
            RenderMgr.PollScreenSettings();
        }
    }
}