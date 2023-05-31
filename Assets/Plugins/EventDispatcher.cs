using System;
using System.Collections;
using System.Collections.Generic;
using BeauUtil;

namespace BeauUtil.Extensions
{
    public sealed class EventDispatcher<TArg>
    {
        #region Types

        /// <summary>
        /// Small struct for holding a queued event.
        /// </summary>
        private struct QueuedEvent
        {
            public readonly StringHash32 Id;
            public readonly TArg Argument;

            public QueuedEvent(StringHash32 id, TArg argument) {
                Id = id;
                Argument = argument;
            }
        }

        /// <summary>
        /// Block of handlers.
        /// </summary>
        private class HandlerBlock
        {
            public StringHash32 EventId;

            private readonly RingBuffer<Handler> m_Handlers = new RingBuffer<Handler>(8, RingBufferMode.Expand);
            private uint[] m_QueuedDeleteMasks = new uint[8];
            private bool m_ForceCleanup;
            private uint m_ExecutionDepth;

            /// <summary>
            /// Returns if this handler block has no handlers.
            /// </summary>
            public bool IsEmpty() {
                return m_Handlers.Count == 0;
            }

            #region Add

            /// <summary>
            /// Registers a handler.
            /// </summary>
            public void Add(Action<TArg> action, UnityEngine.Object binding) {
                m_Handlers.PushBack(new Handler(CastableAction<TArg>.Create(action), binding));
            }

            /// <summary>
            /// Registers a handler.
            /// </summary>
            public void Add(Action action, UnityEngine.Object binding) {
                m_Handlers.PushBack(new Handler(CastableAction<TArg>.Create(action), binding));
            }

            /// <summary>
            /// Registers a handler.
            /// </summary>
            public void Add<U>(Action<U> action, UnityEngine.Object binding) {
                m_Handlers.PushBack(new Handler(CastableAction<TArg>.Create(action), binding));
            }

            #endregion // Add 

            #region Remove

            /// <summary>
            /// Removes all handlers for the given binding.
            /// </summary>
            public void RemoveAll(UnityEngine.Object binding) {
                for (int i = m_Handlers.Count - 1; i >= 0; i--) {
                    if (m_Handlers[i].MatchesBinding(binding)) {
                        if (m_ExecutionDepth > 0) {
                            m_QueuedDeleteMasks[i / 32] |= 1u << (i % 32);
                        }
                        else {
                            m_Handlers.FastRemoveAt(i);
                        }
                    }
                }
            }

            /// <summary>
            /// Removes the first handler for the given method.
            /// </summary>
            public void Remove(Action action) {
                for (int i = 0; i < m_Handlers.Count; i++) {
                    if (m_Handlers[i].MatchesAction(action)) {
                        if (m_ExecutionDepth > 0) {
                            m_QueuedDeleteMasks[i / 32] |= 1u << (i % 32);
                        }
                        else {
                            m_Handlers.FastRemoveAt(i);
                        }
                        return;
                    }
                }
            }

            /// <summary>
            /// Removes the first handler for the given method.
            /// </summary>
            public void Remove(Action<TArg> action) {
                for (int i = 0; i < m_Handlers.Count; i++) {
                    if (m_Handlers[i].MatchesAction(action)) {
                        if (m_ExecutionDepth > 0) {
                            m_QueuedDeleteMasks[i / 32] |= 1u << (i % 32);
                        }
                        else {
                            m_Handlers.FastRemoveAt(i);
                        }
                        return;
                    }
                }
            }

            /// <summary>
            /// Removes the first handler for the given method.
            /// </summary>
            public void Remove<U>(Action<U> action) {
                for (int i = 0; i < m_Handlers.Count; i++) {
                    if (m_Handlers[i].MatchesAction(action)) {
                        if (m_ExecutionDepth > 0) {
                            m_QueuedDeleteMasks[i / 32] |= 1u << (i % 32);
                        }
                        else {
                            m_Handlers.FastRemoveAt(i);
                        }
                        return;
                    }
                }
            }

            #endregion // Remove

            /// <summary>
            /// Invokes the event on all registered handlers.
            /// </summary>
            public void Invoke(TArg context) {
                int count = m_Handlers.Count;

                if (count == 0) {
                    return;
                }

                m_ExecutionDepth++;

                for (int i = 0; i < count; i++) {
                    Handler handler = m_Handlers[i];
                    // if this handler is not queued for removal
                    if ((m_QueuedDeleteMasks[i / 32] & (1 << (i % 32))) == 0) {
                        handler.Invoke(context);
                    }
                }

                if (--m_ExecutionDepth == 1) {
                    if (HasDeletesQueued()) {
                        Cleanup();
                    }
                }
            }

            /// <summary>
            /// Cleans up queued deletes.
            /// </summary>
            public void Cleanup() {
                if (m_ExecutionDepth > 0) {
                    m_ForceCleanup = true;
                }
                else {
                    for (int i = m_Handlers.Count - 1; i >= 0; i--) {
                        if ((m_QueuedDeleteMasks[i / 32] & (1 << (i % 32))) != 0) {
                            m_Handlers.FastRemoveAt(i);
                        }
                    }

                    for (int i = 0; i < m_QueuedDeleteMasks.Length; i++) {
                        m_QueuedDeleteMasks[i] = 0;
                    }
                }
            }

            /// <summary>
            /// Clears all handlers.
            /// </summary>
            public void Clear() {
                if (m_ExecutionDepth > 0) {
                    for (int i = 0; i < m_QueuedDeleteMasks.Length; i++) {
                        m_QueuedDeleteMasks[i] = uint.MaxValue;
                    }
                }
                else {
                    m_Handlers.Clear();
                    m_ForceCleanup = false;
                    for (int i = 0; i < m_QueuedDeleteMasks.Length; i++) {
                        m_QueuedDeleteMasks[i] = 0;
                    }
                }
            }

            /// <summary>
            /// Returns if any deletes are queued on this handler block.
            /// </summary>
            private bool HasDeletesQueued() {
                if (m_ForceCleanup) {
                    return true;
                }

                for (int i = 0; i < m_QueuedDeleteMasks.Length; i++) {
                    if (m_QueuedDeleteMasks[i] != 0) {
                        return true;
                    }
                }

                return false;
            }
        }

        /// <summary>
        /// Struct for holding a method and its binding.
        /// </summary>
        private struct Handler
        {
            private int m_BindingInstanceId;
            private CastableAction<TArg> m_Action;

            public Handler(CastableAction<TArg> action, UnityEngine.Object binding) {
                m_BindingInstanceId = !ReferenceEquals(binding, null) ? binding.GetInstanceID() : 0;
                m_Action = action;
            }

            public void Invoke(TArg context) {
                m_Action.Invoke(context);
            }

            public bool MatchesBinding(UnityEngine.Object binding) {
                return !ReferenceEquals(binding, null) ? m_BindingInstanceId == binding.GetInstanceID() : m_BindingInstanceId == 0;
            }

            public bool MatchesAction(Action action) {
                return m_Action.Equals(action);
            }

            public bool MatchesAction(Action<TArg> action) {
                return m_Action.Equals(action);
            }

            public bool MatchesAction(MulticastDelegate action) {
                return m_Action.Equals(action);
            }
        }

        #endregion // Types

        private readonly Dictionary<StringHash32, HandlerBlock> m_Handlers = new Dictionary<StringHash32, HandlerBlock>();
        private readonly RingBuffer<HandlerBlock> m_FreeHandlers = new RingBuffer<HandlerBlock>();
        private readonly RingBuffer<QueuedEvent> m_QueuedEvents = new RingBuffer<QueuedEvent>(16, RingBufferMode.Expand);
        private readonly List<StringHash32> m_TempEventList = new List<StringHash32>();

        #region Registration

        /// <summary>
        /// Registers an event handler, optionally bound to a given object.
        /// </summary>
        public EventDispatcher<TArg> Register(StringHash32 eventId, Action inAction, UnityEngine.Object binding = null) {
            HandlerBlock block = GetBlock(eventId, true);
            block.Add(inAction, binding);
            return this;
        }

        /// <summary>
        /// Registers an event handler, optionally bound to a given object.
        /// </summary>
        public EventDispatcher<TArg> Register(StringHash32 eventId, Action<TArg> inActionWithContext, UnityEngine.Object binding = null) {
            HandlerBlock block = GetBlock(eventId, true);
            block.Add(inActionWithContext, binding);
            return this;
        }

        /// <summary>
        /// Registers an event handler, optionally bound to a given object.
        /// </summary>
        public EventDispatcher<TArg> Register<U>(StringHash32 eventId, Action<U> inActionWithCastedContext, UnityEngine.Object binding = null) {
            HandlerBlock block = GetBlock(eventId, true);
            block.Add(inActionWithCastedContext, binding);
            return this;
        }

        /// <summary>
        /// Deregisters an event handler.
        /// </summary>
        public EventDispatcher<TArg> Deregister(StringHash32 eventId, Action action) {
            HandlerBlock block;
            if (m_Handlers.TryGetValue(eventId, out block)) {
                block.Remove(action);
                if (block.IsEmpty()) {
                    m_Handlers.Remove(eventId);
                    m_FreeHandlers.PushBack(block);
                }
            }

            return this;
        }

        /// <summary>
        /// Deregisters an event handler.
        /// </summary>
        public EventDispatcher<TArg> Deregister(StringHash32 eventId, Action<TArg> action) {
            HandlerBlock block;
            if (m_Handlers.TryGetValue(eventId, out block)) {
                block.Remove(action);
                if (block.IsEmpty()) {
                    m_Handlers.Remove(eventId);
                    m_FreeHandlers.PushBack(block);
                }
            }

            return this;
        }

        /// <summary>
        /// Deregisters an event handler.
        /// </summary>
        public EventDispatcher<TArg> Deregister<T>(StringHash32 eventId, Action<T> castedAction) {
            HandlerBlock block;
            if (m_Handlers.TryGetValue(eventId, out block)) {
                block.Remove(castedAction);
                if (block.IsEmpty()) {
                    m_Handlers.Remove(eventId);
                    m_FreeHandlers.PushBack(block);
                }
            }

            return this;
        }

        /// <summary>
        /// Deregisters all handlers for the given event.
        /// </summary>
        public EventDispatcher<TArg> DeregisterAll(StringHash32 eventId) {
            HandlerBlock block;
            if (m_Handlers.TryGetValue(eventId, out block)) {
                block.Clear();
                m_Handlers.Remove(eventId);
                m_FreeHandlers.PushBack(block);
            }

            return this;
        }

        /// <summary>
        /// Deregisters all handlers associated with the given binding.
        /// </summary>
        public EventDispatcher<TArg> DeregisterAll(UnityEngine.Object binding) {
            if (binding.IsReferenceNull()) {
                return this;
            }

            m_TempEventList.Clear();
            foreach (var block in m_Handlers.Values) {
                block.RemoveAll(binding);
                if (block.IsEmpty()) {
                    m_FreeHandlers.PushBack(block);
                    m_TempEventList.Add(block.EventId);
                }
            }

            foreach (var tempEventId in m_TempEventList) {
                m_Handlers.Remove(tempEventId);
            }

            m_TempEventList.Clear();
            return this;
        }

        private HandlerBlock GetBlock(StringHash32 eventId, bool createIfNotThere) {
            HandlerBlock block;
            if (!m_Handlers.TryGetValue(eventId, out block) && createIfNotThere) {
                if (m_FreeHandlers.Count > 0) {
                    block = m_FreeHandlers.PopBack();
                }
                else {
                    block = new HandlerBlock();
                }
                block.EventId = eventId;
                m_Handlers.Add(eventId, block);
            }
            return block;
        }

        #endregion // Registration

        #region Operations

        /// <summary>
        /// Dispatches the event to its handlers.
        /// </summary>
        public void Dispatch(StringHash32 eventId, TArg argument = default) {
            HandlerBlock block;
            if (m_Handlers.TryGetValue(eventId, out block)) {
                block.Invoke(argument);
            }
        }

        /// <summary>
        /// Dispatches the event to its handlers the next time FlushQueue() is called.
        /// </summary>
        public void Queue(StringHash32 eventId, TArg argument = default) {
            m_QueuedEvents.PushBack(new QueuedEvent(eventId, argument));
        }

        /// <summary>
        /// Dispatches all events queued up with Queue.
        /// </summary>
        public void FlushQueue() {
            QueuedEvent evt;
            while (m_QueuedEvents.TryPopFront(out evt)) {
                Dispatch(evt.Id, evt.Argument);
            }
        }

        #endregion // Operations

        #region Helpers

        /// <summary>
        /// Coroutine/iterator. Waits for an event with the given id to dispatch before continuing.
        /// </summary>
        public IEnumerator Wait(StringHash32 eventId) {
            return new WaitIterator(this, eventId, null);
        }

        /// <summary>
        /// Coroutine/iterator. Waits for an event with the given id to dispatch before continuing.
        /// </summary>
        public IEnumerator Wait(StringHash32 eventId, Action callback) {
            return new WaitIterator(this, eventId, CastableAction<TArg>.Create(callback));
        }

        /// <summary>
        /// Coroutine/iterator. Waits for an event with the given id to dispatch before continuing.
        /// </summary>
        public IEnumerator Wait(StringHash32 eventId, Action<TArg> callback) {
            return new WaitIterator(this, eventId, CastableAction<TArg>.Create(callback));
        }

        /// <summary>
        /// Coroutine/iterator. Waits for an event with the given id to dispatch before continuing.
        /// </summary>
        public IEnumerator Wait<T>(StringHash32 eventId, Action<T> callback) {
            return new WaitIterator(this, eventId, CastableAction<TArg>.Create(callback));
        }

        private class WaitIterator : IEnumerator, IDisposable
        {
            private const int Phase_Init = 0;
            private const int Phase_Wait = 1;
            private const int Phase_Done = 2;

            private EventDispatcher<TArg> m_Parent;
            private readonly StringHash32 m_EventId;
            private Action<TArg> m_InnerCallback;
            private CastableAction<TArg>? m_CustomCallback;
            private int m_Phase = 0;

            public WaitIterator(EventDispatcher<TArg> parent, StringHash32 eventId, CastableAction<TArg>? customCallback) {
                m_Parent = parent;
                m_EventId = eventId;
                m_InnerCallback = InnerCallback;
                m_CustomCallback = customCallback;
            }

            public object Current { get { return null; } }

            public bool MoveNext() {
                switch (m_Phase) {
                    case Phase_Init: {
                            m_Parent.Register(m_EventId, m_InnerCallback);
                            m_Phase++;
                            return true;
                        }

                    case Phase_Wait: {
                            return true;
                        }

                    case Phase_Done:
                    default: {
                            return false;
                        }
                }
            }

            public void Reset() {
                throw new NotSupportedException();
            }

            public void Dispose() {
                if (m_Phase == Phase_Wait) {
                    m_Parent.Deregister(m_EventId, m_InnerCallback);
                }
                m_Phase = Phase_Done;
                m_CustomCallback = null;
                m_Parent = null;
                m_InnerCallback = null;
            }

            private void InnerCallback(TArg arg) {
                m_Phase = Phase_Done;
                m_Parent.Deregister(m_EventId, m_InnerCallback);
                m_CustomCallback?.Invoke(arg);
                m_InnerCallback = null;
                m_Parent = null;
                m_CustomCallback = null;
            }
        }

        #endregion // Helpers
    }
}