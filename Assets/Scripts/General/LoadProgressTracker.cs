/*
Tracks coarse-grained loading progress during the initial game init.

Progress is derived ONLY from completed work units (phases / loop chunks), never from
elapsed time, so the reported fraction is monotonic non-decreasing and reaches exactly
1.0f when the init chain reports it is done.

Deliberately plain C# (no UnityEngine dependency) so the progress math can be reasoned
about / unit tested outside of the Unity editor.
*/

using System;

namespace ThermoVR
{
    public static class LoadProgressTracker
    {
        private static int s_completed = 0;
        private static int s_total = 1;
        private static bool s_started = false;

        /// <summary>
        /// Begins a load with a fixed denominator of work units. The denominator must be
        /// known up front so that Progress can never move backwards when later phases turn
        /// out to need fewer chunks than their upper bound.
        /// </summary>
        public static void Begin(int totalUnits) {
            s_total = Math.Max(1, totalUnits);
            s_completed = 0;
            s_started = true;
        }

        /// <summary>
        /// Marks additional work units as completed.
        /// </summary>
        public static void Complete(int units = 1) {
            if (units <= 0) {
                // never move backwards
                return;
            }
            s_completed = Math.Min(s_total, s_completed + units);
        }

        /// <summary>
        /// Moves the completed-unit count up to an absolute value. Never decreases it.
        /// Used when a chunked loop finishes early (fewer chunks than its upper bound) and
        /// we want to "catch up" to the budget that was reserved for it.
        /// </summary>
        public static void CompleteTo(int completedUnits) {
            if (completedUnits <= s_completed) {
                return;
            }
            s_completed = Math.Min(s_total, completedUnits);
        }

        /// <summary>
        /// Marks the whole load as finished. Progress becomes exactly 1.0f.
        /// </summary>
        public static void Finish() {
            s_started = true;
            s_completed = s_total;
        }

        /// <summary>
        /// Resets back to "nothing loaded yet".
        /// </summary>
        public static void Reset() {
            s_completed = 0;
            s_total = 1;
            s_started = false;
        }

        /// <summary>
        /// 0..1, monotonic non-decreasing, exactly 1.0f once all units are done.
        /// </summary>
        public static float Progress {
            get {
                if (s_completed >= s_total) {
                    return 1.0f;
                }
                if (s_completed <= 0) {
                    return 0.0f;
                }
                return (float)s_completed / (float)s_total;
            }
        }

        /// <summary>
        /// True once Begin()/Finish() has been called at least once.
        /// </summary>
        public static bool HasStarted {
            get { return s_started; }
        }

        /// <summary>
        /// True only once every reserved work unit has been completed.
        /// </summary>
        public static bool IsComplete {
            get { return s_started && s_completed >= s_total; }
        }

        public static int CompletedUnits {
            get { return s_completed; }
        }

        public static int TotalUnits {
            get { return s_total; }
        }
    }
}
