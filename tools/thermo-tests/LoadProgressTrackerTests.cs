// Behavioural unit tests for ThermoVR.LoadProgressTracker, the pure-C# progress math
// behind ThermoPresent.LoadProgress / GameMgr.LoadProgress / UILoading.LoadProgress.
//
// Acceptance criteria covered here:
//   * progress is derived from completed work units, never from elapsed time
//   * progress is monotonic non-decreasing
//   * progress stays inside [0,1] and reaches exactly 1.0f when init finishes
//   * the denominator is fixed up front by Begin(), so late phases cannot rewind it

using System;
using System.Collections.Generic;
using System.Threading;
using ThermoVR;

namespace ThermoVR.Tests
{
    internal static class LoadProgressTrackerTests
    {
        public static void Register(TestRunner runner) {
            runner.Add("Tracker.Progress_IsZero_BeforeAnyUnitCompletes", Progress_IsZero_BeforeAnyUnitCompletes);
            runner.Add("Tracker.Begin_FixesDenominator", Begin_FixesDenominator);
            runner.Add("Tracker.Begin_ClampsNonPositiveTotalToOne", Begin_ClampsNonPositiveTotalToOne);
            runner.Add("Tracker.Progress_IsFractionOfCompletedUnits", Progress_IsFractionOfCompletedUnits);
            runner.Add("Tracker.Progress_ReachesExactlyOne_WhenAllUnitsDone", Progress_ReachesExactlyOne_WhenAllUnitsDone);
            runner.Add("Tracker.Complete_ClampsAtTotal", Complete_ClampsAtTotal);
            runner.Add("Tracker.Complete_NonPositive_IsNoOp", Complete_NonPositive_IsNoOp);
            runner.Add("Tracker.CompleteTo_CatchesUp_ButNeverRewinds", CompleteTo_CatchesUp_ButNeverRewinds);
            runner.Add("Tracker.Finish_ForcesExactlyOne", Finish_ForcesExactlyOne);
            runner.Add("Tracker.IsComplete_OnlyAfterAllUnitsDone", IsComplete_OnlyAfterAllUnitsDone);
            runner.Add("Tracker.HasStarted_ReflectsBeginAndReset", HasStarted_ReflectsBeginAndReset);
            runner.Add("Tracker.Reset_ReturnsToZeroProgress", Reset_ReturnsToZeroProgress);
            runner.Add("Tracker.Progress_IsNotTimeBased", Progress_IsNotTimeBased);
            runner.Add("Tracker.Progress_IsMonotonic_AcrossRandomSequence", Progress_IsMonotonic_AcrossRandomSequence);
            runner.Add("Tracker.SimulatedInitChain_IsMonotonicAndEndsAtOne", SimulatedInitChain_IsMonotonicAndEndsAtOne);
            runner.Add("Tracker.OverReportingUnits_CannotExceedOne", OverReportingUnits_CannotExceedOne);
        }

        private static void Progress_IsZero_BeforeAnyUnitCompletes() {
            LoadProgressTracker.Reset();
            Check.FloatExact(0.0f, LoadProgressTracker.Progress, "progress after Reset");

            LoadProgressTracker.Begin(10);
            Check.FloatExact(0.0f, LoadProgressTracker.Progress, "progress right after Begin");
            Check.IntEqual(0, LoadProgressTracker.CompletedUnits, "completed units after Begin");
        }

        private static void Begin_FixesDenominator() {
            LoadProgressTracker.Reset();
            LoadProgressTracker.Begin(37);
            Check.IntEqual(37, LoadProgressTracker.TotalUnits, "total units");

            LoadProgressTracker.Complete(10);
            Check.IntEqual(37, LoadProgressTracker.TotalUnits, "total units unchanged by Complete");
            Check.FloatNear(10.0f / 37.0f, LoadProgressTracker.Progress, 1e-6f, "progress at 10/37");
        }

        private static void Begin_ClampsNonPositiveTotalToOne() {
            LoadProgressTracker.Reset();
            LoadProgressTracker.Begin(0);
            Check.IntEqual(1, LoadProgressTracker.TotalUnits, "total units for Begin(0)");
            Check.FloatExact(0.0f, LoadProgressTracker.Progress, "progress for Begin(0)");

            LoadProgressTracker.Begin(-5);
            Check.IntEqual(1, LoadProgressTracker.TotalUnits, "total units for Begin(-5)");
            Check.FloatExact(0.0f, LoadProgressTracker.Progress, "progress for Begin(-5)");

            LoadProgressTracker.Complete();
            Check.FloatExact(1.0f, LoadProgressTracker.Progress, "progress after single unit of 1");
        }

        private static void Progress_IsFractionOfCompletedUnits() {
            LoadProgressTracker.Reset();
            LoadProgressTracker.Begin(4);

            Check.FloatExact(0.0f, LoadProgressTracker.Progress, "0/4");
            LoadProgressTracker.Complete();
            Check.FloatExact(0.25f, LoadProgressTracker.Progress, "1/4");
            LoadProgressTracker.Complete();
            Check.FloatExact(0.5f, LoadProgressTracker.Progress, "2/4");
            LoadProgressTracker.Complete(1);
            Check.FloatExact(0.75f, LoadProgressTracker.Progress, "3/4");
        }

        private static void Progress_ReachesExactlyOne_WhenAllUnitsDone() {
            LoadProgressTracker.Reset();
            // 3 is not representable as an exact binary fraction denominator, so this
            // also guards against 0.9999999 style float drift at the end of the load.
            LoadProgressTracker.Begin(3);
            LoadProgressTracker.Complete();
            LoadProgressTracker.Complete();
            Check.True(LoadProgressTracker.Progress < 1.0f, "progress below 1 at 2/3");
            LoadProgressTracker.Complete();
            Check.FloatExact(1.0f, LoadProgressTracker.Progress, "progress at 3/3");
        }

        private static void Complete_ClampsAtTotal() {
            LoadProgressTracker.Reset();
            LoadProgressTracker.Begin(5);
            LoadProgressTracker.Complete(100);
            Check.IntEqual(5, LoadProgressTracker.CompletedUnits, "completed units clamped to total");
            Check.FloatExact(1.0f, LoadProgressTracker.Progress, "progress clamped to 1");
        }

        private static void Complete_NonPositive_IsNoOp() {
            LoadProgressTracker.Reset();
            LoadProgressTracker.Begin(8);
            LoadProgressTracker.Complete(3);
            float before = LoadProgressTracker.Progress;

            LoadProgressTracker.Complete(0);
            Check.FloatExact(before, LoadProgressTracker.Progress, "progress after Complete(0)");
            LoadProgressTracker.Complete(-4);
            Check.FloatExact(before, LoadProgressTracker.Progress, "progress after Complete(-4)");
            Check.IntEqual(3, LoadProgressTracker.CompletedUnits, "completed units unchanged");
        }

        private static void CompleteTo_CatchesUp_ButNeverRewinds() {
            LoadProgressTracker.Reset();
            LoadProgressTracker.Begin(10);
            LoadProgressTracker.Complete(4);

            LoadProgressTracker.CompleteTo(7);
            Check.IntEqual(7, LoadProgressTracker.CompletedUnits, "CompleteTo catches up to reserved budget");

            LoadProgressTracker.CompleteTo(2);
            Check.IntEqual(7, LoadProgressTracker.CompletedUnits, "CompleteTo cannot rewind");

            LoadProgressTracker.CompleteTo(1000);
            Check.IntEqual(10, LoadProgressTracker.CompletedUnits, "CompleteTo clamps to total");
            Check.FloatExact(1.0f, LoadProgressTracker.Progress, "progress clamped to 1");
        }

        private static void Finish_ForcesExactlyOne() {
            LoadProgressTracker.Reset();
            LoadProgressTracker.Begin(1234);
            LoadProgressTracker.Complete(7);
            Check.True(LoadProgressTracker.Progress < 1.0f, "progress below 1 mid-load");

            LoadProgressTracker.Finish();
            Check.FloatExact(1.0f, LoadProgressTracker.Progress, "progress after Finish");
            Check.True(LoadProgressTracker.IsComplete, "IsComplete after Finish");
        }

        private static void IsComplete_OnlyAfterAllUnitsDone() {
            LoadProgressTracker.Reset();
            Check.False(LoadProgressTracker.IsComplete, "IsComplete before any load started");

            LoadProgressTracker.Begin(3);
            Check.False(LoadProgressTracker.IsComplete, "IsComplete at 0/3");
            LoadProgressTracker.Complete(2);
            Check.False(LoadProgressTracker.IsComplete, "IsComplete at 2/3");
            LoadProgressTracker.Complete();
            Check.True(LoadProgressTracker.IsComplete, "IsComplete at 3/3");
        }

        private static void HasStarted_ReflectsBeginAndReset() {
            LoadProgressTracker.Reset();
            Check.False(LoadProgressTracker.HasStarted, "HasStarted after Reset");

            LoadProgressTracker.Begin(2);
            Check.True(LoadProgressTracker.HasStarted, "HasStarted after Begin");

            LoadProgressTracker.Reset();
            Check.False(LoadProgressTracker.HasStarted, "HasStarted after second Reset");

            LoadProgressTracker.Finish();
            Check.True(LoadProgressTracker.HasStarted, "HasStarted after Finish");
        }

        private static void Reset_ReturnsToZeroProgress() {
            LoadProgressTracker.Reset();
            LoadProgressTracker.Begin(6);
            LoadProgressTracker.Complete(6);
            Check.FloatExact(1.0f, LoadProgressTracker.Progress, "progress before Reset");

            LoadProgressTracker.Reset();
            Check.FloatExact(0.0f, LoadProgressTracker.Progress, "progress after Reset");
            Check.IntEqual(0, LoadProgressTracker.CompletedUnits, "completed units after Reset");
            Check.IntEqual(1, LoadProgressTracker.TotalUnits, "total units after Reset");
            Check.False(LoadProgressTracker.IsComplete, "IsComplete after Reset");
        }

        private static void Progress_IsNotTimeBased() {
            LoadProgressTracker.Reset();
            LoadProgressTracker.Begin(20);
            LoadProgressTracker.Complete(5);

            float before = LoadProgressTracker.Progress;
            Thread.Sleep(60);
            Check.FloatExact(before, LoadProgressTracker.Progress,
                "progress must not advance with wall-clock time");

            // and it still only advances when work units are reported
            LoadProgressTracker.Complete();
            Check.True(LoadProgressTracker.Progress > before, "progress advances on completed work");
        }

        private static void Progress_IsMonotonic_AcrossRandomSequence() {
            // Deterministic seed: a failure here must be reproducible.
            Random rng = new Random(20251103);

            for (int trial = 0; trial < 50; trial++) {
                int total = 1 + rng.Next(200);
                LoadProgressTracker.Reset();
                LoadProgressTracker.Begin(total);

                float last = LoadProgressTracker.Progress;
                Check.FloatExact(0.0f, last, "trial start progress");

                for (int step = 0; step < 400; step++) {
                    int roll = rng.Next(4);
                    if (roll == 0) {
                        LoadProgressTracker.Complete(rng.Next(-3, 4));
                    }
                    else if (roll == 1) {
                        LoadProgressTracker.CompleteTo(rng.Next(-5, total + 5));
                    }
                    else {
                        LoadProgressTracker.Complete();
                    }

                    float now = LoadProgressTracker.Progress;
                    Check.True(now >= last, "progress must never decrease (was "
                        + last + ", now " + now + ")");
                    Check.True(now >= 0.0f && now <= 1.0f, "progress must stay in [0,1], got " + now);
                    last = now;
                }

                LoadProgressTracker.Finish();
                Check.True(LoadProgressTracker.Progress >= last, "Finish must not decrease progress");
                Check.FloatExact(1.0f, LoadProgressTracker.Progress, "progress at end of trial");
            }
        }

        // Mirrors the real init chain's unit accounting:
        //   ThermoPresent.CountInitWorkUnits() = 3 fixed (findObjects, state.reset,
        //      visualize_state) + CountGenMeshWorkUnits()
        //   CountGenMeshWorkUnits() = surface chunks + cull chunks + 3 fixed
        //      (dome, tri_index, mesh_upload)
        //   GameMgr adds 3 more (world.init, thermo_present.reset, dispatch)
        // The denominator is computed once, before Begin(), so the reported fraction can
        // only ever move forwards.
        private static void SimulatedInitChain_IsMonotonicAndEndsAtOne() {
            int[] cull_chunk_counts = new int[] { 1, 3, 9, 64 };
            const int samples = 350;
            const int surface_rows_per_chunk = 25;
            int surface_chunks = (samples + surface_rows_per_chunk - 1) / surface_rows_per_chunk;
            Check.IntEqual(14, surface_chunks, "surface chunk count for samples=350, rows=25");

            foreach (int cull_chunks in cull_chunk_counts) {
                int gen_mesh_units = surface_chunks + cull_chunks + 3;
                int thermo_units = 3 + gen_mesh_units;
                int total_units = thermo_units + 3;

                LoadProgressTracker.Reset();
                LoadProgressTracker.Begin(total_units);
                Check.IntEqual(total_units, LoadProgressTracker.TotalUnits, "fixed denominator");

                List<float> samples_seen = new List<float>();
                samples_seen.Add(LoadProgressTracker.Progress);

                // ThermoPresent.InitRoutine
                LoadProgressTracker.Complete();                 // findObjects
                samples_seen.Add(LoadProgressTracker.Progress);

                for (int i = 0; i < surface_chunks; i++) {      // chunked surface sweep
                    LoadProgressTracker.Complete();
                    samples_seen.Add(LoadProgressTracker.Progress);
                }
                LoadProgressTracker.Complete();                 // dome
                LoadProgressTracker.Complete();                 // tri_index
                samples_seen.Add(LoadProgressTracker.Progress);

                for (int i = 0; i < cull_chunks; i++) {         // chunked cull
                    LoadProgressTracker.Complete();
                    samples_seen.Add(LoadProgressTracker.Progress);
                }
                LoadProgressTracker.Complete();                 // mesh_upload
                LoadProgressTracker.Complete();                 // state.reset
                LoadProgressTracker.Complete();                 // visualize_state
                samples_seen.Add(LoadProgressTracker.Progress);

                // GameMgr.InitRoutine
                LoadProgressTracker.Complete();                 // world.init
                LoadProgressTracker.Complete();                 // thermo_present.reset
                samples_seen.Add(LoadProgressTracker.Progress);

                Check.True(LoadProgressTracker.Progress < 1.0f,
                    "progress must still be below 1 before dispatch (cull_chunks=" + cull_chunks + ")");
                Check.False(LoadProgressTracker.IsComplete, "not complete before dispatch");

                LoadProgressTracker.Finish();                   // dispatch
                samples_seen.Add(LoadProgressTracker.Progress);

                for (int i = 1; i < samples_seen.Count; i++) {
                    Check.True(samples_seen[i] >= samples_seen[i - 1],
                        "init chain progress must be monotonic at step " + i
                        + " (cull_chunks=" + cull_chunks + ")");
                    Check.True(samples_seen[i] >= 0.0f && samples_seen[i] <= 1.0f,
                        "init chain progress in [0,1] at step " + i);
                }

                Check.FloatExact(0.0f, samples_seen[0], "init chain starts at 0");
                Check.FloatExact(1.0f, samples_seen[samples_seen.Count - 1], "init chain ends at exactly 1");
                Check.True(LoadProgressTracker.IsComplete, "IsComplete at end of init chain");
            }
        }

        // If a chunked loop reports more chunks than were reserved (e.g. the cull loop
        // sweeps more triangles than the up-front estimate), progress must saturate at 1
        // rather than overshoot.
        private static void OverReportingUnits_CannotExceedOne() {
            LoadProgressTracker.Reset();
            LoadProgressTracker.Begin(12);

            float last = 0.0f;
            for (int i = 0; i < 100; i++) {
                LoadProgressTracker.Complete();
                float now = LoadProgressTracker.Progress;
                Check.True(now >= last, "monotonic while over-reporting");
                Check.True(now <= 1.0f, "progress never exceeds 1 while over-reporting, got " + now);
                last = now;
            }
            Check.FloatExact(1.0f, LoadProgressTracker.Progress, "saturated progress");
            Check.IntEqual(12, LoadProgressTracker.CompletedUnits, "completed units saturate at total");
        }
    }
}
