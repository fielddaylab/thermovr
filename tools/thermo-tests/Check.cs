// Minimal assertion helpers + test runner.
//
// This repository has no test framework available on the build host (no network, no
// NuGet cache), so the suite is a plain console app. Each registered test is a named
// Action; a thrown AssertionException marks it failed, and the process exit code is
// non-zero if any test failed.

using System;
using System.Collections.Generic;
using System.Globalization;

namespace ThermoVR.Tests
{
    internal sealed class AssertionException : Exception
    {
        public AssertionException(string message) : base(message) { }
    }

    internal static class Check
    {
        public static void True(bool condition, string what) {
            if (!condition) {
                throw new AssertionException("expected true: " + what);
            }
        }

        public static void False(bool condition, string what) {
            if (condition) {
                throw new AssertionException("expected false: " + what);
            }
        }

        public static void IntEqual(int expected, int actual, string what) {
            if (expected != actual) {
                throw new AssertionException(string.Format(CultureInfo.InvariantCulture,
                    "{0}: expected {1}, got {2}", what, expected, actual));
            }
        }

        // Bit-exact float comparison, used where the card demands exactly 1.0f.
        public static void FloatExact(float expected, float actual, string what) {
            if (!expected.Equals(actual)) {
                throw new AssertionException(string.Format(CultureInfo.InvariantCulture,
                    "{0}: expected exactly {1}, got {2}", what, expected, actual));
            }
        }

        public static void FloatNear(float expected, float actual, float tol, string what) {
            if (Math.Abs(expected - actual) > tol) {
                throw new AssertionException(string.Format(CultureInfo.InvariantCulture,
                    "{0}: expected {1} (+/- {2}), got {3}", what, expected, tol, actual));
            }
        }
    }
}
