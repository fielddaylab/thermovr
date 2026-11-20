// Minimal stand-ins for the Unity APIs that the production thermodynamics sources touch.
//
// The real project only ever compiles inside Unity, so there is no assembly we can reference
// from a plain `dotnet test` run. Instead the csproj compiles
// Assets/Scripts/General/IAPWS97.cs and Assets/Scripts/General/ThermoMath.cs directly and
// this file supplies the handful of UnityEngine symbols they need.
//
// UnityEngine.Debug additionally *records* everything it is handed, so tests can assert on the
// diagnostics a failed solve is required to emit (instead of the call being swallowed silently).

using System;
using System.Collections.Generic;

// `using ThermoVR;` appears at the top of ThermoMath.cs.
namespace ThermoVR
{
    internal static class ThermoVRShimAnchor { }
}

// `using ThermoVR.State;` appears at the top of ThermoMath.cs.
namespace ThermoVR.State
{
    internal static class ThermoVRStateShimAnchor { }
}

// `using BeauUtil;` appears at the top of ThermoMath.cs.
namespace BeauUtil
{
    internal static class BeauUtilShimAnchor { }
}

// IAPWS97.cs and ThermoMath.cs carry `[Il2CppSetOption(Option.NullChecks, false)]` attributes.
// These are IL2CPP build hints with no runtime behavior, so empty stand-ins are sufficient.
namespace Unity.IL2CPP.CompilerServices
{
    public enum Option
    {
        NullChecks = 1,
        ArrayBoundsChecks = 2,
        DivideByZeroChecks = 3,
    }

    [AttributeUsage(AttributeTargets.Class | AttributeTargets.Method | AttributeTargets.Property
                    | AttributeTargets.Struct | AttributeTargets.Constructor, AllowMultiple = true)]
    public class Il2CppSetOptionAttribute : Attribute
    {
        public Option Option { get; private set; }
        public object Value { get; private set; }
        public Il2CppSetOptionAttribute(Option option, object value) { Option = option; Value = value; }
    }
}

namespace UnityEngine
{
    public enum LogLevel
    {
        Log,
        Warning,
        Error,
        Exception,
        Assert
    }

    public struct LogEntry
    {
        public LogLevel Level;
        public string Message;

        public override string ToString()
        {
            return Level + ": " + Message;
        }
    }

    public static class Debug
    {
        private static readonly List<LogEntry> s_entries = new List<LogEntry>();

        public static IReadOnlyList<LogEntry> Entries
        {
            get { return s_entries; }
        }

        public static void ClearEntries()
        {
            s_entries.Clear();
        }

        public static List<LogEntry> EntriesOfLevel(LogLevel level)
        {
            List<LogEntry> result = new List<LogEntry>();
            foreach (LogEntry entry in s_entries)
            {
                if (entry.Level == level) { result.Add(entry); }
            }
            return result;
        }

        private static void Record(LogLevel level, object message)
        {
            s_entries.Add(new LogEntry { Level = level, Message = message == null ? "" : message.ToString() });
        }

        public static void Log(object message) { Record(LogLevel.Log, message); }
        public static void Log(object message, object context) { Record(LogLevel.Log, message); }
        public static void LogWarning(object message) { Record(LogLevel.Warning, message); }
        public static void LogWarning(object message, object context) { Record(LogLevel.Warning, message); }
        public static void LogError(object message) { Record(LogLevel.Error, message); }
        public static void LogError(object message, object context) { Record(LogLevel.Error, message); }
        public static void LogException(Exception exception) { Record(LogLevel.Exception, exception); }
        public static void LogException(Exception exception, object context) { Record(LogLevel.Exception, exception); }

        public static void LogFormat(string format, params object[] args) { Record(LogLevel.Log, string.Format(format, args)); }
        public static void LogWarningFormat(string format, params object[] args) { Record(LogLevel.Warning, string.Format(format, args)); }
        public static void LogErrorFormat(string format, params object[] args) { Record(LogLevel.Error, string.Format(format, args)); }

        public static void Assert(bool condition) { if (!condition) { Record(LogLevel.Assert, "assertion failed"); } }
        public static void Assert(bool condition, object message) { if (!condition) { Record(LogLevel.Assert, message); } }
        public static void Break() { }
    }

    public static class Mathf
    {
        public const float PI = (float)Math.PI;
        public const float Infinity = float.PositiveInfinity;
        public const float NegativeInfinity = float.NegativeInfinity;
        public const float Deg2Rad = (float)(Math.PI / 180.0);
        public const float Rad2Deg = (float)(180.0 / Math.PI);
        public static readonly float Epsilon = float.Epsilon;

        public static float Abs(float v) { return Math.Abs(v); }
        public static int Abs(int v) { return Math.Abs(v); }
        public static float Sign(float v) { return v < 0f ? -1f : 1f; }
        public static float Pow(float a, float b) { return (float)Math.Pow(a, b); }
        public static float Sqrt(float v) { return (float)Math.Sqrt(v); }
        public static float Exp(float v) { return (float)Math.Exp(v); }
        public static float Log(float v) { return (float)Math.Log(v); }
        public static float Log(float v, float b) { return (float)Math.Log(v, b); }
        public static float Log10(float v) { return (float)Math.Log10(v); }
        public static float Sin(float v) { return (float)Math.Sin(v); }
        public static float Cos(float v) { return (float)Math.Cos(v); }
        public static float Tan(float v) { return (float)Math.Tan(v); }
        public static float Atan(float v) { return (float)Math.Atan(v); }
        public static float Atan2(float y, float x) { return (float)Math.Atan2(y, x); }
        public static float Floor(float v) { return (float)Math.Floor(v); }
        public static float Ceil(float v) { return (float)Math.Ceiling(v); }
        public static float Round(float v) { return (float)Math.Round(v); }
        public static int FloorToInt(float v) { return (int)Math.Floor(v); }
        public static int CeilToInt(float v) { return (int)Math.Ceiling(v); }
        public static int RoundToInt(float v) { return (int)Math.Round(v); }
        public static float Min(float a, float b) { return Math.Min(a, b); }
        public static int Min(int a, int b) { return Math.Min(a, b); }
        public static float Max(float a, float b) { return Math.Max(a, b); }
        public static int Max(int a, int b) { return Math.Max(a, b); }
        public static float Clamp(float v, float lo, float hi) { return v < lo ? lo : (v > hi ? hi : v); }
        public static int Clamp(int v, int lo, int hi) { return v < lo ? lo : (v > hi ? hi : v); }
        public static float Clamp01(float v) { return Clamp(v, 0f, 1f); }
        public static float Lerp(float a, float b, float t) { return a + (b - a) * Clamp01(t); }
        public static float LerpUnclamped(float a, float b, float t) { return a + (b - a) * t; }
        public static float InverseLerp(float a, float b, float v) { return a == b ? 0f : Clamp01((v - a) / (b - a)); }
        public static bool Approximately(float a, float b) { return Math.Abs(b - a) < Math.Max(1e-6f * Math.Max(Math.Abs(a), Math.Abs(b)), float.Epsilon * 8f); }
    }
}
