using System;
using System.Collections.Generic;

namespace ThermoVR.Tests
{
    internal sealed class TestRunner
    {
        private readonly List<KeyValuePair<string, Action>> m_tests =
            new List<KeyValuePair<string, Action>>();

        public void Add(string name, Action body) {
            m_tests.Add(new KeyValuePair<string, Action>(name, body));
        }

        public int Run() {
            int passed = 0;
            List<string> failures = new List<string>();

            foreach (KeyValuePair<string, Action> test in m_tests) {
                try {
                    test.Value();
                    passed++;
                    Console.WriteLine("PASS  " + test.Key);
                }
                catch (Exception e) {
                    failures.Add(test.Key + ": " + e.Message);
                    Console.WriteLine("FAIL  " + test.Key + "  --  " + e.Message);
                }
            }

            Console.WriteLine();
            Console.WriteLine(passed + " passed, " + failures.Count + " failed, "
                + m_tests.Count + " total");

            return failures.Count == 0 ? 0 : 1;
        }
    }
}
