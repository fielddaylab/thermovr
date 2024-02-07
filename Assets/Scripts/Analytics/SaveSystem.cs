using UnityEngine;
using BeauUtil.Debugger;
using BeauData;
using System;
using BeauRoutine;
using System.Collections;

namespace ThermoVR {
    public sealed class SaveSystem : MonoBehaviour {

        [SerializeField] private string m_ServerAddress = null;

        private const string LastUserNameKey = "settings/last-known-profile";

        [SerializeField] private float m_SaveRetryDelay = 10f;

        private PlayerData m_CurrentData;
        [NonSerialized] private string m_UserCode;

        [NonSerialized] private string m_LastKnownProfile;

        private Future<bool> m_SaveOperation;

        public string LastProfileName() { return m_LastKnownProfile; }

        private void Start() {
            OGD.Core.Configure(m_ServerAddress, "THERMOVR");

            // uncomment if we want to retrieve the last known code
            // m_LastKnownProfile = PlayerPrefs.GetString(LastUserNameKey, string.Empty);
        }

        public bool IsServerSave() {
            return !string.IsNullOrEmpty(m_UserCode);
        }

        public void NewLocalSave() {
            PlayerData data = new PlayerData();
            
            DeclareSave(data, null);
        }

        private void DeclareSave(PlayerData data, string userCode) {
            m_CurrentData = data;
            m_UserCode = userCode;
            SetLastKnownSave(userCode);
            // Game.Events.Queue(GameEvents.SaveDeclared, data);
        }

        #region Save

        public Future<bool> NewServerSave(string userCode) {
            if (string.IsNullOrEmpty(userCode)) {
                return Future.Failed<bool>();
            }

            return Future.CreateLinked<bool, string>(ServerDeclareRoutine, userCode, this);
        }

        private IEnumerator ServerDeclareRoutine(Future<bool> future, string userCode) {
            using(var localFuture = Future.Create())
            using(var request = OGD.Player.ClaimId(userCode, null, localFuture.Complete, (f) => localFuture.Fail(f))) {
                yield return localFuture;

                if (localFuture.IsComplete()) {
                    Log.Msg("[DataService] Profile name '{0}' declared to server", userCode);
                } else {
                    string warnMsg = "[DataService] Failed to declare profile name to server: " + localFuture.GetFailure();
                    Log.Warn(warnMsg);
                    // Game.Events.Dispatch(GameEvents.TitleErrorReceived, warnMsg);
                    future.Fail(localFuture.GetFailure());
                    yield break;
                }
            }

            PlayerData newData = new PlayerData();
            string serialized = Serializer.Write(newData, OutputOptions.None, Serializer.Format.Binary);

            using(var localFuture = Future.Create())
            using(var saveRequest = OGD.GameState.PushState(userCode, serialized, localFuture.Complete, (f) => localFuture.Fail(f))) {
                Log.Msg("[SaveSystem] Attempting declare starting save data for id '{0}'", userCode);
                yield return localFuture;

                if (localFuture.IsComplete()) {
                    Log.Msg("[SaveSystem] Saved to server!");
                    DeclareSave(newData, userCode);
                    future.Complete(true);
                } else {
                    string warnMsg = "[SaveSystem] Server save failed.";
                    Log.Warn(warnMsg);
                    // Game.Events.Dispatch(GameEvents.TitleErrorReceived, warnMsg);
                    future.Fail(localFuture.GetFailure());
                }
            }
        }

        public Future<bool> WriteServerSave() {
            if (string.IsNullOrEmpty(m_UserCode)) {
                return Future.Failed<bool>();
            }

            string binarySave = Serializer.Write(m_CurrentData, OutputOptions.None, Serializer.Format.Binary);

            m_SaveOperation?.Cancel();
            return m_SaveOperation = Future.CreateLinked<bool, string, string>(ServerWriteRoutine, m_UserCode, binarySave, this);
        }

        private IEnumerator ServerWriteRoutine(Future<bool> future, string userCode, string saveData) {
            bool saved = false;
            while(!saved) {
                using(var localFuture = Future.Create())
                using(var saveRequest = OGD.GameState.PushState(userCode, saveData, localFuture.Complete, (f) => localFuture.Fail(f))) {
                    Log.Msg("[SaveSystem] Attempting server save with user code '{0}'", userCode);
                    yield return localFuture;

                    if (localFuture.IsComplete()) {
                        Log.Msg("[SaveSystem] Saved to server!");
                        SetLastKnownSave(userCode);
                        saved = true;
                    } else {
                        string warnMsg = "[SaveSystem] Server save failed! Trying again in " + m_SaveRetryDelay + " seconds...";
                        Log.Warn(warnMsg);
                        // Game.Events.Dispatch(GameEvents.TitleErrorReceived, warnMsg);
                        yield return m_SaveRetryDelay;
                    }
                }
            }

            future.Complete(true);
            m_SaveOperation = null;
        }

        public Future<PlayerData> ReadServerSave(string userCode) {
            if (string.IsNullOrEmpty(userCode)) {
                return Future.Failed<PlayerData>();
            }

            return Future.CreateLinked<PlayerData, string>(ServerReadRoutine, userCode, this);
        }

        private IEnumerator ServerReadRoutine(Future<PlayerData> data, string userCode) {
            using(var localFuture = Future.Create<string>())
            using(var loadRequest = OGD.GameState.RequestLatestState(userCode, localFuture.Complete, (f) => localFuture.Fail(f))) {
                Log.Msg("[SaveSystem] Attempting server load with user code '{0}'", userCode);
                yield return localFuture;

                if (localFuture.IsComplete()) {
                    Log.Msg("[SaveSystem] Loaded save from server!");
                    PlayerData deserialized = Serializer.Read<PlayerData>(localFuture.Get());
                    DeclareSave(deserialized, userCode);
                    data.Complete(deserialized);
                } else {
                    string warnMsg = "[SaveSystem] Failed to find profile on server: " + localFuture.GetFailure();
                    Log.Warn(warnMsg);
                    // Game.Events.Dispatch(GameEvents.TitleErrorReceived, warnMsg);
                    data.Fail(localFuture.GetFailure());
                }
            }
        }

        private void SetLastKnownSave(string userCode) {
            m_LastKnownProfile = userCode;
            Log.Msg("[DataService] Profile name {0} set as last known name", m_LastKnownProfile);
            PlayerPrefs.SetString(LastUserNameKey, m_LastKnownProfile ?? string.Empty);
            PlayerPrefs.Save();
        }

        #endregion // Save
    }
}