using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace ThermoVR.Lab
{
    public class CartridgeMgr : MonoBehaviour
    {
        [SerializeField] private Cartridge[] m_initialCartridges;
        [SerializeField] private GameObject m_cartridgePrefab;
        [SerializeField] private float m_cartridgeSpacing;

        private int m_numCartridges;

        private void Awake() {
            m_numCartridges = m_initialCartridges.Length;

            GameMgr.Events.Register<LabInfo>(GameEvents.LabLoaded, HandleLabLoaded);
        }

        #region Handlers

        private void HandleLabLoaded(LabInfo labInfo) {
            GameObject newCartridge = Instantiate(m_cartridgePrefab.gameObject, this.transform);

            newCartridge.transform.localPosition += new Vector3(m_cartridgeSpacing * m_numCartridges * newCartridge.transform.localScale.x, 0, 0);

            newCartridge.GetComponent<Cartridge>()?.SetInfo(labInfo);

            // TODO: save the "home" position?

            m_numCartridges++;
        }

        #endregion // Handlers
    }
}
