using System.Collections;
using System.Collections.Generic;
using TMPro;
using UnityEngine;

public class GraphTextMgr : MonoBehaviour
{
    [SerializeField] private TMP_Text p_min_text;
    [SerializeField] private TMP_Text p_max_text;
    [SerializeField] private TMP_Text v_min_text;
    [SerializeField] private TMP_Text v_max_text;
    [SerializeField] private TMP_Text t_min_text;
    [SerializeField] private TMP_Text t_max_text;


    public void Start() {
        p_min_text.SetText(string.Format("{0:#.0e+0} " + Units.Pressure, ThermoMath.p_min));
        p_max_text.SetText(string.Format("{0:#.0e+0} " + Units.Pressure, ThermoMath.p_max));
        v_min_text.SetText(string.Format("{0:#.0e+0} " + Units.Volume, ThermoMath.v_min));
        v_max_text.SetText(string.Format("{0:#.0e+0} \n" + Units.Volume, ThermoMath.v_max));
        t_min_text.SetText(string.Format("{0:#.0e+0} " + Units.TemperatureK, ThermoMath.t_min));
        t_max_text.SetText(string.Format("{0:#.0e+0} " + Units.TemperatureK, ThermoMath.t_max));
    }
}
