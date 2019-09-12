using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class ThermoMath : MonoBehaviour
{
  //state
  public double pressure_p;
  public double temperature_k;
  public double specificvolume_q;
  public double entropy;
  public double enthalpy;

  //constraints
  public double content_moles;
  public double radius_m;
  public double minstop_m;
  public double maxstop_m;
  public double weight_g;
  public double lift_g;
  public double flame_k;
  public double coolant_k;

  //derived
  public double pistonheight_m;
  public double contentvolume_m3;

  //vessel
  GameObject vessel;
  GameObject container;
  GameObject contents;
  GameObject piston;
  GameObject minstop;
  GameObject maxstop;
  GameObject flame;
  GameObject coolant;
  GameObject weights;
  GameObject lifts;

  // Start is called before the first frame update
  void Start()
  {
    reset();
    derive();
    findObjects();
    dotransform();
  }

  void reset()
  {
    //state
    pressure_p = 0;
    temperature_k = 0;
    specificvolume_q = 0;
    entropy = 0;
    enthalpy = 0;

    //constraints
    content_moles = 0;
    radius_m = 0;
    minstop_m = 0;
    maxstop_m = 0;
    weight_g = 0;
    lift_g = 0;
    flame_k = 0;
    coolant_k = 0;
  }

  void findObjects()
  {
  /*
    vessel    = GameObject.find("Vessel");
    container = GameObject.find("Container");
    contents  = GameObject.find("Contents");
    piston    = GameObject.find("Piston");
    minstop   = GameObject.find("Minstop");
    maxstop   = GameObject.find("Maxstop");
    flame     = GameObject.find("Flame");
    coolant   = GameObject.find("Coolant");
    weights   = GameObject.find("Weights");
    lifts     = GameObject.find("Lifts");
  */
  }

  void derive()
  {
    pistonheight_m = 0;
    contentvolume_m3 = 0;
  }

  void dotransform()
  {
    /*
    vessel = Scene.findObjects();
    container = Scene.findObjects();
    contents = Scene.findObjects();
    piston = Scene.findObjects();
    minstop = Scene.findObjects();
    maxstop = Scene.findObjects();
    flame = Scene.findObjects();
    coolant = Scene.findObjects();
    weights = Scene.findObjects();
    lifts = Scene.findObjects();
    */
  }

  // Update is called once per frame
  void Update()
  {

  }
}
