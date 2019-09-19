using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class ThermoMath : MonoBehaviour
{
  //math limits ; xYz = pVt
  double p_min = 0.0035;
  double p_max = 50;
  double v_min = 0;
  double v_max = 1000;
  double t_min = 300;
  double t_max = 2000;

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

  //mesh
  public GameObject pt_prefab;
  GameObject mesh_pts;

  /*
  //MATH API
  public static double rhomass_Tp(double T, double p)     // Get the mass density [kg/m^3] as a function of T [K] and p [Pa]
  public static double hmass_Tp(double T, double p)       // Get the mass enthalpy [J/kg] as a function of T [K] and p [Pa]
  public static double smass_Tp(double T, double p)       // Get the mass entropy [J/kg/K] as a function of T [K] and p [Pa]
  public static double umass_Tp(double T, double p)       // Get the mass internal energy [J/kg] as a function of T [K] and p [Pa]
  public static double cpmass_Tp(double T, double p)      // Get the mass constant-pressure specific heat [J/kg/K] as a function of T [K] and p [Pa]
  public static double cvmass_Tp(double T, double p)      // Get the mass constant-volume specific heat [J/kg/K] as a function of T [K] and p [Pa]
  public static double speed_sound_Tp(double T, double p) // Get the speed of sound [m/s] as a function of T [K] and p [Pa]
  public static double drhodp_Tp(double T, double p)      // Get the [d(rho)/d(p)]T [kg/mï¿½/Pa] as a function of T [K] and p [Pa]
  */

  // Start is called before the first frame update
  void Start()
  {
    genMesh();
    findObjects();
    reset();
    derive();
    dotransform();
    IF97.print_tables();
  }

  double Lerpd(double a, double b, double t) { return (b-a)*t+a; }
  float NormalizeMeshPt(double min, double max, double val) { return (float)((val-min)/(max-min)); }

  void genMesh()
  {
    mesh_pts = GameObject.Find("MeshPts");

    GameObject[] pt_groups;
    GameObject pt;
    Vector3[] pt_positions;
    Vector3 pt_pos;

    int n_tsamples = 100;
    int n_psamples = 100;
    int n_pts = n_tsamples*n_psamples;
    int n_pts_per_group = 1000;
    int n_groups = (int)Mathf.Ceil(n_pts / n_pts_per_group);
    pt_groups = new GameObject[n_groups];
    pt = (GameObject)Instantiate(pt_prefab);
    float pt_size = 0.01f;

    //gen positions
    pt_positions = new Vector3[n_pts];
    for(int x = 0; x < n_tsamples; x++)
    {
      for(int y = 0; y < n_psamples; y++)
      {
        double p = Lerpd(p_min,p_max,((double)y/(n_psamples-1)));
        double t = Lerpd(t_min,t_max,((double)x/(n_tsamples-1)));
        //Debug.LogFormat("t:{0} p:{1}",t,p);
        double v = IF97.rhomass_Tp(t,p);
        pt_positions[n_psamples*y+x] = new Vector3(NormalizeMeshPt(t_min,t_max,t),NormalizeMeshPt(v_min,v_max,v),NormalizeMeshPt(p_min,p_max,p));
      }
    }

    //gen assets
    int n_pts_remaining = n_pts;
    int n_pts_this_group = n_pts_per_group;
    for(int i = 0; i < n_groups; i++)
    {
      n_pts_this_group = Mathf.Min(n_pts_per_group, n_pts_remaining);
      CombineInstance[] combine = new CombineInstance[n_pts_this_group];

      for(int j = 0; j < n_pts_this_group; j++)
      {
        pt_pos = pt_positions[i * n_pts_per_group + j];
        pt.transform.position = pt_pos;
        pt.transform.localScale = new Vector3(pt_size, pt_size, pt_size);

        combine[j].mesh = pt.transform.GetComponent<MeshFilter>().mesh;
        combine[j].transform = pt.transform.localToWorldMatrix;
      }

      pt_groups[i] = (GameObject)Instantiate(pt_prefab);
      //pt_groups[i].transform.parent = ptsscale.transform;
      pt_groups[i].transform.localPosition = new Vector3(0, 0, 0);
      pt_groups[i].transform.localRotation = Quaternion.Euler(0, 0, 0);
      pt_groups[i].transform.localScale = new Vector3(1, 1, 1);
      pt_groups[i].transform.GetComponent<MeshFilter>().mesh = new Mesh();
      pt_groups[i].transform.GetComponent<MeshFilter>().mesh.CombineMeshes(combine);

      n_pts_remaining -= n_pts_this_group;
    }
    Destroy(pt, 0f);

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
    vessel    = GameObject.Find("Vessel");
    container = GameObject.Find("Container");
    contents  = GameObject.Find("Contents");
    piston    = GameObject.Find("Piston");
    minstop   = GameObject.Find("Minstop");
    maxstop   = GameObject.Find("Maxstop");
    flame     = GameObject.Find("Flame");
    coolant   = GameObject.Find("Coolant");
    weights   = GameObject.Find("Weights");
    lifts     = GameObject.Find("Lifts");
  }

  void derive()
  {
    pistonheight_m = 0;
    contentvolume_m3 = 0;
  }

  void dotransform()
  {

  }

  // Update is called once per frame
  void Update()
  {

  }
}
