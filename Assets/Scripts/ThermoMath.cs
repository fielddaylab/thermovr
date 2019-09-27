using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class ThermoMath : MonoBehaviour
{
  //math limits ; xYz = vPt
  //Pa
  double p_min = IF97.get_Pmin()*1000000.0; // 0.000611213
  double p_max = IF97.get_Pmax()*1000000.0; // 100.0
  //Kg/M^3
  double v_min = 1.0/3000;
  double v_max = 1.0/0.001;
  //K
  double t_min = IF97.get_Tmin(); // 273.15
  double t_max = IF97.get_Tmax(); // 1073.15

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
  GameObject graph;
  GameObject[] axis;
  GameObject[,] axis_markers;
  GameObject[] plot_markers;
  GameObject[] graph_bits;
  public GameObject pt_prefab;

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
    IF97.initRegions();
    findObjects();
    genMesh();
    reset();
    derive();
    dotransform();
    //IF97.print_tables();

    sample_lbase_prev = sample_lbase;
    plot_lbase_prev = plot_lbase;
  }

  double Lerpd(double a, double b, double t) { return (b-a)*t+a; }
  double Clampd(double v, double min, double max) { if(v < min) return min; if(v > max) return max; return v; } //v,min,max ordering mirrors Mathf.Clamp

  //sample bias- "graph density"
  [Range(0.001f,1000)]
  public double sample_lbase = 10.0f;
  double sample_lbase_prev = 0.0f;
  double sampleP(double pt, double tt) { return Math.Pow(pt,sample_lbase); }
  double sampleV(double vt, double tt) { return Math.Pow(vt,sample_lbase); }
  double sampleT(double pt, double tt) { return Math.Pow(tt,sample_lbase); }

  //plot bias- "graph zoom"
  [Range(0.001f,1000)]
  public double plot_lbase = 10.0f;
  double plot_lbase_prev = 0.0f;
  float log_plot(double min, double max, double val) { return (float)((Math.Log(val,plot_lbase)-Math.Log(min,plot_lbase))/(Math.Log(max,plot_lbase)-Math.Log(min,plot_lbase))); }

  float p_plot(double min, double max, double val) { return log_plot(min,max,val); }
  float v_plot(double min, double max, double val) { return log_plot(min,max,val); }
  float t_plot(double min, double max, double val) { return log_plot(min,max,val); }

  void genMesh()
  {
    int n_psamples = 100;
    int n_vsamples = 100;
    int n_tsamples = 100;
    int n_pts = n_tsamples*n_psamples;
    int n_pts_per_group = 1000;
    int n_groups = (int)Mathf.Ceil(n_pts / n_pts_per_group);
    float pt_size = 0.005f;

    Vector3[] pt_positions;

//*IF97
    //gen positions
    pt_positions = new Vector3[n_pts];
    for(int y = 0; y < n_psamples; y++)
    {
      double pt = ((double)y/(n_psamples-1));
      for(int z = 0; z < n_tsamples; z++)
      {
        double tt = ((double)z/(n_tsamples-1));
        double pst = sampleP(pt,tt);
        double tst = sampleT(pt,tt);
        double p = Lerpd(p_min,p_max,pst);
        double t = Lerpd(t_min,t_max,tst);
        double v = 1.0/IF97.rhomass_Tp(t,p/1000000.0); //expects p:MPa, v:M^3/Kg, t:K

        //Debug.LogFormat("p:{0}Pa, v:{1}Kg/M^3, t:{2}K",p,v,t);
        float pplot = p_plot(p_min,p_max,p);
        float vplot = v_plot(v_min,v_max,v);
        float tplot = t_plot(t_min,t_max,t);

        int i = n_tsamples*y+z;

        pt_positions[i] = new Vector3(vplot,pplot,tplot);
      }
    }
//*/

/*IAPWS95
    //gen positions
    pt_positions = new Vector3[n_pts];
    for(int x = 0; x < n_vsamples; x++)
    {
      double vt = ((double)x/(n_vsamples-1));
      for(int z = 0; z < n_tsamples; z++)
      {
        double tt = ((double)z/(n_tsamples-1));
        double vst = sampleV(vt,tt);
        double tst = sampleT(vt,tt);
        double v = Lerpd(v_min,v_max,vst);
        double t = Lerpd(t_min,t_max,tst);
        double p = IAPWS95.IAPWS95_pressure(v,t);

        Debug.LogFormat("p:{0}Pa, v:{1}Kg/M^3, t:{2}K",p,v,t);
        float pplot = p_plot(p_min,p_max,p);
        float vplot = v_plot(v_min,v_max,v);
        float tplot = t_plot(t_min,t_max,t);

        int i = n_tsamples*x+z;

        //if(Double.IsNaN(vplot)) Debug.LogFormat("vi{0} {1} {2} {3}",i,vplot,pplot,tplot);
        //if(Double.IsNaN(pplot)) Debug.LogFormat("pi{0} {1} {2} {3}",i,vplot,pplot,tplot);
        //if(Double.IsNaN(tplot)) Debug.LogFormat("ti{0} {1} {2} {3}",i,vplot,pplot,tplot);

             if(i <  4200 &&  Double.IsNaN(pplot)) Debug.LogFormat("pi{0} {1} {2} {3}",i,vplot,pplot,tplot);
        else if(i >= 4200 && !Double.IsNaN(pplot)) Debug.LogFormat("pi{0} {1} {2} {3}",i,vplot,pplot,tplot);

        pt_positions[i] = new Vector3(vplot,pplot,tplot);
      }
    }
//*/





/*
    //x
    axis_markers[0,0].transform.position = new Vector3(p_plot(0,1,1.0/2),0,0);
    axis_markers[0,1].transform.position = new Vector3(p_plot(0,1,1.0/4),0,0);
    axis_markers[0,2].transform.position = new Vector3(p_plot(0,1,1.0/8),0,0);
    //y
    axis_markers[1,0].transform.position = new Vector3(0,v_plot(0,1,1.0/2),0);
    axis_markers[1,1].transform.position = new Vector3(0,v_plot(0,1,1.0/4),0);
    axis_markers[1,2].transform.position = new Vector3(0,v_plot(0,1,1.0/8),0);
    //z
    axis_markers[2,0].transform.position = new Vector3(0,0,t_plot(0,1,1.0/2));
    axis_markers[2,1].transform.position = new Vector3(0,0,t_plot(0,1,1.0/4));
    axis_markers[2,2].transform.position = new Vector3(0,0,t_plot(0,1,1.0/8));
*/

/*
    //gen assets
    graph_bits = new GameObject[n_groups];
    GameObject ptfab = (GameObject)Instantiate(pt_prefab);
    int n_pts_remaining = n_pts;
    int n_pts_this_group = n_pts_per_group;
    for(int i = 0; i < n_groups; i++)
    {
      n_pts_this_group = Mathf.Min(n_pts_per_group, n_pts_remaining);
      CombineInstance[] combine = new CombineInstance[n_pts_this_group];

      for(int j = 0; j < n_pts_this_group; j++)
      {
        ptfab.transform.position = pt_positions[i * n_pts_per_group + j];
        ptfab.transform.localScale = new Vector3(pt_size, pt_size, pt_size);

        combine[j].mesh = ptfab.transform.GetComponent<MeshFilter>().mesh;
        combine[j].transform = ptfab.transform.localToWorldMatrix;
      }

      graph_bits[i] = (GameObject)Instantiate(pt_prefab);
      graph_bits[i].transform.parent = graph.transform;
      graph_bits[i].transform.localPosition = new Vector3(0, 0, 0);
      graph_bits[i].transform.localRotation = Quaternion.Euler(0, 0, 0);
      graph_bits[i].transform.localScale = new Vector3(1, 1, 1);
      graph_bits[i].transform.GetComponent<MeshFilter>().mesh = new Mesh();
      graph_bits[i].transform.GetComponent<MeshFilter>().mesh.CombineMeshes(combine);

      n_pts_remaining -= n_pts_this_group;
    }
    Destroy(ptfab, 0f);
*/

    //HACK
    //gen assets
    graph_bits = new GameObject[n_groups*n_pts_per_group];
    for(int i = 0; i < n_groups*n_pts_per_group; i++)
    {
      graph_bits[i] = (GameObject)Instantiate(pt_prefab);
      graph_bits[i].transform.parent = graph.transform;
      /*
      if(Double.IsNaN(pt_positions[i].x)) Debug.LogFormat("xi{0}",i);
      if(Double.IsNaN(pt_positions[i].y)) Debug.LogFormat("yi{0}",i);
      if(Double.IsNaN(pt_positions[i].z)) Debug.LogFormat("zi{0}",i);
      */
      graph_bits[i].transform.position = pt_positions[i];
      graph_bits[i].transform.localScale = new Vector3(pt_size, pt_size, pt_size);
    }

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
    graph     = GameObject.Find("Graph");
    axis = new GameObject[3];
    axis[0]   = graph.transform.Find("x").gameObject;
    axis[1]   = graph.transform.Find("y").gameObject;
    axis[2]   = graph.transform.Find("z").gameObject;
    axis_markers = new GameObject[3,3];
    for(int i = 0; i < 3; i++)
    {
      axis_markers[i,0] = axis[i].transform.Find("mark_2").gameObject;
      axis_markers[i,1] = axis[i].transform.Find("mark_4").gameObject;
      axis_markers[i,2] = axis[i].transform.Find("mark_8").gameObject;
    }
    plot_markers = new GameObject[5];
    for(int i = 0; i < 5; i++)
      plot_markers[i] = graph.transform.Find("plots/"+i.ToString()).gameObject;
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
    bool modified = false;
    modified = ((plot_lbase != plot_lbase_prev) || (sample_lbase != sample_lbase_prev));
    sample_lbase_prev = sample_lbase;
    plot_lbase_prev = plot_lbase;
    if(modified)
    {
      //delete old
      for(int i = 0; i < graph_bits.Length; i++)
        Destroy(graph_bits[i]);
      genMesh();
    }
  }
}
