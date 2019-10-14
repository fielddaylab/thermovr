using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class ThermoMath : MonoBehaviour
{
  //math limits ; xYz = vPt
  //Pa
  double p_min = IF97.get_Pmin()*1000000.0; // 611.213
  double p_max = IF97.get_Pmax()*1000000.0; // 100000000
  //Pa
  double psat_min = IF97.get_ptrip()*1000000.0; // ???
  double psat_max = IF97.get_pcrit()*1000000.0; // ???
  //M^3/kg
  double v_min = 1.0/3000;  //0.0003_
  double v_max = 1.0/0.001; //1000
  //K
  double t_min = IF97.get_Tmin(); // 273.15
  double t_max = IF97.get_Tmax(); // 1073.15

  int samples = 20;

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
  GameObject[] graph_bits;
  public Material graph_material;
  public GameObject pt_prefab;

  /*
  //IF97 API
  public static double rhomass_Tp(double T, double p)     // Get the mass density [kg/m^3] as a function of T [K] and p [Pa]
  public static double hmass_Tp(double T, double p)       // Get the mass enthalpy [J/kg] as a function of T [K] and p [Pa]
  public static double smass_Tp(double T, double p)       // Get the mass entropy [J/kg/K] as a function of T [K] and p [Pa]
  public static double umass_Tp(double T, double p)       // Get the mass internal energy [J/kg] as a function of T [K] and p [Pa]
  public static double cpmass_Tp(double T, double p)      // Get the mass constant-pressure specific heat [J/kg/K] as a function of T [K] and p [Pa]
  public static double cvmass_Tp(double T, double p)      // Get the mass constant-volume specific heat [J/kg/K] as a function of T [K] and p [Pa]
  public static double speed_sound_Tp(double T, double p) // Get the speed of sound [m/s] as a function of T [K] and p [Pa]
  public static double drhodp_Tp(double T, double p)      // Get the [d(rho)/d(p)]T [kg/mï¿½/Pa] as a function of T [K] and p [Pa]
  //IF97 Paper verified units:
  rhomass_Tp(T,p) | 1.0/rhomass_Tp(300, 3) = 0.00100215 | p:MPa, v:M^3/Kg, T:K | expects:K,MPa returns Kg/M^3
  */

  /*
  //IAPWS95 API
  public static double IAPWS95_pressure(double rho, double T);                                   //Input: rho in kg/m3, T in K, Output: Pa
  public static double IAPWS95_internal_energy(double rho, double T);                            //Input: rho in kg/m3, T in K, Output: Pa
  public static double IAPWS95_entropy(double rho, double T);                                    //Input: rho in kg/m3, T in K, Output: kJ/kg-K
  public static double IAPWS95_enthalpy(double rho, double T);                                   //Input: rho in kg/m3, T in K, Output: kJ/kg
  public static double IAPWS95_isochoric_heat_capacity(double rho, double T);                    //Input: rho in kg/m3, T in K, Output: kJ/kg-K
  public static double IAPWS95_isobaric_heat_capacity(double rho, double T);                     //Input: rho in kg/m3, T in K, Output: kJ/kg-K
  public static double IAPWS95_speed_of_sound(double rho, double T);                             //Input: rho in kg/m3, T in K, Output: m/s
  public static double IAPWS95_joule_thompson_coefficient(double rho, double T);                 //Input: rho in kg/m3, T in K
  public static double IAPWS95_isothermal_throttling_coefficient(double rho, double T);          //Input: rho in kg/m3, T in K
  public static double IAPWS95_isentropic_temperature_pressure_coefficent(double rho, double T); //Input: rho in kg/m3, T in K
  //IAPWS95 Paper verified units:
  IAPWS95_pressure(rho, T) | IAPWS95_pressure(999.887406, 275.0)/1000 = 0.0006982125 | p:MPa, v:Kg/M^3, t:K | expects:Kg/M^3,K returns KPa
  */

  // Start is called before the first frame update
  void Start()
  {
    sample_lbase_prev = sample_lbase;
    plot_lbase_prev = plot_lbase;

    IF97.initRegions();

    findObjects();
    genMesh();
    reset();
    derive();
    dotransform();
    //IF97.print_tables();
    //IAPWS95.print_tables();
    //compare_impls();
  }

  void compare_impls()
  {
    /*
    //passing test
    double t = 373.15; //K
    double v = 17.1969045; //M^3/Kg
    double p = IAPWS95.IAPWS95_pressure(1.0/v,t)*1000; //expects:Kg/M^3,K returns KPa
    Debug.LogFormat("{0,3:E} {1,3:E} {2,3:E}\n", p, v, t);
    v = 1.0/IF97.rhomass_Tp(t,p/1000000); //expects:K,MPa returns Kg/M^3
    Debug.LogFormat("{0,3:E} {1,3:E} {2,3:E}\n", p, v, t);
    */


    //*
    //IF97 primary
    int n_psamples = samples;
    int n_tsamples = samples;
    for(int y = 0; y < n_psamples; y++)
    {
      double pt = ((double)y/(n_psamples-1));
      for(int z = 0; z < n_tsamples; z++)
      {
        double tt = ((double)z/(n_tsamples-1));
        double pst = sample(pt);
        double tst = sample(tt);
        double p = Lerpd(p_min,p_max,pst);
        double t = Lerpd(t_min,t_max,tst);
        double v = 1.0/IF97.rhomass_Tp(t,p/1000000.0); //expects:K,MPa returns Kg/M^3
        //pvt in Pa, M^3/Kg, K
        double _p = IAPWS95.IAPWS95_pressure(1.0/v,t)*1000.0; //expects:Kg/M^3,K returns KPa

        Debug.LogFormat("error:{0} p:{1}Pa ({2}Pa), v:{3}M^3/Kg, t:{4}K",p-_p,p,_p,v,t);
      }
    }
    //*/

    /*
    //IAPWS95 primary
    int n_vsamples = samples;
    int n_tsamples = samples;
    for(int x = 0; x < n_vsamples; x++)
    {
      double vt = ((double)x/(n_vsamples-1));
      for(int z = 0; z < n_tsamples; z++)
      {
        double tt = ((double)z/(n_tsamples-1));
        double vst = sample(vt);
        double tst = sample(tt);
        double v = Lerpd(v_min,v_max,vst);
        double t = Lerpd(t_min,t_max,tst);
        double p = IAPWS95.IAPWS95_pressure(1.0/v,t)*1000.0; //expects:Kg/M^3,K returns KPa
        //pvt in Pa, M^3/Kg, K
        double _v = 1.0/IF97.rhomass_Tp(t,p/1000000.0); //expects:K,MPa returns Kg/M^3

        Debug.LogFormat("error:{0} p:{1}Pa, v:{2}M^3/Kg ({3}M^3/Kg), t:{4}K",v-_v,p,v,_v,t);
      }
    }
    //*/

  }

  double Lerpd(double a, double b, double t) { return (b-a)*t+a; }
  double Clampd(double v, double min, double max) { if(v < min) return min; if(v > max) return max; return v; } //v,min,max ordering mirrors Mathf.Clamp

  //sample bias- "graph density"
  [Range(0.001f,1000)]
  public double sample_lbase = 10.0f;
  double sample_lbase_prev = 0.0f;
  double sample(double t) { return Math.Pow(t,sample_lbase); }

  //plot bias- "graph zoom"
  [Range(0.001f,1000)]
  public double plot_lbase = 10.0f;
  double plot_lbase_prev = 0.0f;
  float log_plot(double min, double max, double val) { return (float)((Math.Log(val,plot_lbase)-Math.Log(min,plot_lbase))/(Math.Log(max,plot_lbase)-Math.Log(min,plot_lbase))); }

  float plot(double min, double max, double val) { return log_plot(min,max,val); }

  void genMesh()
  {
    var rand = new System.Random();
    int n_psamples = samples;
    int n_tsamples = samples;
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
        double pst = sample(pt);
        double tst = sample(tt);
        double p = Lerpd(p_min,p_max,pst);
        double t = Lerpd(t_min,t_max,tst);
        double v = 1.0/IF97.rhomass_Tp(t,p/1000000.0); //expects:K,MPa returns Kg/M^3
        //pvt in Pa, M^3/Kg, K

        //Debug.LogFormat("p:{0}Pa, v:{1}M^3/Kg, t:{2}K",p,v,t);
        float pplot = plot(p_min,p_max,p);
        float vplot = plot(v_min,v_max,v);
        float tplot = plot(t_min,t_max,t);

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
        double vst = sample(vt);
        double tst = sample(tt);
        double v = Lerpd(v_min,v_max,vst);
        double t = Lerpd(t_min,t_max,tst);
        double p = IAPWS95.IAPWS95_pressure(1.0/v,t)*1000.0; //expects:Kg/M^3,K returns KPa
        //pvt in Pa, M^3/Kg, K
        float pplot = plot(p_min,p_max,p);
        float vplot = plot(v_min,v_max,v);
        float tplot = plot(t_min,t_max,t);

        int i = n_tsamples*x+z;
        pt_positions[i] = new Vector3(vplot,pplot,tplot);
      }
    }
//*/

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
//*/

/*
    //HACK
    //gen assets
    graph_bits = new GameObject[n_groups*n_pts_per_group];
    for(int i = 0; i < n_groups*n_pts_per_group; i++)
    {
      graph_bits[i] = (GameObject)Instantiate(pt_prefab);
      graph_bits[i].transform.parent = graph.transform;
      graph_bits[i].transform.localPosition = pt_positions[i];
      graph_bits[i].transform.localScale = new Vector3(pt_size, pt_size, pt_size);
    }
//*/

//*
    //MESH
    List<Vector3> mesh_positions;
    List<Vector3> mesh_normals;
    List<int> mesh_triangles;

    mesh_positions = new List<Vector3>(pt_positions);
    mesh_normals = new List<Vector3>(n_pts);

    int vi = 0;
    int ni = 0;
    mesh_triangles = new List<int>((n_psamples-1)*(n_tsamples-1)*6);
    for(int y = 0; y < n_psamples-1; y++)
    {
      for(int z = 0; z < n_tsamples-1; z++)
      {
        vi = n_tsamples*y+z;
        mesh_triangles.Add(vi           +0); ni++;
        mesh_triangles.Add(vi+n_tsamples+0); ni++;
        mesh_triangles.Add(vi+n_tsamples+1); ni++;
        mesh_triangles.Add(vi           +0); ni++;
        mesh_triangles.Add(vi+n_tsamples+1); ni++;
        mesh_triangles.Add(vi           +1); ni++;

        mesh_normals.Add(Vector3.Cross(Vector3.Normalize(mesh_positions[vi+n_tsamples+0]-mesh_positions[vi+0]),Vector3.Normalize(mesh_positions[vi+1]-mesh_positions[vi+0])));
      }
      mesh_normals.Add(Vector3.Cross(Vector3.Normalize(mesh_positions[vi+n_tsamples+0]-mesh_positions[vi+1]),Vector3.Normalize(mesh_positions[vi+1]-mesh_positions[vi+n_tsamples+1])));
    }

    {
      //fill out the rest of the normals
      int y = n_psamples-1;
      for(int z = 0; z < n_tsamples-1; z++)
      {
        vi = n_tsamples*y+z;
        mesh_normals.Add(Vector3.Cross(Vector3.Normalize(mesh_positions[vi-n_tsamples+0]-mesh_positions[vi+0]),Vector3.Normalize(mesh_positions[vi+1]-mesh_positions[vi+0])));
      }
      mesh_normals.Add(Vector3.Cross(Vector3.Normalize(mesh_positions[vi-n_tsamples+0]-mesh_positions[vi+1]),Vector3.Normalize(mesh_positions[vi+1]-mesh_positions[vi-n_tsamples+1])));
    }

    for(var i = 0; i < mesh_triangles.Count; i+=3)
    {
      if(winding(mesh_positions[mesh_triangles[i+0]],mesh_positions[mesh_triangles[i+1]],mesh_positions[mesh_triangles[i+2]]) == -1)
        Debug.Log("GAH");
    }

    int concentrated_samples = samples*2;
    for(int i = 0; i < 2; i++)
    {
      for(int y = 0; y < concentrated_samples; y++)
      {
        double pt = ((double)y/(concentrated_samples-1));
        double pst = sample(pt);
        double p = Lerpd(psat_min,psat_max,pst);
        double v = 0.0;
             if(i == 0) v = 1.0/IF97.rholiq_p(p/1000000.0); //expects:MPa returns Kg/M^3
        else if(i == 1) v = 1.0/IF97.rhovap_p(p/1000000.0); //expects:MPa returns Kg/M^3
        double t = IF97.Tsat97(p/1000000.0);
        //pvt in Pa, M^3/Kg, K

        //Debug.LogFormat("p:{0}Pa, v:{1}M^3/Kg, t:{2}K",p,v,t);
        float pplot = plot(p_min,p_max,p);
        float vplot = plot(v_min,v_max,v);
        float tplot = plot(t_min,t_max,t);
        Vector3 point = new Vector3(vplot,pplot,tplot);
        merge_pt(point, mesh_positions, mesh_normals, mesh_triangles);
      }
    }

    Mesh mesh = new Mesh();
    mesh.vertices = mesh_positions.ToArray();
    mesh.normals = mesh_normals.ToArray();
    mesh.triangles = mesh_triangles.ToArray();

    GameObject gameObject = new GameObject("graph_mesh", typeof(MeshFilter), typeof(MeshRenderer));
    gameObject.transform.parent = graph.transform;
    gameObject.transform.localPosition = new Vector3(0.0f,0.0f,0.0f);
    gameObject.transform.localScale = new Vector3(1.0f,1.0f,1.0f);
    gameObject.GetComponent<MeshFilter>().mesh = mesh;
    gameObject.GetComponent<MeshRenderer>().material = graph_material;
//*/


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

  float cross_fv2z(Vector3 a, Vector3 b) { return (a.z*b.y)-(a.y*b.z); } //z of cross_fv3 with zs set to 0 (good for '2d pt in tri')
  float len_fv2(Vector3 v) { return Mathf.Sqrt(v.z*v.z+v.y*v.y); }
  float norm_f(float f) { if(f < 0.0f) return -1.0f; if(f > 0.0f) return 1.0f; return 0.0f; }

  bool pt_in_triangle_inclusive(Vector3 a, Vector3 b, Vector3 c, Vector3 p)
  {
    Vector3 subab = a-b;
    Vector3 subac = a-c;
    Vector3 subba = b-a;
    Vector3 subbc = b-c;
    Vector3 subca = c-a;
    Vector3 subcb = c-b;
    float one   = norm_f(cross_fv2z(subba,p-a));
    float two   = norm_f(cross_fv2z(subba,subca));
    float three = norm_f(cross_fv2z(subcb,p-b));
    float four  = norm_f(cross_fv2z(subcb,subab));
    float five  = norm_f(cross_fv2z(subac,p-c));
    float six   = norm_f(cross_fv2z(subac,subbc));

    return (
      (one == 0.0f || two == 0.0f || one == two)  &&
      (three == 0.0f || four == 0.0f || three == four)  &&
      (five == 0.0f || six == 0.0f || five == six)
    );
  }

  bool pt_in_triangle(Vector3 a, Vector3 b, Vector3 c, Vector3 p)
  {
    Vector3 subab = a-b;
    Vector3 subac = a-c;
    Vector3 subba = b-a;
    Vector3 subbc = b-c;
    Vector3 subca = c-a;
    Vector3 subcb = c-b;
    float one   = norm_f(cross_fv2z(subba,p-a));
    float two   = norm_f(cross_fv2z(subba,subca));
    float three = norm_f(cross_fv2z(subcb,p-b));
    float four  = norm_f(cross_fv2z(subcb,subab));
    float five  = norm_f(cross_fv2z(subac,p-c));
    float six   = norm_f(cross_fv2z(subac,subbc));

    return (
      one == two &&
      three == four &&
      five == six
    );
  }

  void delaunay_massage_zy(int i, List<Vector3> vbuff, List<int> ibuff)
  {
    //whoo boy.this is a bad (but straightforward...ish?) algorithm.
    //find and execute delaunay flips. if you flip, start over
    Vector3 circum_center = new Vector3();
    float circum_rad = 0.0f;
    int t1_ai = 0; Vector3 t1_a = new Vector3();
    int t1_bi = 0; Vector3 t1_b = new Vector3();
    int t1_ci = 0; Vector3 t1_c = new Vector3();
    int t2_ai = 0; Vector3 t2_a = new Vector3();
    int t2_bi = 0; Vector3 t2_b = new Vector3();
    int t2_ci = 0; Vector3 t2_c = new Vector3();
    int sai = 0;
    int sbi = 0;
    int sci = 0;
    int sdi = 0;

    //this is kinda clever (uh oh). ii starts at three triangle indexes before the end (max num created by adding vertex/splitting),
    //but if a flip is made it researches by setting ii = -3.
    for(int ii = i-9; ii < i; ii+=3)
    {
      t1_ai = ibuff[ii+0]; t1_a = vbuff[t1_ai];
      t1_bi = ibuff[ii+1]; t1_b = vbuff[t1_bi];
      t1_ci = ibuff[ii+2]; t1_c = vbuff[t1_ci];

      //don't *really* understand... from https://www.ics.uci.edu/~eppstein/junkyard/circumcenter.html
      float d = (t1_a.z - t1_c.z) * (t1_b.y - t1_c.y) - (t1_b.z - t1_c.z) * (t1_a.y - t1_c.y);
      circum_center.z = (((t1_a.z - t1_c.z) * (t1_a.z + t1_c.z) + (t1_a.y - t1_c.y) * (t1_a.y + t1_c.y)) / 2 * (t1_b.y - t1_c.y) - ((t1_b.z - t1_c.z) * (t1_b.z + t1_c.z) + (t1_b.y - t1_c.y) * (t1_b.y + t1_c.y)) / 2 * (t1_a.y - t1_c.y)) / d;
      circum_center.y = (((t1_b.z - t1_c.z) * (t1_b.z + t1_c.z) + (t1_b.y - t1_c.y) * (t1_b.y + t1_c.y)) / 2 * (t1_a.z - t1_c.z) - ((t1_a.z - t1_c.z) * (t1_a.z + t1_c.z) + (t1_a.y - t1_c.y) * (t1_a.y + t1_c.y)) / 2 * (t1_b.z - t1_c.z)) / d;
      circum_rad = Mathf.Sqrt(Mathf.Pow(t1_c.z-circum_center.z,2.0f) + Mathf.Pow(t1_c.y-circum_center.y,2.0f));

      int flip = 0;
      //does first prioritized iteration. at every ii, the entire ii list before has been internally exhausted.
      //add iii, test against all prev ii's, maintain internal exhaustion. (used in tandem with clever bit above- uh oh)
      for(int iii = 0; iii < ii && flip == 0; iii+=3)
      {
        int n_matching = 0;
        int t1a_matching = 0;
        int t1b_matching = 0;
        int t1c_matching = 0;
        int t2a_matching = 0;
        int t2b_matching = 0;
        int t2c_matching = 0;

        t2_ai = ibuff[iii+0]; t2_a = vbuff[t2_ai];
        t2_bi = ibuff[iii+1]; t2_b = vbuff[t2_bi];
        t2_ci = ibuff[iii+2]; t2_c = vbuff[t2_ci];

        if(t1_ai == t2_ai) { t1a_matching = 1; t2a_matching = 1; n_matching++; }
        if(t1_ai == t2_bi) { t1a_matching = 1; t2b_matching = 1; n_matching++; }
        if(t1_ai == t2_ci) { t1a_matching = 1; t2c_matching = 1; n_matching++; }
        if(t1_bi == t2_ai) { t1b_matching = 1; t2a_matching = 1; n_matching++; }
        if(t1_bi == t2_bi) { t1b_matching = 1; t2b_matching = 1; n_matching++; }
        if(t1_bi == t2_ci) { t1b_matching = 1; t2c_matching = 1; n_matching++; }
        if(t1_ci == t2_ai) { t1c_matching = 1; t2a_matching = 1; n_matching++; }
        if(t1_ci == t2_bi) { t1c_matching = 1; t2b_matching = 1; n_matching++; }
        if(t1_ci == t2_ci) { t1c_matching = 1; t2c_matching = 1; n_matching++; }

        if(n_matching == 2)
        {
          int unused_t1i = 0;
          if(t1a_matching == 0) unused_t1i = t1_ai;
          if(t1b_matching == 0) unused_t1i = t1_bi;
          if(t1c_matching == 0) unused_t1i = t1_ci;

               if(t2a_matching == 0 && (d != 0 && len_fv2(circum_center-t2_a)/circum_rad <= 0.5f))
          {
            flip = 1;
            sai = t2_ai;
            sbi = t2_bi;
            sci = t2_ci;
            sdi = unused_t1i;
          }
          else if(t2b_matching == 0 && (d != 0 && len_fv2(circum_center-t2_b)/circum_rad <= 0.5f))
          {
            flip = 1;
            sci = t2_ai;
            sai = t2_bi;
            sbi = t2_ci;
            sdi = unused_t1i;
          }
          else if(t2c_matching == 0 && (d != 0 && len_fv2(circum_center-t2_c)/circum_rad <= 0.5f))
          {
            flip = 1;
            sbi = t2_ai;
            sci = t2_bi;
            sai = t2_ci;
            sdi = unused_t1i;
          }

          if(flip != 0)
          {
            //  //cur layout, needs to be flipped
            //
            //  sdi sbi
            //     /
            //  sci sai

            ibuff[ii+0] = sdi;
            ibuff[ii+1] = sci;
            ibuff[ii+2] = sai;

            ibuff[iii+0] = sdi;
            ibuff[iii+1] = sai;
            ibuff[iii+2] = sbi;

            ii = -3; //reset search! (ouch)
          }
        }
      }
    }
  }
  void delaunay_massage_xz(int i, List<Vector3> vbuff, List<int> ibuff)
  {
    //whoo boy.this is a bad (but straightforward...ish?) algorithm.
    //find and execute delaunay flips. if you flip, start over
    Vector3 circum_center = new Vector3();
    float circum_rad = 0.0f;
    int t1_ai = 0; Vector3 t1_a = new Vector3();
    int t1_bi = 0; Vector3 t1_b = new Vector3();
    int t1_ci = 0; Vector3 t1_c = new Vector3();
    int t2_ai = 0; Vector3 t2_a = new Vector3();
    int t2_bi = 0; Vector3 t2_b = new Vector3();
    int t2_ci = 0; Vector3 t2_c = new Vector3();
    int sai = 0;
    int sbi = 0;
    int sci = 0;
    int sdi = 0;

    //this is kinda clever (uh oh). ii starts at three triangle indexes before the end (max num created by adding vertex/splitting),
    //but if a flip is made it researches by setting ii = -3.
    for(int ii = i-9; ii < i; ii+=3)
    {
      t1_ai = ibuff[ii+0]; t1_a = vbuff[t1_ai];
      t1_bi = ibuff[ii+1]; t1_b = vbuff[t1_bi];
      t1_ci = ibuff[ii+2]; t1_c = vbuff[t1_ci];

      //don't *really* understand... from https://www.ics.uci.edu/~eppstein/junkyard/circumcenter.html
      float d = (t1_a.x - t1_c.x) * (t1_b.z - t1_c.z) - (t1_b.x - t1_c.x) * (t1_a.z - t1_c.z);
      circum_center.x = (((t1_a.x - t1_c.x) * (t1_a.x + t1_c.x) + (t1_a.z - t1_c.z) * (t1_a.z + t1_c.z)) / 2 * (t1_b.z - t1_c.z) - ((t1_b.x - t1_c.x) * (t1_b.x + t1_c.x) + (t1_b.z - t1_c.z) * (t1_b.z + t1_c.z)) / 2 * (t1_a.z - t1_c.z)) / d;
      circum_center.z = (((t1_b.x - t1_c.x) * (t1_b.x + t1_c.x) + (t1_b.z - t1_c.z) * (t1_b.z + t1_c.z)) / 2 * (t1_a.x - t1_c.x) - ((t1_a.x - t1_c.x) * (t1_a.x + t1_c.x) + (t1_a.z - t1_c.z) * (t1_a.z + t1_c.z)) / 2 * (t1_b.x - t1_c.x)) / d;
      circum_rad = Mathf.Sqrt(Mathf.Pow(t1_c.x-circum_center.x,2.0f) + Mathf.Pow(t1_c.z-circum_center.z,2.0f));

      int flip = 0;
      //does first prioritized iteration. at every ii, the entire ii list before has been internally exhausted.
      //add iii, test against all prev ii's, maintain internal exhaustion. (used in tandem with clever bit above- uh oh)
      for(int iii = 0; iii < ii && flip == 0; iii+=3)
      {
        int n_matching = 0;
        int t1a_matching = 0;
        int t1b_matching = 0;
        int t1c_matching = 0;
        int t2a_matching = 0;
        int t2b_matching = 0;
        int t2c_matching = 0;

        t2_ai = ibuff[iii+0]; t2_a = vbuff[t2_ai];
        t2_bi = ibuff[iii+1]; t2_b = vbuff[t2_bi];
        t2_ci = ibuff[iii+2]; t2_c = vbuff[t2_ci];

        if(t1_ai == t2_ai) { t1a_matching = 1; t2a_matching = 1; n_matching++; }
        if(t1_ai == t2_bi) { t1a_matching = 1; t2b_matching = 1; n_matching++; }
        if(t1_ai == t2_ci) { t1a_matching = 1; t2c_matching = 1; n_matching++; }
        if(t1_bi == t2_ai) { t1b_matching = 1; t2a_matching = 1; n_matching++; }
        if(t1_bi == t2_bi) { t1b_matching = 1; t2b_matching = 1; n_matching++; }
        if(t1_bi == t2_ci) { t1b_matching = 1; t2c_matching = 1; n_matching++; }
        if(t1_ci == t2_ai) { t1c_matching = 1; t2a_matching = 1; n_matching++; }
        if(t1_ci == t2_bi) { t1c_matching = 1; t2b_matching = 1; n_matching++; }
        if(t1_ci == t2_ci) { t1c_matching = 1; t2c_matching = 1; n_matching++; }

        if(n_matching == 2)
        {
          int unused_t1i = 0;
          if(t1a_matching == 0) unused_t1i = t1_ai;
          if(t1b_matching == 0) unused_t1i = t1_bi;
          if(t1c_matching == 0) unused_t1i = t1_ci;

               if(t2a_matching == 0 && (d != 0 && len_fv2(circum_center-t2_a)/circum_rad <= 0.5f))
          {
            flip = 1;
            sai = t2_ai;
            sbi = t2_bi;
            sci = t2_ci;
            sdi = unused_t1i;
          }
          else if(t2b_matching == 0 && (d != 0 && len_fv2(circum_center-t2_b)/circum_rad <= 0.5f))
          {
            flip = 1;
            sci = t2_ai;
            sai = t2_bi;
            sbi = t2_ci;
            sdi = unused_t1i;
          }
          else if(t2c_matching == 0 && (d != 0 && len_fv2(circum_center-t2_c)/circum_rad <= 0.5f))
          {
            flip = 1;
            sbi = t2_ai;
            sci = t2_bi;
            sai = t2_ci;
            sdi = unused_t1i;
          }

          if(flip != 0)
          {
            //  //cur layout, needs to be flipped
            //
            //  sdi sbi
            //     /
            //  sci sai

            ibuff[ii+0] = sdi;
            ibuff[ii+1] = sci;
            ibuff[ii+2] = sai;

            ibuff[iii+0] = sdi;
            ibuff[iii+1] = sai;
            ibuff[iii+2] = sbi;

            ii = -3; //reset search! (ouch)
          }
        }
      }
    }
  }

  int winding(Vector3 a, Vector3 b, Vector3 c)
  {
      a.x = 0.0f;
      b.x = 0.0f;
      c.x = 0.0f;
      Vector3 cross = Vector3.Cross(a-c,b-c);
      if(cross.x > 0) return 1;
      if(cross.x < 0) return -1;
      return 0;
  }

  double amt_off_line(double ax, double ay, double bx, double by, double px, double py)
  {
    if(bx == ax) return Math.Abs(bx-px);
    if(by == ay) return Math.Abs(by-py);
    double xoff = (px-ax)/(bx-ax);
    double yoff = (py-ay)/(by-ay);
    return Math.Abs(xoff-yoff);
  }

  void merge_pt(Vector3 pt, List<Vector3> vbuff, List<Vector3> nbuff, List<int> ibuff)
  {
    for(int ti = 0; ti < vbuff.Count; ti++)
      if(pt == vbuff[ti]) return;

    int i = ibuff.Count;
    int start = vbuff.Count;

    //add new pt
    vbuff.Add(pt);

    //for every vert, find+split bounding triangle into three
    for(int v = start; v < vbuff.Count; v++)
    {

      bool found_tri = false;
      for(int ii = 0; ii < i && !found_tri; ii+=3)
      {
        int ai = ibuff[ii+0];
        int bi = ibuff[ii+1];
        int ci = ibuff[ii+2];
        int pi = v;

        Vector3 a = vbuff[ai];
        Vector3 b = vbuff[bi];
        Vector3 c = vbuff[ci];
        Vector3 p = vbuff[pi];

        if(pt_in_triangle(a,b,c,p))
        {
          found_tri = true;
          nbuff.Add(nbuff[ai]);

          //split bounding triangle into 3 by
          //swap cur triangle with last (to keep all new tris near end)
          ibuff[ii+0] = ibuff[i-3];
          ibuff[ii+1] = ibuff[i-2];
          ibuff[ii+2] = ibuff[i-1];
          if(winding(vbuff[ibuff[i-3]],vbuff[ibuff[i-2]],vbuff[ibuff[i-1]]) == -1) Debug.Log("Z");
          //then hijack existing triangle in ibuff
          ibuff[i-3] = ai;
          ibuff[i-2] = bi;
          ibuff[i-1] = pi;
          if(winding(vbuff[ai],vbuff[bi],vbuff[pi]) == -1) Debug.Log("A");
          //then adding two more triangle buffs
          ibuff.Add(bi); i++;
          ibuff.Add(ci); i++;
          ibuff.Add(pi); i++;
          if(winding(vbuff[bi],vbuff[ci],vbuff[pi]) == -1) Debug.Log("B");

          ibuff.Add(ci); i++;
          ibuff.Add(ai); i++;
          ibuff.Add(pi); i++;
          if(winding(vbuff[ci],vbuff[ai],vbuff[pi]) == -1) Debug.Log("C");

        }
        else if(pt_in_triangle_inclusive(a,b,c,p))
        {
          Debug.Log("inclusive!");
          found_tri = true;
          nbuff.Add(nbuff[ai]);

          double off_ab = amt_off_line(a.z,a.y,b.z,b.y,p.z,p.y);
          double off_bc = amt_off_line(b.z,b.y,c.z,c.y,p.z,p.y);
          double off_ca = amt_off_line(c.z,c.y,a.z,a.y,p.z,p.y);
          double min_off = 0.001f;
          int which = 0;
          if(off_ab < min_off) { min_off = off_ab; which = 1; }
          if(off_bc < min_off) { min_off = off_bc; which = 2; }
          if(off_ca < min_off) { min_off = off_ca; which = 3; }
          if(min_off == 0.0)
          {
                 if(off_ab == off_bc && off_ab == 0.0) which = 4; //pt ON b
            else if(off_bc == off_ca && off_ab == 0.0) which = 5; //pt ON c
            else if(off_ca == off_ab && off_ab == 0.0) which = 6; //pt ON a
          }

          if(which > 3) //ON a point
          {
            switch(which)
            {
              case 4: //on B
              {
                int tmp = ai;
                ai = ci;
                ci = bi;
                bi = tmp;
              }
              break;
              case 5: //on C
              {
                //good to go!
              }
              break;
              case 6: //on A
              {
                int tmp = ai;
                ai = bi;
                bi = ci;
                ci = tmp;
              }
              break;
            }
            //pt now on C

            if(vbuff[ci].x > pt.x)
              vbuff[ci] = new Vector3(vbuff[ci].x,vbuff[ci].y,vbuff[ci].z+0.00001f);
            else
              vbuff[ci] = new Vector3(vbuff[ci].x,vbuff[ci].y,vbuff[ci].z-0.00001f);

            ii -= 3;
            continue;
          }

          //BETWEEN points
          {
            switch(which)
            {
              case 0: //not on any particular line
              {
                Debug.Log("uh oh b");
              }
              break;
              case 1: //between a and b
              {
                int tmp = ai;
                ai = ci;
                ci = bi;
                bi = tmp;
              }
              break;
              case 2: //between b and c
              {
                //good to go!
              }
              break;
              case 3: //between c and a
              {
                int tmp = ai;
                ai = bi;
                bi = ci;
                ci = tmp;
              }
              break;
            }
            //position p now between b and c

            //find triangle that shares conflicting edge (if it exists)
            bool found_line = false;
            for(int iii = ii+3; iii < i && !found_line; iii+=3)
            {
              if(ibuff[iii+0] == ci && ibuff[iii+1] == bi)
              {
                found_line = true;
                //first add one triangle buff
                ibuff.Add(ibuff[iii+2]); i++;
                ibuff.Add(ci); i++;
                ibuff.Add(pi); i++;
                if(winding(vbuff[ibuff[iii+2]],vbuff[ci],vbuff[pi]) == -1)
                {
                  Debug.Log("(fix)A'");
                  //fix
                  ibuff[i-2] = pi;
                  ibuff[i-1] = ci;
                }

                //then split second in place
                ibuff[iii+0] = pi;
                ibuff[iii+1] = bi;
                ibuff[iii+2] = ibuff[iii+2]; //leave the same
                if(winding(vbuff[pi],vbuff[bi],vbuff[ibuff[iii+2]]) == -1) Debug.Log("B'");
              }
              else if(ibuff[iii+1] == ci && ibuff[iii+2] == bi)
              {
                found_line = true;
                //first add one triangle buff
                ibuff.Add(ibuff[iii+0]); i++;
                ibuff.Add(ci); i++;
                ibuff.Add(pi); i++;
                if(winding(vbuff[ibuff[iii+0]],vbuff[ci],vbuff[pi]) == -1)
                {
                  Debug.Log("(fix)A''");
                  //fix
                  ibuff[i-2] = pi;
                  ibuff[i-1] = ci;
                }

                //then split second in place
                ibuff[iii+0] = ibuff[iii+0]; //leave the same
                ibuff[iii+1] = pi;
                ibuff[iii+2] = bi;
                if(winding(vbuff[ibuff[iii+0]],vbuff[pi],vbuff[bi]) == -1) Debug.Log("B''");
              }
              else if(ibuff[iii+2] == ci && ibuff[iii+0] == bi)
              {
                found_line = true;
                //first add one triangle buff
                ibuff.Add(ibuff[iii+1]); i++;
                ibuff.Add(ci); i++;
                ibuff.Add(pi); i++;
                if(winding(vbuff[ibuff[iii+1]],vbuff[ci],vbuff[pi]) == -1) Debug.Log("A'''");

                //then split second in place
                ibuff[iii+0] = bi;
                ibuff[iii+1] = ibuff[iii+1]; //leave the same
                ibuff[iii+2] = pi;
                if(winding(vbuff[bi],vbuff[ibuff[iii+1]],vbuff[pi]) == -1) Debug.Log("B'''");
              }
              
            }
            if(!found_line) Debug.Log("UH OH");

            //split bounding triangle into 2 by
            //hijacking existing triangle in ibuff
            ibuff[ii+0] = ai;
            ibuff[ii+1] = bi;
            ibuff[ii+2] = pi;
            if(winding(vbuff[ai],vbuff[bi],vbuff[pi]) == -1)
            {
              Debug.Log("(fix)C'");
              //fix
              ibuff[ii+1] = pi;
              ibuff[ii+2] = bi;
            }

            //then adding one more triangle buffs
            ibuff.Add(ai); i++;
            ibuff.Add(pi); i++;
            ibuff.Add(ci); i++;
            if(winding(vbuff[ai],vbuff[pi],vbuff[ci]) == -1)
            {
              Debug.Log("(fix)D'");
              //fix
              ibuff[i-2] = ci;
              ibuff[i-1] = pi;
            }
          }
        }
      }

      delaunay_massage_zy(i, vbuff, ibuff);
      //delaunay_massage_xz(i, vbuff, ibuff);
    }
  }


}

