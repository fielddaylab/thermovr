/*
DOCUMENTATION- phil, 12/16/19
This is the stateful wrapper around ThermoMath.
It represents a "current" state of water. After initialization, the state must always remain consistent.
For this reason, the API consists of applying deltas to some assumed consistent state.
It should be safe to assume that after any method call, the state remains consistent.

This is also responsible for applying itself visually to the game objects in the scene, (ie, position of the piston, % of water/steam, etc...) including generating the 3d graph
(^ this could be abstracted out into _another_ wrapper, but I don't see the value added [this is not a general-purpose thermomathstate class; it is built for this one VR application])
*/

using System;
using System.IO;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using TMPro;
using System.Collections.Specialized;

//One-Off class used for ordering points in graphgen zipper phase
class GRAPHPTCMP : IComparer<int>
{
  public List<Vector3> mesh_positions;
  public GRAPHPTCMP(List<Vector3> _mesh_positions)
  {
    mesh_positions = _mesh_positions;
  }

  public int Compare(int ai, int bi)
  {
    Vector3 a = mesh_positions[ai];
    Vector3 b = mesh_positions[bi];
    if(a.y > b.y) return 1;
    if(a.y < b.y) return -1;
    if(a.z > b.z) return 1;
    if(a.z < b.z) return -1;
    return 0;
  }
}

public class ThermoState : MonoBehaviour
{
  bool debug_write = false;
  StreamWriter debug_file;

  int samples = 350;

  //state
  //xyz corresponds to vpt (Y = "up")
  public double pressure;       //p //pascals
  public double temperature;    //t //°kelvin
  public double volume;         //v //M³/kg
  public double internalenergy; //u //J/kg
  public double entropy;        //s //J/kgK
  public double enthalpy;       //h //J/kg
  public double quality;        //x //%
  public int region;            //0 subcooled liquid, 1 two-phase, 2 superheated vapor
  double prev_pressure;
  double prev_temperature;
  double prev_volume;
  double prev_internalenergy;
  double prev_entropy;
  double prev_enthalpy;
  double prev_quality;
  int prev_region;

  //static properties of system
  public double mass = 1; //kg
  public double radius = 0.05; //M
  //public double surfacearea = Math.Pow(3.141592*radius,2.0); //M^2 //hardcoded answer below
  public double surfacearea = 0.024674011; //M^2 //hardcoded answer to eqn above
  public double surfacearea_insqr = 38.2447935395871; //in^2 //hardcoded conversion from m^2 to in^2

  //vessel
  GameObject vessel;
  GameObject container;
  GameObject piston;
  float piston_min_y;
  float piston_max_y;
  GameObject contents;
  float contents_min_h; //h = "height", not "enthalpy"
  float contents_max_h; //h = "height", not "enthalpy"
  GameObject water;
  GameObject steam;

  //mesh
  GameObject graph;
  GameObject state_dot;
  public Material graph_material;
  public Material graph_material_lit;
  TextMeshPro text_pressure;
  TextMeshPro text_temperature;
  TextMeshPro text_volume;
  TextMeshPro text_internalenergy;
  TextMeshPro text_entropy;
  TextMeshPro text_enthalpy;
  TextMeshPro text_quality;
  TextMeshPro text_region;

  public Flasher error_flasher;
  public TextMeshProUGUI error_message;

  public float size_p;

  void Awake()
  {
    ThermoMath.Init();
  }

  // Start is called before the first frame update
  void Start()
  {
    if(debug_write) debug_file = File.CreateText("debug.txt");

    //(these are just used to detect editor deltas on a frame boundary)
    sample_lbase_prev = sample_lbase;
    plot_lbase_prev = plot_lbase;

    findObjects();
    genMesh();
    reset_state();
    visualize_state();

    HideError();
  }

  public void Reset()
  {
    this.reset_state();
    clamp_state();
    visualize_state();
    this.HideError();
  }

  //sample bias- "graph density"
  [Range(0.001f,20)]
  public double sample_lbase = 1.6f;
  double sample_lbase_prev = 0f;
  double sample(double t) { return Math.Pow(t,sample_lbase); }

  //plot bias- "graph zoom"
  [Range(0.001f,10)]
  public double plot_lbase = 10f;
  double plot_lbase_prev = 0f;
  public float plot_dimension(double min, double max, double val) { double lval = Math.Log(val,plot_lbase); double lmax = Math.Log(max,plot_lbase); double lmin = Math.Log(min,plot_lbase); return (float)((lval-lmin)/(lmax-lmin)); }
  public float invplot_dimension(double min, double max, double res)
  {
    double lmax = Math.Log(max,plot_lbase);
    double lmin = Math.Log(min,plot_lbase);
    //return (float)Math.Pow(plot_lbase, 1.0/((res*(lmax-lmin))+lmin));
    return (float)Math.Pow(plot_lbase, (res*(lmax-lmin))+lmin);

/*
    double lval = ;
    double lmax = Math.Log(max,plot_lbase);
    double lmin = Math.Log(min,plot_lbase);
    (res*(lmax-lmin))+lmin = Math.Log(val,plot_lbase);
    plot_lbase^
*/

  }

  public Vector3 plot(double pressure, double volume, double temperature)
  {
    float pplot = plot_dimension(ThermoMath.p_min,ThermoMath.p_max,pressure);
    float vplot = plot_dimension(ThermoMath.v_min,ThermoMath.v_max,volume);
    float tplot = plot_dimension(ThermoMath.t_min,ThermoMath.t_max,temperature);
    return new Vector3(vplot,pplot,tplot);
  }
  public Vector3 invplot(double pplot, double vplot, double tplot)
  {
    float pressure    = invplot_dimension(ThermoMath.p_min,ThermoMath.p_max,pplot);
    float volume      = invplot_dimension(ThermoMath.v_min,ThermoMath.v_max,vplot);
    float temperature = invplot_dimension(ThermoMath.t_min,ThermoMath.t_max,tplot);
    return new Vector3(volume,pressure,temperature);
  }
  public void UpdateErrorState()
  {
    if (ThermoMath.got_error)
    {
      NotifyError();
    }
    else
    {
      HideError();
    }
  }

  private void NotifyError()
  {
    if (ThermoMath.got_error)
    {
      error_flasher.Flash();
      error_message.enabled = true;
    }
    else
    {
      Debug.Log("ThermoState was signaled to notify of a math state instability, but ThermoMath does not indicate an error occurred.");
    }
  }

  private void HideError()
  {
    ThermoMath.got_error = false;
    error_flasher.Stop();
    error_message.enabled = false;
  }

  //generates points from thermomath api, and stitches them together into a mesh
  //the "only reason" this is complex is:
  // we generate a "biased", "zoomed" grid of the mesh looked at from one axis ("looking at yz graph").
  // then we stitch this uniform (uniform other than bias/zoom, which can be "ignored") graph together.
  // however, there is a region of the generated graph ("the vapor dome") which is "constant z" (so invisible to this perspective).
  // so we detect triangles that span this "invisible" region, and cut them out of the stitching.
  // we then generate the vapor dome points _independently_, and create a very nice mesh of the region across the "xy" plane, which by design fits right into the cutaway stitching.
  // the final step then, is to "zip" together the two meshes.
  // this is done by walking the sorted list of "orphaned" points (<- I could have come up with a better name for that...), which corresponds to the list of points disconnected by the cutting of the grid mesh
  // and simultaneously walking the sorted list of the vapor dome region points, zig-zagging triangles to fill the space
  //the good news: any complexity from the generation of the mesh is pretty well isolated to this one function
  // NOTE: Phil didn't really say *what* this mesh is, just "a mesh". This function generates the graph object's mesh.

  //required for distance checks
  List<Vector3> mesh_positions;
  int position_dome_region;
  void genMesh()
  {
    GameObject old_gm = GameObject.Find("graph_mesh");
    if(old_gm != null) Destroy(old_gm);

    int n_pts = samples*samples;
    int n_pts_per_group = 1000;
    int n_groups = (int)Mathf.Ceil(n_pts / n_pts_per_group);

    //gen positions
    Vector3[] pt_positions;
    pt_positions = new Vector3[n_pts];
    for(int y = 0; y < samples; y++)
    {
      double pt = ((double)y/(samples-1));
      for(int z = 0; z < samples; z++)
      {
        double tt = ((double)z/(samples-1));
        double pst = sample(pt);
        double tst = sample(tt);
        double p = ThermoMath.p_given_percent(pst);
        double t = ThermoMath.t_given_percent(tst);
        double v = ThermoMath.v_given_pt(p,t, region);
        //pvt in Pa, M³/Kg, K

        //Debug.LogFormat("p:{0}Pa, v:{1}M³/Kg, t:{2}K",p,v,t);
        int i = samples*y+z;
        pt_positions[i] = plot(p,v,t);
      }
    }

    //MESH
    List<Vector3> mesh_normals;
    List<int> mesh_triangles;

    mesh_positions = new List<Vector3>(pt_positions);

    int vi = 0;
    int ni = 0;
    mesh_triangles = new List<int>((samples-1)*(samples-1)*6);
    for(int y = 0; y < samples-1; y++)
    {
      for(int z = 0; z < samples-1; z++)
      {
        vi = samples*y+z;
        mesh_triangles.Add(vi        +0); ni++;
        mesh_triangles.Add(vi+samples+0); ni++;
        mesh_triangles.Add(vi+samples+1); ni++;
        mesh_triangles.Add(vi        +0); ni++;
        mesh_triangles.Add(vi+samples+1); ni++;
        mesh_triangles.Add(vi        +1); ni++;
      }
    }

    int concentrated_samples = samples*2;
    position_dome_region = mesh_positions.Count;
    float highest_y = 0f;
    int highest_y_i = 0;
    for(int y = 0; y < concentrated_samples; y++)
    {
      double pt = ((double)y/(concentrated_samples-1));
      double pst = sample(pt);
      double p = ThermoMath.psat_given_percent(pst);
      double t = ThermoMath.tsat_given_p(p, region);
      //pvt in Pa, M³/Kg, K

      //Debug.LogFormat("p:{0}Pa, v:{1}M³/Kg, t:{2}° K",p,v,t);
      float pplot = plot_dimension(ThermoMath.p_min,ThermoMath.p_max,p);
      if(pplot > highest_y) { highest_y = pplot; highest_y_i = mesh_positions.Count; }
      float tplot = plot_dimension(ThermoMath.t_min,ThermoMath.t_max,t);

      double v;
      float vplot;
      Vector3 point;

      v = ThermoMath.vliq_given_p(p, region);
      vplot = plot_dimension(ThermoMath.v_min,ThermoMath.v_max,v);
      point = new Vector3(vplot,pplot,tplot);
      mesh_positions.Add(point);

      v = ThermoMath.vvap_given_p(p, region);
      vplot = plot_dimension(ThermoMath.v_min,ThermoMath.v_max,v);
      point = new Vector3(vplot,pplot,tplot);
      mesh_positions.Add(point);
    }
    highest_y = Mathf.Lerp(highest_y,1f,0.01f); //extra nudge up

    //kill spanning triangles; gather orphans
    //"ladder"/"rung" terminology a bit arbitrary- attempts to keep track of each side of a "zipper" for each seam ("left" seam, "right" seam, each have own ladder/rung)
    List<int> left_orphans = new List<int>();
    List<int> right_orphans = new List<int>();
    int left_ladder_i = position_dome_region;
    Vector3 left_ladder = mesh_positions[left_ladder_i];
    Vector3 left_rung = mesh_positions[left_ladder_i+2];
    int right_ladder_i = left_ladder_i+1;
    Vector3 right_ladder = mesh_positions[right_ladder_i];
    Vector3 right_rung = mesh_positions[right_ladder_i+2];
    for(var i = 0; i < mesh_triangles.Count; i+=3)
    {
      int ai = mesh_triangles[i+0];
      int bi = mesh_triangles[i+1];
      int ci = mesh_triangles[i+2];
      Vector3 a = mesh_positions[ai];
      Vector3 b = mesh_positions[bi];
      Vector3 c = mesh_positions[ci];

      if((left_rung.y  < a.y || left_rung.y  < b.y || left_rung.y  < c.y) && left_ladder_i+4  < mesh_positions.Count) { left_ladder_i  += 2; left_ladder  = mesh_positions[left_ladder_i];  left_rung  = mesh_positions[left_ladder_i+2];  }
      if((right_rung.y < a.y || right_rung.y < b.y || right_rung.y < c.y) && right_ladder_i+4 < mesh_positions.Count) { right_ladder_i += 2; right_ladder = mesh_positions[right_ladder_i]; right_rung = mesh_positions[right_ladder_i+2]; }

      float x_cmp = (left_ladder.x+right_ladder.x)/2f;
      if(
        (a.y < highest_y || b.y < highest_y || c.y < highest_y) &&
        (a.x < x_cmp || b.x < x_cmp || c.x < x_cmp) &&
        (a.x > x_cmp || b.x > x_cmp || c.x > x_cmp)
      )
      {
        mesh_triangles.RemoveAt(i+2);
        mesh_triangles.RemoveAt(i+1);
        mesh_triangles.RemoveAt(i+0);
        i -= 3;

        if(a.x < x_cmp && b.x < x_cmp)
        {
          left_orphans.Add(ai);
          left_orphans.Add(bi);
          right_orphans.Add(ci);
        }
        else if(b.x < x_cmp && c.x < x_cmp)
        {
          left_orphans.Add(bi);
          left_orphans.Add(ci);
          right_orphans.Add(ai);
        }
        else if(c.x < x_cmp && a.x < x_cmp)
        {
          left_orphans.Add(ci);
          left_orphans.Add(ai);
          right_orphans.Add(bi);
        }
        else if(a.x < x_cmp)
        {
          right_orphans.Add(bi);
          right_orphans.Add(ci);
          left_orphans.Add(ai);
        }
        else if(b.x < x_cmp)
        {
          right_orphans.Add(ai);
          right_orphans.Add(ci);
          left_orphans.Add(bi);
        }
        else if(c.x < x_cmp)
        {
          right_orphans.Add(ai);
          right_orphans.Add(bi);
          left_orphans.Add(ci);
        }
        else
        {
          Debug.Log("NOOOO");
        }
      }
    }

    //sort orphans
    GRAPHPTCMP cmp = new GRAPHPTCMP(mesh_positions);

    left_orphans.Sort(cmp);
    for(int i = 1; i < left_orphans.Count; i++)
    { if(left_orphans[i-1] == left_orphans[i]) { left_orphans.RemoveAt(i); i--; } }

    right_orphans.Sort(cmp);
    for(int i = 1; i < right_orphans.Count; i++)
    { if(right_orphans[i-1] == right_orphans[i]) { right_orphans.RemoveAt(i); i--; } }

    //stitch orphans
    int left_orphan_i = 0;
    int right_orphan_i = 0;
    {
      int triangle_stitch_region = mesh_triangles.Count;
      List<int> orphans;
      int ladder_i;
      Vector3 ladder;
      Vector3 rung;
      int orphan_i;
      Vector3 orphan;
      Vector3 orung;
      int ai = 0;
      int bi = 0;
      int ci = 0;

      //left
      orphans = left_orphans;
      orphan_i = 0;
      orphan = mesh_positions[orphans[orphan_i]];
      ladder_i = position_dome_region;
      ladder = mesh_positions[ladder_i];
      rung = mesh_positions[ladder_i+2];
      mesh_triangles.Add(ladder_i);
      mesh_triangles.Add(orphans[orphan_i]);
      orphan_i++;
      orphan = mesh_positions[orphans[orphan_i]];
      orung = mesh_positions[orphans[orphan_i+1]];
      mesh_triangles.Add(orphans[orphan_i]);
      orphan = mesh_positions[orphans[orphan_i]];
      while(ladder_i+2 < mesh_positions.Count)
      {
        while(orung.z <= rung.z && orung.y <= rung.y && orphan_i+1 < orphans.Count)
        { //increment orphan
          ai = ladder_i;
          bi = orphans[orphan_i];
          ci = orphans[orphan_i+1];
          mesh_triangles.Add(ai);
          mesh_triangles.Add(bi);
          mesh_triangles.Add(ci);

          orphan_i++;
          orphan = mesh_positions[orphans[orphan_i]];
          if(orphan_i+1 < orphans.Count) orung = mesh_positions[orphans[orphan_i+1]]; //yes, both this AND previous line need +1 (+1 for advance, +1 for orung)
        }
        if(ladder_i+2 < mesh_positions.Count)
        { //increment ladder
          ai = ladder_i;
          bi = orphans[orphan_i];
          ci = ladder_i+2;
          mesh_triangles.Add(ai);
          mesh_triangles.Add(bi);
          mesh_triangles.Add(ci);

          ladder_i += 2;
          ladder = mesh_positions[ladder_i];
          if(ladder_i+2 < mesh_positions.Count) rung = mesh_positions[ladder_i+2]; //yes, both this AND previous line need +2 (+2 for advance, +2 for rung)
        }
      }
      left_orphan_i = orphan_i;

      //right
      orphans = right_orphans;
      orphan_i = 0;
      orphan = mesh_positions[orphans[orphan_i]];
      orung = mesh_positions[orphans[orphan_i+1]];
      ladder_i = position_dome_region+1;
      ladder = mesh_positions[ladder_i];
      rung = mesh_positions[ladder_i+2];
      mesh_triangles.Add(orphans[orphan_i]);
      mesh_triangles.Add(ladder_i);
      ladder_i += 2;
      ladder = mesh_positions[ladder_i];
      mesh_triangles.Add(ladder_i);
      while(ladder_i+2 < mesh_positions.Count)
      {
        while((ladder.y > orung.y || rung.z > orung.z) && orphan_i+1 < orphans.Count)
        { //increment orphan
          ai = orphans[orphan_i];
          bi = ladder_i;
          ci = orphans[orphan_i+1];
          mesh_triangles.Add(ai);
          mesh_triangles.Add(bi);
          mesh_triangles.Add(ci);

          orphan_i++;
          orphan = mesh_positions[orphans[orphan_i]];
          if(orphan_i+1 < orphans.Count) orung = mesh_positions[orphans[orphan_i+1]]; //yes, both this AND previous line need +1 (+1 for advance, +1 for orung)
        }
        if(ladder_i+2 < mesh_positions.Count)
        { //increment ladder
          ai = orphans[orphan_i];
          bi = ladder_i;
          ci = ladder_i+2;
          mesh_triangles.Add(ai);
          mesh_triangles.Add(bi);
          mesh_triangles.Add(ci);

          ladder_i += 2;
          ladder = mesh_positions[ladder_i];
          if(ladder_i+2 < mesh_positions.Count) rung = mesh_positions[ladder_i+2]; //yes, both this AND previous line need +2 (+2 for advance, +2 for rung)
        }
      }
      right_orphan_i = orphan_i;
    }

    //fan missing top
    for(int i = left_orphan_i+1; i < left_orphans.Count; i++)
    {
      mesh_triangles.Add(left_orphans[i-1]);
      mesh_triangles.Add(left_orphans[i]);
      mesh_triangles.Add(highest_y_i);
    }
    for(int i = right_orphan_i+1; i < right_orphans.Count; i++)
    {
      mesh_triangles.Add(right_orphans[i]);
      mesh_triangles.Add(right_orphans[i-1]);
      mesh_triangles.Add(highest_y_i);
    }
    mesh_triangles.Add(left_orphans[left_orphans.Count-1]);
    mesh_triangles.Add(right_orphans[right_orphans.Count-1]);
    mesh_triangles.Add(highest_y_i);

    //fill in dome
    int triangle_inner_dome_region = mesh_triangles.Count;
    int position_dome_inner_region = mesh_positions.Count;
    for(int i = position_dome_region; i < position_dome_inner_region; i++) //duplicate inner positions so each can have own normal at seam
    {
      mesh_positions.Add(mesh_positions[i]);
    }
    for(int y = 0; y < concentrated_samples-1; y++)
    {
      mesh_triangles.Add(position_dome_inner_region+y*2+0);
      mesh_triangles.Add(position_dome_inner_region+y*2+2);
      mesh_triangles.Add(position_dome_inner_region+y*2+1);
      mesh_triangles.Add(position_dome_inner_region+y*2+1);
      mesh_triangles.Add(position_dome_inner_region+y*2+2);
      mesh_triangles.Add(position_dome_inner_region+y*2+3);
    }

    //set normals
    mesh_normals = new List<Vector3>(new Vector3[mesh_positions.Count]);
    for(int i = 0; i < triangle_inner_dome_region; i+=3)
    {
      int ai = mesh_triangles[i+0];
      int bi = mesh_triangles[i+1];
      int ci = mesh_triangles[i+2];
      Vector3 a = mesh_positions[ai];
      Vector3 b = mesh_positions[bi];
      Vector3 c = mesh_positions[ci];
      Vector3 n = Vector3.Cross(Vector3.Normalize(b-a),Vector3.Normalize(c-a));
      mesh_normals[ai] = n;
      mesh_normals[bi] = n;
      mesh_normals[ci] = n;
    }

    for(int i = triangle_inner_dome_region; i < mesh_triangles.Count; i+=3)
    {
      int ai = mesh_triangles[i+0];
      int bi = mesh_triangles[i+1];
      int ci = mesh_triangles[i+2];
      Vector3 a = mesh_positions[ai];
      Vector3 b = mesh_positions[bi];
      Vector3 c = mesh_positions[ci];
      Vector3 n = Vector3.Cross(Vector3.Normalize(b-a),Vector3.Normalize(c-a));
      mesh_normals[ai] = n;
      mesh_normals[bi] = n;
      mesh_normals[ci] = n;
    }

    Mesh mesh = new Mesh();
    mesh.indexFormat = UnityEngine.Rendering.IndexFormat.UInt32;
    mesh.vertices = mesh_positions.ToArray();
    mesh.normals = mesh_normals.ToArray();
    mesh.triangles = mesh_triangles.ToArray();

    GameObject graphObject = new GameObject("graph_mesh", typeof(MeshFilter), typeof(MeshRenderer), typeof(Lightable));
    graphObject.transform.parent = graph.transform;
    graphObject.transform.localPosition = new Vector3(0f,0f,0f);
    graphObject.transform.localScale = new Vector3(1f,1f,1f);
    graphObject.transform.localRotation = Quaternion.identity;
    graphObject.GetComponent<MeshFilter>().mesh = mesh;
    graphObject.GetComponent<MeshRenderer>().material = graph_material;
    graphObject.GetComponent<MeshRenderer>().shadowCastingMode = UnityEngine.Rendering.ShadowCastingMode.Off;
    var lightable = graphObject.GetComponent<Lightable>();
    lightable.use_custom_mats = true;
    lightable.base_mat = graph_material;
    lightable.lit_mat  = graph_material_lit;
  }

  void reset_state()
  {
    //ensure consistent state
    pressure       = ThermoMath.p_neutral[region];
    temperature    = ThermoMath.t_neutral[region];
    //from this point, the rest should be derived!
    volume         = ThermoMath.v_given_pt(pressure,temperature,region);
    internalenergy = ThermoMath.u_given_pt(pressure,temperature,region);
    enthalpy       = ThermoMath.h_given_vt(volume,temperature,region);
    entropy        = ThermoMath.s_given_vt(volume,temperature,region);
    quality        = ThermoMath.x_neutral[region];
    region         = ThermoMath.region_given_pvt(pressure,volume,temperature); //should certainly stay the same, as bases were calculated from assumed region

    prev_pressure       = -1;
    prev_temperature    = -1;
    prev_volume         = -1;
    prev_internalenergy = -1;
    prev_entropy        = -1;
    prev_enthalpy       = -1;
    prev_quality        = -1;
    prev_region         = -1;
  }

  double Clampd(double v, double min, double max) { if(v < min) return min;  if(v > max) return max; return v; }
  float Clampf(float v, float min, float max) { if(v < min) return min;  if(v > max) return max; return v; }
  void clamp_state()
  {
    if(Double.IsNaN(pressure))       pressure       = prev_pressure;
    if(Double.IsNaN(temperature))    temperature    = prev_temperature;
    if(Double.IsNaN(volume))         volume         = prev_volume;
    if(Double.IsNaN(internalenergy)) internalenergy = prev_internalenergy;
    if(Double.IsNaN(entropy))        entropy        = prev_entropy;
    if(Double.IsNaN(enthalpy))       enthalpy       = prev_enthalpy;
    if(Double.IsNaN(quality))        quality        = prev_quality;
    if(region == -1)                 region         = prev_region;

    //DEBUGGING! [COMMENT OUT IN PROD]
    /*
    double npressure       = Clampd(pressure,       ThermoMath.p_min,ThermoMath.p_max);
    double nvolume         = Clampd(volume,         ThermoMath.v_min,ThermoMath.v_max);
    double ntemperature    = Clampd(temperature,    ThermoMath.t_min,ThermoMath.t_max);
    double ninternalenergy = Clampd(internalenergy, ThermoMath.u_min,ThermoMath.u_max);
    double nentropy        = Clampd(entropy,        ThermoMath.s_min,ThermoMath.s_max);
    double nenthalpy       = Clampd(enthalpy,       ThermoMath.h_min,ThermoMath.h_max);
    double nquality        = Clampd(quality,        ThermoMath.x_min,ThermoMath.x_max);

    if(npressure       != pressure)       Debug.LogFormat("pressure!       {0} clamped to {1}",pressure,npressure);
    if(nvolume         != volume)         Debug.LogFormat("volume!         {0} clamped to {1}",volume,nvolume);
    if(ntemperature    != temperature)    Debug.LogFormat("temperature!    {0} clamped to {1}",temperature,ntemperature);
    if(ninternalenergy != internalenergy) Debug.LogFormat("internalenergy! {0} clamped to {1}",internalenergy,ninternalenergy);
    if(nentropy        != entropy)        Debug.LogFormat("entropy!        {0} clamped to {1}",entropy,nentropy);
    if(nenthalpy       != enthalpy)       Debug.LogFormat("enthalpy!       {0} clamped to {1}",enthalpy,nenthalpy);
    if(nquality        != quality)        Debug.LogFormat("quality!        {0} clamped to {1}",quality,nquality);
    */
    //END DEBUGGING

    //TODO: do this! should fix many errors!
    //snap to known region using region-stable perspectives
    if(prev_region != region)
    { //make sure new state snapped to new region!
      switch(region)
      {
        case 0:  //subcooled liquid
          //use P,T to fix the rest
          break;
        case 1:  //two-phase region
          //use P,V to fix the rest
          //temperature = ThermoMath.t_given_pv();
          break;
        case 2:  //superheated vapor
          //use P,T to fix the rest
          break;
      }
    }

    pressure       = Clampd(pressure,       ThermoMath.p_min,ThermoMath.p_max);
    volume         = Clampd(volume,         ThermoMath.v_min,ThermoMath.v_max);
    temperature    = Clampd(temperature,    ThermoMath.t_min,ThermoMath.t_max);
    internalenergy = Clampd(internalenergy, ThermoMath.u_min,ThermoMath.u_max);
    entropy        = Clampd(entropy,        ThermoMath.s_min,ThermoMath.s_max);
    enthalpy       = Clampd(enthalpy,       ThermoMath.h_min,ThermoMath.h_max);
    quality        = Clampd(quality,        ThermoMath.x_min,ThermoMath.x_max);
  }

  void findObjects()
  {
    vessel    = GameObject.Find("Vessel");
    container = GameObject.Find("Container");
    piston    = GameObject.Find("Piston");
    piston_min_y = piston.transform.localPosition.y;
    piston_max_y = piston_min_y+0.17f; //experimentally derived...
    contents = GameObject.Find("Contents");
    contents_min_h = contents.transform.localScale.y;
    contents_max_h = contents_min_h+0.17f; //experimentally derived...
    water     = GameObject.Find("Water");
    steam     = GameObject.Find("Steam");

    graph     = GameObject.Find("gmodel");
    state_dot = GameObject.Find("gstate");

    text_pressure       = GameObject.Find("text_pressure").GetComponent<TextMeshPro>();
    text_temperature    = GameObject.Find("text_temperature").GetComponent<TextMeshPro>();
    text_volume         = GameObject.Find("text_volume").GetComponent<TextMeshPro>();
    text_internalenergy = GameObject.Find("text_internalenergy").GetComponent<TextMeshPro>();
    text_entropy        = GameObject.Find("text_entropy").GetComponent<TextMeshPro>();
    text_enthalpy       = GameObject.Find("text_enthalpy").GetComponent<TextMeshPro>();
    text_quality        = GameObject.Find("text_quality").GetComponent<TextMeshPro>();
    text_region         = GameObject.Find("text_region").GetComponent<TextMeshPro>();
  }

  public void debug_deltas()
  {
    if(!debug_write) return;
    debug_file.WriteLine("pressure {0} changed to {1} (delta {2})",prev_pressure,pressure,pressure-prev_pressure);
    debug_file.WriteLine("temperature {0} changed to {1} (delta {2})",prev_temperature,temperature,temperature-prev_temperature);
    debug_file.WriteLine("volume {0} changed to {1} (delta {2})",prev_volume,volume,volume-prev_volume);
    debug_file.WriteLine("internalenergy {0} changed to {1} (delta {2})",prev_internalenergy,internalenergy,internalenergy-prev_internalenergy);
    debug_file.WriteLine("entropy {0} changed to {1} (delta {2})",prev_entropy,entropy,entropy-prev_entropy);
    debug_file.WriteLine("enthalpy {0} changed to {1} (delta {2})",prev_enthalpy,enthalpy,enthalpy-prev_enthalpy);
    debug_file.WriteLine("quality {0} changed to {1} (delta {2})",prev_quality,quality,quality-prev_quality);
  }

  //assume starting/ending point consistent for whole API!

  public Vector3 guessPlot(double guess_t, double pcube, double vcube)
  {
    //attempt to do it "the right way"
    /* //not worth the complexity
        Vector3 plot;
        plot.x = Clampf(invplot_dimension(ThermoMath.v_min, ThermoMath.v_max, vcube), (float)ThermoMath.v_min, (float)ThermoMath.v_max);
        plot.y = Clampf(invplot_dimension(ThermoMath.p_min, ThermoMath.p_max, pcube), (float)ThermoMath.p_min, (float)ThermoMath.p_max);
        plot.z = (float)ThermoMath.iterate_t_given_pv(guess_t, 10.0, (double)plot.y, (double)plot.x); //gives off warnings because now these math functions are stateful (bad)
        int region = ThermoMath.region_given_pvt(pressure,volume,temperature);
        if (region != 1)
        {
          plot.y = (float)ThermoMath.p_given_vt((double)plot.x, (double)plot.z);
          plot.x = (float)ThermoMath.v_given_pt((double)plot.y, (double)plot.z);
        }
        else
        {
           //TODO: !!
        }
        if (plot.z == 0.0) plot = new Vector3((float)volume, (float)pressure, (float)temperature);
        return plot;
    */
    return new Vector3(0, 0, 0);
  }

  public Vector3 guessMeshPlot(double vcube, double pcube, double tcube)
  {
    Vector3 pt = new Vector3((float)vcube, (float)pcube, (float)tcube);
    float d = Vector3.SqrMagnitude(pt - mesh_positions[0]);
    int closest = 0;
    for (int i = 1; i < mesh_positions.Count; i++)
    {
      float od = Vector3.SqrMagnitude(pt - mesh_positions[i]);
      if(od < d)
      {
        d = od;
        closest = i;
      }
    }

    Vector3 plot;
    if(closest < position_dome_region)
    {
      plot.x = Clampf(invplot_dimension(ThermoMath.v_min, ThermoMath.v_max, mesh_positions[closest].x), (float)ThermoMath.v_min, (float)ThermoMath.v_max);
      plot.y = Clampf(invplot_dimension(ThermoMath.p_min, ThermoMath.p_max, mesh_positions[closest].y), (float)ThermoMath.p_min, (float)ThermoMath.p_max);
      plot.z = Clampf(invplot_dimension(ThermoMath.t_min, ThermoMath.t_max, mesh_positions[closest].z), (float)ThermoMath.t_min, (float)ThermoMath.t_max);
    }
    else
    {
      plot.x = Clampf(invplot_dimension(ThermoMath.v_min, ThermoMath.v_max, vcube),                     (float)ThermoMath.v_min, (float)ThermoMath.v_max);
      plot.y = Clampf(invplot_dimension(ThermoMath.p_min, ThermoMath.p_max, mesh_positions[closest].y), (float)ThermoMath.p_min, (float)ThermoMath.p_max);
      plot.z = Clampf(invplot_dimension(ThermoMath.t_min, ThermoMath.t_max, mesh_positions[closest].z), (float)ThermoMath.t_min, (float)ThermoMath.t_max);
    }
    return plot;
  }

  public void warp_pv(double p, double v, double t)
  {
    try
    {
      pressure = p;
      volume = v;
      temperature = t;

      region = ThermoMath.region_given_pvt(pressure, volume, temperature);
      entropy = ThermoMath.s_given_vt(volume, temperature, region);
      enthalpy = ThermoMath.h_given_vt(volume, temperature, region);

      switch (region)
      {
        case 0: quality = 0; internalenergy = ThermoMath.u_given_vt(volume, temperature, region); break; //subcooled liquid
        case 1: quality = ThermoMath.x_given_pv(pressure, volume, region); internalenergy = ThermoMath.u_given_px(pressure, quality, region); break; //two-phase region
        case 2: quality = 1; internalenergy = ThermoMath.u_given_vt(volume, temperature, region); break; //superheated vapor
      }
    }
    catch (Exception e)
    {
      temperature = prev_temperature;
      pressure = prev_pressure;
      volume = prev_volume;
      region = prev_region;
      entropy = prev_entropy;
      enthalpy = prev_enthalpy;
      internalenergy = prev_internalenergy;
    }
    clamp_state();
    visualize_state();
  }

  public void add_heat_constant_p(double j)
  {
    try
    {
      double new_h = enthalpy+j;
      
      switch(region)
      {
        case 1:
          double new_x = ThermoMath.x_given_ph(pressure,new_h, region);
          //at this point, we have enough internal state to derive the rest
          enthalpy = new_h;
          quality = new_x;
          volume = ThermoMath.v_given_px(pressure,new_x, region);
          temperature = ThermoMath.tsat_given_p(pressure, region);
          entropy = ThermoMath.s_given_px(pressure,new_x, region);
          internalenergy = ThermoMath.u_given_px(pressure,new_x, region);
          break;
        case 0:
        case 2:
          //at this point, we have enough internal state to derive the rest
          enthalpy = new_h;
          volume = ThermoMath.v_given_ph(pressure, new_h, region);
          temperature = ThermoMath.t_given_ph(pressure, new_h, region);
          entropy = ThermoMath.s_given_vt(volume,temperature, region);
          internalenergy = ThermoMath.u_given_vt(volume, temperature, region);
          break;
      }

      region = ThermoMath.region_given_pvt(pressure,volume,temperature);
      switch(region)
      {
        case 0: quality = 0;                                       break; //subcooled liquid
        case 1: quality = ThermoMath.x_given_pv(pressure, volume, region); break; //two-phase region
        case 2: quality = 1;                                       break; //superheated vapor
      }
    }
    catch(Exception e) {}

    if(debug_write)
    {
      debug_file.WriteLine("add_heat_constant_p({0})",j);
      debug_deltas();
    }

    clamp_state();
    visualize_state();
  }

  public void add_heat_constant_v(double j)
  {
    try
    {
      double new_u = internalenergy+j;
      double new_t = ThermoMath.iterate_t_given_v_verify_u(temperature,volume,new_u, region);

      //at this point, we have enough internal state to derive the rest
      internalenergy = new_u;
      temperature = new_t;
      pressure = ThermoMath.p_given_vt(volume,temperature, region);
      enthalpy = ThermoMath.h_given_vt(volume,temperature, region);
      entropy = ThermoMath.s_given_vt(volume,temperature, region);

      region = ThermoMath.region_given_pvt(pressure,volume,temperature);
      switch(region)
      {
        case 0: quality = 0;                                       break; //subcooled liquid
        case 1: quality = ThermoMath.x_given_pv(pressure, volume, region); break; //two-phase region
        case 2: quality = 1;                                       break; //superheated vapor
      }
    }
    catch(Exception e) {}

    clamp_state();
    visualize_state();
  }

  public void add_pressure_uninsulated(double p)
  {
    try
    {
      double new_p = pressure+p;

      //default guess
      double new_u = internalenergy;
      double new_v = volume;

      //already done!
      new_v = ThermoMath.v_given_pt(new_p,temperature, region);
      new_u = ThermoMath.u_given_pt(new_p,temperature, region);
      //at this point, we have enough internal state to derive the rest
      pressure = new_p;
      volume = new_v;
      internalenergy = new_u;
      enthalpy = ThermoMath.h_given_vt(volume,temperature, region);
      entropy = ThermoMath.s_given_vt(volume,temperature, region);
      region = ThermoMath.region_given_pvt(pressure,volume,temperature);
    }
    catch(Exception e) {}

    clamp_state();
    visualize_state();
  }

  public void add_pressure_insulated(double p)
  {
    try
    {
      double new_p = pressure+p;

      switch(region)
      {
        case 0: //subcooled liquid
        case 1: //two-phase region
        {
          //AVOID THESE SCENARIOS
        }
        break;
        case 2: //superheated vapor
        {
          //default guess
          double new_t = temperature;
          double new_u = internalenergy;
          double new_v = volume;

          double k = 1.27;
          new_v = volume*Math.Pow(pressure/new_p,1.0/k);
          new_u = internalenergy-((new_p*new_v-pressure*volume)/(1-k));
          new_t = ThermoMath.iterate_t_given_p_verify_u(temperature,pressure,new_u, region);

          //at this point, we have enough internal state to derive the rest
          pressure = new_p;
          volume = new_v;
          temperature = new_t;
          internalenergy = new_u;
          enthalpy = ThermoMath.h_given_vt(volume,temperature, region);
          entropy = ThermoMath.s_given_vt(volume,temperature, region);
          region = ThermoMath.region_given_pvt(pressure,volume,temperature);
        }
        break;
      }
    }
    catch(Exception e) {}

    clamp_state();
    visualize_state();
  }

  void visualize_state()
  {
    state_dot.transform.localPosition = plot(pressure,volume,temperature);

    float height =  (float)(volume/surfacearea); //M
    float reductionFactor = 0.02f; //hack to reduce overall volume; do all calculations assuming 1.0kg of water, do all visualizations assuming reductionFactor*1.0kg of water. Assuming volume is linearly proportional to molarity (I _think_ it is?), we should be good. If not, we're still better than altering the behavior of vapor at an arbitrary threshhold
    height *= reductionFactor;
    float size_p = height/((float)radius*2f); //"max height" is approx 2x diameter, so this sets size_p to essentially "%_contents_size"
    Vector3 piston_lt = piston.transform.localPosition;
    piston_lt.y = piston_min_y+size_p*(piston_max_y-piston_min_y);
    piston.transform.localPosition = piston_lt;

    Vector3 contents_lt = contents.transform.localScale;
    contents_lt.y = contents_min_h+size_p*(contents_max_h-contents_min_h);
    contents.transform.localScale = contents_lt;

    Vector3 water_lt = water.transform.localScale;
    water_lt.y = 1f-(float)quality;
    water.transform.localScale = water_lt;
    //Vector3 steam_lt = steam.transform.localScale;
    //steam_lt.y = (float)quality;
    //steam.transform.localScale = -1f*steam_lt;
  }

  string regionToName(int region) //0 subcooled liquid, 1 two-phase, 2 superheated vapor
  {
    switch(region)
    {
      case 0: return "Subcooled Liquid";
      case 1: return "Two-Phase";
      case 2: return "Superheated Vapor";
    }
    return "Undefined";
  }

  // Update is called once per frame
  void Update()
  {
    //detect editor graphgen modifications
    bool modified = false;
    modified = ((plot_lbase != plot_lbase_prev) || (sample_lbase != sample_lbase_prev));
    sample_lbase_prev = sample_lbase;
    plot_lbase_prev = plot_lbase;
    if(modified) genMesh();

    if(Math.Abs(pressure - prev_pressure)              > ThermoMath.p_smallstep) text_pressure.SetText(      string.Format("P: {0:#.##E+0} kPa",   (float)pressure/1000f));
    if(Math.Abs(temperature - prev_temperature)        > ThermoMath.t_smallstep) text_temperature.SetText(                 "T: {0:3}°K ({1:3}°C)", (float)temperature, (float)temperature-273.15f);
    if(Math.Abs(volume - prev_volume)                  > ThermoMath.v_smallstep) text_volume.SetText(        string.Format("v: {0:#.##E+0} M³/kg", (float)volume));
    if(Math.Abs(internalenergy - prev_internalenergy)  > ThermoMath.u_smallstep) text_internalenergy.SetText(string.Format("u: {0:#.##E+0} kJ/kg", (float)internalenergy/1000f));
    if(Math.Abs(entropy - prev_entropy)                > ThermoMath.s_smallstep) text_entropy.SetText(                     "s: {0:3} kJ/kgK",      (float)entropy/1000f);
    if(Math.Abs(enthalpy - prev_enthalpy)              > ThermoMath.h_smallstep) text_enthalpy.SetText(      string.Format("h: {0:#.##E+0} kJ/kg", (float)enthalpy/1000f));
    if(region == 1 && Math.Abs(quality - prev_quality) > ThermoMath.x_smallstep) text_quality.SetText(                     "x: {0:3}%", (float)(quality * 100f));
    if(region != prev_region)
    {
      text_region.SetText("region: " + regionToName(region));
      if(region == 1)                                                            text_quality.SetText(                     "x: {0:3}%", (float)(quality * 100f));
      else                                                                       text_quality.SetText(                     "x: Undefined");
    }

    prev_pressure       = pressure;
    prev_temperature    = temperature;
    prev_volume         = volume;
    prev_internalenergy = internalenergy;
    prev_entropy        = entropy;
    prev_enthalpy       = enthalpy;
    prev_quality        = quality;
  }

}

