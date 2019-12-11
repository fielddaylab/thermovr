﻿using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using TMPro;

class CMP : IComparer<int>
{
  public List<Vector3> mesh_positions;
  public CMP(List<Vector3> _mesh_positions)
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
  //math limits ; xYz = vPt
  int samples = 350;

  //state
  public double pressure; //pascals
  public double temperature; //kalvin
  public double volume; //M^3/kg
  public double internalenergy; //J
  public double entropy; //?
  public double enthalpy; //?
  public double quality = 1; //%
  double prev_pressure;
  double prev_temperature;
  double prev_volume;
  double prev_internalenergy;
  double prev_entropy;
  double prev_enthalpy;
  double prev_quality;
  public double mass = 1; //kg
  public double radius = 0.05; //M
  public double surfacearea = 1.0; //Math.Pow(3.141592*radius,2.0);//M^2

  //vessel
  GameObject vessel;
  GameObject container;
  GameObject piston;
  float piston_min_y;
  float piston_max_y;
  GameObject contents;
  float contents_min_h;
  float contents_max_h;
  GameObject water;
  GameObject steam;

  //mesh
  GameObject graph;
  GameObject[] graph_bits;
  GameObject state;
  public Material graph_material;
  public GameObject pt_prefab;
  TextMeshPro text_pressure;
  TextMeshPro text_temperature;
  TextMeshPro text_volume;
  TextMeshPro text_entropy;
  TextMeshPro text_enthalpy;

  void Awake()
  {
    ThermoMath.Init();
  }

  // Start is called before the first frame update
  void Start()
  {
    sample_lbase_prev = sample_lbase;
    plot_lbase_prev = plot_lbase;

    findObjects();
    genMesh();
    //genHackMesh();

    reset();
    pressure = ThermoMath.p_given_percent(0.1);
    temperature = ThermoMath.t_given_percent(0.9);
    volume = ThermoMath.v_given_pt(pressure,temperature);
    /*
    Debug.LogFormat("{0}",pressure);
    Debug.LogFormat("{0}",volume);
    Debug.LogFormat("{0}",temperature);
    */
    dotransform();
  }

  //sample bias- "graph density"
  [Range(0.001f,20)]
  public double sample_lbase = 1.6f;
  double sample_lbase_prev = 0.0f;
  double sample(double t) { return Math.Pow(t,sample_lbase); }

  //plot bias- "graph zoom"
  [Range(0.001f,10)]
  public double plot_lbase = 10.0f;
  double plot_lbase_prev = 0.0f;
  float log_plot(double min, double max, double val) { double lval = Math.Log(val,plot_lbase); double lmax = Math.Log(max,plot_lbase); double lmin = Math.Log(min,plot_lbase); return (float)((lval-lmin)/(lmax-lmin)); }

  float plot(double min, double max, double val) { return log_plot(min,max,val); }
  String pv(Vector3 v) { return String.Format("{0}, {1}, {2}",v.x.ToString(".################"),v.y.ToString(".################"),v.z.ToString(".################")); }

  void genMesh()
  {
    int n_pts = samples*samples;
    int n_pts_per_group = 1000;
    int n_groups = (int)Mathf.Ceil(n_pts / n_pts_per_group);

    Vector3[] pt_positions;

    //gen positions
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
        double v = ThermoMath.v_given_pt(p,t);
        //pvt in Pa, M^3/Kg, K

        //Debug.LogFormat("p:{0}Pa, v:{1}M^3/Kg, t:{2}K",p,v,t);
        float pplot = plot(ThermoMath.p_min,ThermoMath.p_max,p);
        float vplot = plot(ThermoMath.v_min,ThermoMath.v_max,v);
        float tplot = plot(ThermoMath.t_min,ThermoMath.t_max,t);

        int i = samples*y+z;
        pt_positions[i] = new Vector3(vplot,pplot,tplot);
      }
    }

    //MESH
    List<Vector3> mesh_positions;
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
    int position_dome_region = mesh_positions.Count;
    float highest_y = 0.0f;
    int highest_y_i = 0;
    for(int y = 0; y < concentrated_samples; y++)
    {
      double pt = ((double)y/(concentrated_samples-1));
      double pst = sample(pt);
      double p = ThermoMath.psat_given_percent(pst);
      double t = ThermoMath.tsat_given_p(p);
      //pvt in Pa, M^3/Kg, K

      //Debug.LogFormat("p:{0}Pa, v:{1}M^3/Kg, t:{2}K",p,v,t);
      float pplot = plot(ThermoMath.p_min,ThermoMath.p_max,p);
      if(pplot > highest_y) { highest_y = pplot; highest_y_i = mesh_positions.Count; }
      float tplot = plot(ThermoMath.t_min,ThermoMath.t_max,t);

      double v;
      float vplot;
      Vector3 point;

      v = ThermoMath.vliq_given_p(p);
      vplot = plot(ThermoMath.v_min,ThermoMath.v_max,v);
      point = new Vector3(vplot,pplot,tplot);
      mesh_positions.Add(point);

      v = ThermoMath.vvap_given_p(p);
      vplot = plot(ThermoMath.v_min,ThermoMath.v_max,v);
      point = new Vector3(vplot,pplot,tplot);
      mesh_positions.Add(point);
    }
    highest_y = Mathf.Lerp(highest_y,1.0f,0.01f); //extra nudge up

    //kill spanning triangles; gather orphans
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

      float x_cmp = (left_ladder.x+right_ladder.x)/2.0f;
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
    CMP cmp = new CMP(mesh_positions);

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

    GameObject gameObject = new GameObject("graph_mesh", typeof(MeshFilter), typeof(MeshRenderer));
    gameObject.transform.parent = graph.transform;
    gameObject.transform.localPosition = new Vector3(0.0f,0.0f,0.0f);
    gameObject.transform.localScale = new Vector3(1.0f,1.0f,1.0f);
    gameObject.transform.localRotation = Quaternion.identity;
    gameObject.GetComponent<MeshFilter>().mesh = mesh;
    gameObject.GetComponent<MeshRenderer>().material = graph_material;
  }

  void genHackMesh()
  {
    int n_pts = samples*samples;
    int n_pts_per_group = 1000;
    int n_groups = (int)Mathf.Ceil(n_pts / n_pts_per_group);
    float pt_size = 0.005f;

    Vector3[] pt_positions;

//* derived v
    //gen positions
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
        double v = ThermoMath.v_given_pt(p,t);
        //pvt in Pa, M^3/Kg, K

        //Debug.LogFormat("p:{0}Pa, v:{1}M^3/Kg, t:{2}K",p,v,t);
        float pplot = plot(ThermoMath.p_min,ThermoMath.p_max,p);
        float vplot = plot(ThermoMath.v_min,ThermoMath.v_max,v);
        float tplot = plot(ThermoMath.t_min,ThermoMath.t_max,t);

        int i = samples*y+z;
        pt_positions[i] = new Vector3(vplot,pplot,tplot);
      }
    }
//*/

/* derived p
    //gen positions
    pt_positions = new Vector3[n_pts];
    for(int x = 0; x < samples; x++)
    {
      double vt = ((double)x/(samples-1));
      for(int z = 0; z < samples; z++)
      {
        double tt = ((double)z/(samples-1));
        double vst = sample(vt);
        double tst = sample(tt);
        double v = ThermoMath.v_given_percent(vst);
        double t = ThermoMath.t_given_percent(tst);
        double p = ThermoMath.p_given_vt(v,t);
        //pvt in Pa, M^3/Kg, K
        float pplot = plot(ThermoMath.p_min,ThermoMath.p_max,p);
        float vplot = plot(ThermoMath.v_min,ThermoMath.v_max,v);
        float tplot = plot(ThermoMath.t_min,ThermoMath.t_max,t);

        int i = samples*x+z;
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

//*
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
  }

  void reset()
  {
    //state
    pressure    = 0;
    temperature = 0;
    volume      = 0;
    entropy     = 0;
    enthalpy    = 0;
    prev_pressure    = -1;
    prev_temperature = -1;
    prev_volume      = -1;
    prev_entropy     = -1;
    prev_enthalpy    = -1;
  }

  void findObjects()
  {
    vessel    = GameObject.Find("Vessel");
    container = GameObject.Find("Container");
    piston    = GameObject.Find("Piston");
    piston_min_y = piston.transform.localPosition.y;
    piston_max_y = piston_min_y+0.1f; //experimentally derived...
    contents = GameObject.Find("Contents");
    contents_min_h = contents.transform.localScale.y;
    contents_max_h = contents_min_h+0.1f; //experimentally derived...
    water     = GameObject.Find("Water");
    steam     = GameObject.Find("Steam");

    graph     = GameObject.Find("gmodel");
    state     = GameObject.Find("gstate");

    text_pressure    = GameObject.Find("text_pressure").GetComponent<TextMeshPro>();
    text_temperature = GameObject.Find("text_temperature").GetComponent<TextMeshPro>();
    text_volume      = GameObject.Find("text_volume").GetComponent<TextMeshPro>();
    text_entropy     = GameObject.Find("text_entropy").GetComponent<TextMeshPro>();
    text_enthalpy    = GameObject.Find("text_enthalpy").GetComponent<TextMeshPro>();
  }

  /*
  //assume starting point consistent:
  pressure; //pascals
  temperature; //kalvin
  volume; //M^3/kg
  internalenergy; //J
  entropy; //?
  enthalpy; //?
  quality; //%
  */

  public void add_heat_constant_p(double j)
  {
    //newie = q - p(newv-oldv) + ie;
  /*
    int region = 0;
    switch(region)
    {
      case 0: //subcooled liquid
      {
        double q = j*t; //watts
        newinternalenergy = internalenergy+(q*t/mass) - pressure*(newvolume-volume);
        newvolume = get_volume_given(newinternalenergy,pressure);
        //at this point, we have enough internal state to derive the rest
      }
      break;
      case 1: //two-phase region
      {
        double q = j*t; //watts
        newvolume = get_volume_given(newinternalenergy,pressure);
        internalenergy = internalenergy+(q*t/mass) - pressure*(newvolume-volume);
        //at this point, we have enough internal state to derive the rest
      }
      break;
      case 2: //superheated vapor
      {
        double q = j*t; //watts
        newvolume = get_volume_given(newinternalenergy,pressure);
        internalenergy = internalenergy+(q*t/mass) - pressure*(newvolume-volume);
        //at this point, we have enough internal state to derive the rest
      }
      break;
    }
  */
    dotransform();
  }

  public void add_heat_constant_v(double j)
  {
    //newie = q - p(newv-oldv) + ie;
  /*
    int region = 0;
    switch(region)
    {
      case 0: //subcooled liquid
      {
        double q = j*t; //watts
        newinternalenergy = internalenergy+(q*t/mass) - pressure*(newvolume-volume);
        //at this point, we have enough internal state to derive the rest
      }
      break;
      case 1: //two-phase region
      {
        double q = j*t; //watts
        newinternalenergy = internalenergy+(q*t/mass);
        //at this point, we have enough internal state to derive the rest
      }
      break;
      case 2: //superheated vapor
      {
        double q = j*t; //watts
        newinternalenergy = internalenergy+(q*t/mass);
        //at this point, we have enough internal state to derive the rest
      }
      break;
    }
  */
    dotransform();
  }

  public void add_pressure(double w)
  {
  /*
    int region = 0;
    switch(region)
    {
      case 0: //subcooled liquid
      {
        double p = w*blah;
        newpressure = pressure+p;
        //at this point, we have enough internal state to derive the rest
      }
      break;
      case 1: //two-phase region
      {
        double p = w*blah;
        newpressure = pressure+p;
        //at this point, we have enough internal state to derive the rest
      }
      break;
      case 2: //superheated vapor
      {
        double p = w*blah;
        newpressure = pressure+p;
      }
      break;
    }
  */
    dotransform();
  }

  void dotransform()
  {
    float pplot = plot(ThermoMath.p_min,ThermoMath.p_max,pressure);
    float vplot = plot(ThermoMath.v_min,ThermoMath.v_max,volume);
    float tplot = plot(ThermoMath.t_min,ThermoMath.t_max,temperature);
    state.transform.localPosition = new Vector3(vplot,pplot,tplot);

    //ACTUALLY DO THIS CALCULATION
    float size_p = (float)ThermoMath.percent_given_v(volume);
    Vector3 piston_lt = piston.transform.localPosition;
    piston_lt.y = piston_min_y+size_p*(piston_max_y-piston_min_y);
    piston.transform.localPosition = piston_lt;

    Vector3 contents_lt = contents.transform.localScale;
    contents_lt.y = contents_min_h+size_p*(contents_max_h-contents_min_h);
    contents.transform.localScale = contents_lt;

    Vector3 water_lt = water.transform.localScale;
    water_lt.y = (float)quality;
    water.transform.localScale = water_lt;
    Vector3 steam_lt = steam.transform.localScale;
    steam_lt.y = 1.0f-(float)quality;
    steam.transform.localScale = steam_lt;
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
      /*
      //delete old
      for(int i = 0; i < graph_bits.Length; i++)
        Destroy(graph_bits[i]);
      genHackMesh();
      */
      Destroy(GameObject.Find("graph_mesh"));
      genMesh();
    }

    if(Math.Abs(pressure    - prev_pressure)    > 0.001) text_pressure.SetText(   "P: {0:3}KP",     (float)pressure/1000.0f);
    if(Math.Abs(temperature - prev_temperature) > 0.001) text_temperature.SetText("T: {0:3}K",      (float)temperature);
    if(Math.Abs(volume      - prev_volume)      > 0.001) text_volume.SetText(     "v: {0:3}M^3/kg", (float)volume);
    if(Math.Abs(entropy     - prev_entropy)     > 0.001) text_entropy.SetText(    "s: {0:3}J",      (float)entropy);
    if(Math.Abs(enthalpy    - prev_enthalpy)    > 0.001) text_enthalpy.SetText(   "h: {0:3}J",      (float)enthalpy);

    prev_pressure    = pressure;
    prev_temperature = temperature;
    prev_volume      = volume;
    prev_entropy     = entropy;
    prev_enthalpy    = enthalpy;
  }

}

