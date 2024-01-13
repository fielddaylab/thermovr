/*
DOCUMENTATION- phil
This is responsible for applying a ThermoSTate visually to the game objects in the scene, (ie, position of the piston, % of water/steam, etc...) including generating the 3d graph
*/

using System;
using System.IO;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using TMPro;
using System.Collections.Specialized;
using ThermoVR.Tools;
using ThermoVR;
using ThermoVR.State;

//One-Off class used for ordering points in graphgen zipper phase
class GRAPHPTCMP : IComparer<int>
{
    public List<Vector3> mesh_positions;
    public GRAPHPTCMP(List<Vector3> _mesh_positions) {
        mesh_positions = _mesh_positions;
    }

    public int Compare(int ai, int bi) {
        Vector3 a = mesh_positions[ai];
        Vector3 b = mesh_positions[bi];
        if (a.y > b.y) return 1;
        if (a.y < b.y) return -1;
        if (a.z > b.z) return 1;
        if (a.z < b.z) return -1;
        return 0;
    }
}

[RequireComponent(typeof(ThermoState))]
public class ThermoPresent : MonoBehaviour
{
    public static float max_height_log = 15; // log of max real-world height given 1 kg water and radius 0.15 m or 0.05 m

    bool debug_write = false;
    StreamWriter debug_file;

    int samples = 350;

    [SerializeField] ThermoState state;

    //vessel
    // GameObject vessel;
    // GameObject container;
    [SerializeField] private PistonController piston;
    [SerializeField] private ContentsController contents;

    //mesh
    GameObject graph;
    GameObject state_dot;
    [SerializeField] private GameObject dot_x_tracker;
    [SerializeField] private GameObject dot_y_tracker;
    [SerializeField] private GameObject dot_z_tracker;
    public Material graph_material;
    public Material graph_material_lit;

    public Flasher error_flasher;
    public TextMeshProUGUI error_message;

    public float size_p;

    private double pressure_range;
    private double temperature_range;
    private double volume_range;
    private double internal_energy_range;
    private double enthalpy_range;
    private double entropy_range;
    private double quality_range;

    void Awake() {
        ThermoMath.Init();
        state = this.GetComponent<ThermoState>();
        state.reset();
    }

    // Start is called before the first frame update
    public void Init() {
        if (debug_write) debug_file = File.CreateText("debug.txt");

        //(these are just used to detect editor deltas on a frame boundary)
        sample_lbase_prev = sample_lbase;
        plot_lbase_prev = plot_lbase;

        pressure_range = ThermoMath.p_max - ThermoMath.p_min;
        temperature_range = ThermoMath.t_max - ThermoMath.t_min;
        volume_range = ThermoMath.v_max - ThermoMath.v_min;
        internal_energy_range = 3670000f;
        enthalpy_range = 4130000f;
        entropy_range = 11500f;
        quality_range = 1;

        findObjects();
        genMesh();
        state.reset();
        visualize_state();

        HideError();
    }

    public void Reset() {
        state.reset();
        visualize_state();
        HideError();
    }

    //sample bias- "graph density"
    [Range(0.001f, 20)]
    public double sample_lbase = 1.6f;
    double sample_lbase_prev = 0f;
    double sample(double t) { return Math.Pow(t, sample_lbase); }

    //plot bias- "graph zoom"
    [Range(0.001f, 10)]
    public double plot_lbase = 10f;
    double plot_lbase_prev = 0f;
    public float plot_dimension(double min, double max, double val) { double lval = Math.Log(val, plot_lbase); double lmax = Math.Log(max, plot_lbase); double lmin = Math.Log(min, plot_lbase); return (float)((lval - lmin) / (lmax - lmin)); }
    public float invplot_dimension(double min, double max, double res) {
        double lmax = Math.Log(max, plot_lbase);
        double lmin = Math.Log(min, plot_lbase);
        //return (float)Math.Pow(plot_lbase, 1.0/((res*(lmax-lmin))+lmin));
        return (float)Math.Pow(plot_lbase, (res * (lmax - lmin)) + lmin);

        /*
            double lval = ;
            double lmax = Math.Log(max,plot_lbase);
            double lmin = Math.Log(min,plot_lbase);
            (res*(lmax-lmin))+lmin = Math.Log(val,plot_lbase);
            plot_lbase^
        */

    }

    public Vector3 plot(double pressure, double volume, double temperature) {
        float pplot = plot_dimension(ThermoMath.p_min, ThermoMath.p_max, pressure);
        float vplot = plot_dimension(ThermoMath.v_min, ThermoMath.v_max, volume);
        float tplot = plot_dimension(ThermoMath.t_min, ThermoMath.t_max, temperature);
        return new Vector3(vplot, pplot, tplot);
    }
    public Vector3 invplot(double pplot, double vplot, double tplot) {
        float pressure = invplot_dimension(ThermoMath.p_min, ThermoMath.p_max, pplot);
        float volume = invplot_dimension(ThermoMath.v_min, ThermoMath.v_max, vplot);
        float temperature = invplot_dimension(ThermoMath.t_min, ThermoMath.t_max, tplot);
        return new Vector3(volume, pressure, temperature);
    }
    public void UpdateErrorState() {
        if (ThermoMath.got_error) {
            NotifyError();
        }
        else {
            HideError();
        }
    }

    private void NotifyError() {
        if (ThermoMath.got_error) {
            error_flasher.Flash();
            error_message.enabled = true;
        }
        else {
            Debug.Log("ThermoState was signaled to notify of a math state instability, but ThermoMath does not indicate an error occurred.");
        }
    }

    private void HideError() {
        ThermoMath.got_error = false;
        error_flasher.Stop();
        if (error_message)
        {
            error_message.enabled = false;
        }
    }

    //generates points from thermomath api, and stitches them together into a mesh (the "graph")
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

    //required for distance checks
    List<Vector3> mesh_positions;
    int position_dome_region;

    GameObject genTestMeshFromPts(Vector3[] pt_positions, string name) {
        List<Vector3> mesh_positions = new List<Vector3>(pt_positions);
        List<int> mesh_triangles = new List<int>((samples - 1) * (samples - 1) * 6);

        int vi = 0;
        int ni = 0;
        for (int y = 0; y < samples - 1; y++) {
            for (int z = 0; z < samples - 1; z++) {
                vi = samples * y + z;
                mesh_triangles.Add(vi + 0); ni++;
                mesh_triangles.Add(vi + samples + 0); ni++;
                mesh_triangles.Add(vi + samples + 1); ni++;
                mesh_triangles.Add(vi + 0); ni++;
                mesh_triangles.Add(vi + samples + 1); ni++;
                mesh_triangles.Add(vi + 1); ni++;
            }
        }

        List<Vector3> mesh_normals = new List<Vector3>(new Vector3[mesh_positions.Count]);
        for (int i = 0; i < mesh_triangles.Count; i += 3) {
            int ai = mesh_triangles[i + 0];
            int bi = mesh_triangles[i + 1];
            int ci = mesh_triangles[i + 2];
            Vector3 a = mesh_positions[ai];
            Vector3 b = mesh_positions[bi];
            Vector3 c = mesh_positions[ci];
            Vector3 n = Vector3.Cross(Vector3.Normalize(b - a), Vector3.Normalize(c - a));
            mesh_normals[ai] = n;
            mesh_normals[bi] = n;
            mesh_normals[ci] = n;
        }

        Mesh mesh = new Mesh();
        mesh.indexFormat = UnityEngine.Rendering.IndexFormat.UInt32;
        mesh.vertices = mesh_positions.ToArray();
        mesh.normals = mesh_normals.ToArray();
        mesh.triangles = mesh_triangles.ToArray();

        GameObject graphObject = new GameObject(name, typeof(MeshFilter), typeof(MeshRenderer));
        graphObject.transform.parent = graph.transform;
        graphObject.transform.localPosition = new Vector3(0f, 0f, 0f);
        graphObject.transform.localScale = new Vector3(1f, 1f, 1f);
        graphObject.transform.localRotation = Quaternion.identity;
        graphObject.GetComponent<MeshFilter>().mesh = mesh;
        graphObject.GetComponent<MeshRenderer>().material = graph_material;
        graphObject.GetComponent<MeshRenderer>().shadowCastingMode = UnityEngine.Rendering.ShadowCastingMode.Off;

        return graphObject;
    }

    public double get_pressure() {
        return state.pressure;
    }

    public double get_temperature() {
        return state.temperature;
    }

    public double get_quality() {
        return state.quality;
    }

    public int get_region() {
        return state.region;
    }

    /// <summary>
    ///  Used for Reach State lab questions
    /// </summary>
    /// <param name="id"></param>
    /// <returns></returns>
    public double get_state_var(VarID id) {
        switch (id) {
            case VarID.Region:
                return state.region;
            case VarID.Pressure:
                return state.pressure / 1000; // in kPa
            case VarID.Temperature:
                return state.temperature;
            case VarID.Volume:
                return state.volume;
            case VarID.InternalEnergy:
                return state.internalenergy;
            case VarID.Entropy:
                return state.entropy;
            case VarID.Enthalpy:
                return state.enthalpy;
            case VarID.Quality:
                return state.quality;
            default:
                // VolumeStop
                return -1;
        }
    }

    void genMesh() {
        GameObject old_gm = GameObject.Find("graph_mesh");
        if (old_gm != null) Destroy(old_gm);

        int n_pts = samples * samples;
        int n_pts_per_group = 1000;
        int n_groups = (int)Mathf.Ceil(n_pts / n_pts_per_group);

        //gen positions
        Vector3[] pt_positions;
        pt_positions = new Vector3[n_pts];
        for (int y = 0; y < samples; y++) {
            double pt = ((double)y / (samples - 1));
            for (int z = 0; z < samples; z++) {
                double tt = ((double)z / (samples - 1));
                double pst = sample(pt);
                double tst = sample(tt);
                double p = ThermoMath.p_given_percent(pst);
                double t = ThermoMath.t_given_percent(tst);
                double v = ThermoMath.v_given_pt(p, t, state.region);
                //pvt in Pa, M³/Kg, K

                //Debug.LogFormat("p:{0}Pa, v:{1}M³/Kg, t:{2}K",p,v,t);
                int i = samples * y + z;
                pt_positions[i] = plot(p, v, t);
            }
        }

        //MESH
        List<Vector3> mesh_normals;
        List<int> mesh_triangles;

        mesh_positions = new List<Vector3>(pt_positions);

        int vi = 0;
        int ni = 0;
        mesh_triangles = new List<int>((samples - 1) * (samples - 1) * 6);
        for (int y = 0; y < samples - 1; y++) {
            for (int z = 0; z < samples - 1; z++) {
                vi = samples * y + z;
                mesh_triangles.Add(vi + 0); ni++;
                mesh_triangles.Add(vi + samples + 0); ni++;
                mesh_triangles.Add(vi + samples + 1); ni++;
                mesh_triangles.Add(vi + 0); ni++;
                mesh_triangles.Add(vi + samples + 1); ni++;
                mesh_triangles.Add(vi + 1); ni++;
            }
        }

        int concentrated_samples = samples * 2;
        position_dome_region = mesh_positions.Count;
        float highest_y = 0f;
        int highest_y_i = 0;
        for (int y = 0; y < concentrated_samples; y++) {
            double pt = ((double)y / (concentrated_samples - 1));
            double pst = sample(pt);
            double p = ThermoMath.psat_given_percent(pst);
            double t = ThermoMath.tsat_given_p(p, state.region);
            //pvt in Pa, m³/Kg, K

            //Debug.LogFormat("p:{0}Pa, v:{1}m³/Kg, t:{2}° K",p,v,t);
            float pplot = plot_dimension(ThermoMath.p_min, ThermoMath.p_max, p);
            if (pplot > highest_y) { highest_y = pplot; highest_y_i = mesh_positions.Count; }
            float tplot = plot_dimension(ThermoMath.t_min, ThermoMath.t_max, t);

            double v;
            float vplot;
            Vector3 point;

            v = ThermoMath.vliq_given_p(p, state.region);
            vplot = plot_dimension(ThermoMath.v_min, ThermoMath.v_max, v);
            point = new Vector3(vplot, pplot, tplot);
            mesh_positions.Add(point);

            v = ThermoMath.vvap_given_p(p, state.region);
            vplot = plot_dimension(ThermoMath.v_min, ThermoMath.v_max, v);
            point = new Vector3(vplot, pplot, tplot);
            mesh_positions.Add(point);
        }
        highest_y = Mathf.Lerp(highest_y, 1f, 0.01f); //extra nudge up

        //kill spanning triangles; gather orphans
        //"ladder"/"rung" terminology a bit arbitrary- attempts to keep track of each side of a "zipper" for each seam ("left" seam, "right" seam, each have own ladder/rung)
        List<int> left_orphans = new List<int>();
        List<int> right_orphans = new List<int>();
        int left_ladder_i = position_dome_region;
        Vector3 left_ladder = mesh_positions[left_ladder_i];
        Vector3 left_rung = mesh_positions[left_ladder_i + 2];
        int right_ladder_i = left_ladder_i + 1;
        Vector3 right_ladder = mesh_positions[right_ladder_i];
        Vector3 right_rung = mesh_positions[right_ladder_i + 2];
        for (var i = 0; i < mesh_triangles.Count; i += 3) {
            int ai = mesh_triangles[i + 0];
            int bi = mesh_triangles[i + 1];
            int ci = mesh_triangles[i + 2];
            Vector3 a = mesh_positions[ai];
            Vector3 b = mesh_positions[bi];
            Vector3 c = mesh_positions[ci];

            if ((left_rung.y < a.y || left_rung.y < b.y || left_rung.y < c.y) && left_ladder_i + 4 < mesh_positions.Count) { left_ladder_i += 2; left_ladder = mesh_positions[left_ladder_i]; left_rung = mesh_positions[left_ladder_i + 2]; }
            if ((right_rung.y < a.y || right_rung.y < b.y || right_rung.y < c.y) && right_ladder_i + 4 < mesh_positions.Count) { right_ladder_i += 2; right_ladder = mesh_positions[right_ladder_i]; right_rung = mesh_positions[right_ladder_i + 2]; }

            float x_cmp = (left_ladder.x + right_ladder.x) / 2f;
            if (
              (a.y < highest_y || b.y < highest_y || c.y < highest_y) &&
              (a.x < x_cmp || b.x < x_cmp || c.x < x_cmp) &&
              (a.x > x_cmp || b.x > x_cmp || c.x > x_cmp)
            ) {
                mesh_triangles.RemoveAt(i + 2);
                mesh_triangles.RemoveAt(i + 1);
                mesh_triangles.RemoveAt(i + 0);
                i -= 3;

                if (a.x < x_cmp && b.x < x_cmp) {
                    left_orphans.Add(ai);
                    left_orphans.Add(bi);
                    right_orphans.Add(ci);
                }
                else if (b.x < x_cmp && c.x < x_cmp) {
                    left_orphans.Add(bi);
                    left_orphans.Add(ci);
                    right_orphans.Add(ai);
                }
                else if (c.x < x_cmp && a.x < x_cmp) {
                    left_orphans.Add(ci);
                    left_orphans.Add(ai);
                    right_orphans.Add(bi);
                }
                else if (a.x < x_cmp) {
                    right_orphans.Add(bi);
                    right_orphans.Add(ci);
                    left_orphans.Add(ai);
                }
                else if (b.x < x_cmp) {
                    right_orphans.Add(ai);
                    right_orphans.Add(ci);
                    left_orphans.Add(bi);
                }
                else if (c.x < x_cmp) {
                    right_orphans.Add(ai);
                    right_orphans.Add(bi);
                    left_orphans.Add(ci);
                }
                else {
                    Debug.Log("NOOOO");
                }
            }
        }

        //sort orphans
        GRAPHPTCMP cmp = new GRAPHPTCMP(mesh_positions);

        left_orphans.Sort(cmp);
        for (int i = 1; i < left_orphans.Count; i++) { if (left_orphans[i - 1] == left_orphans[i]) { left_orphans.RemoveAt(i); i--; } }

        right_orphans.Sort(cmp);
        for (int i = 1; i < right_orphans.Count; i++) { if (right_orphans[i - 1] == right_orphans[i]) { right_orphans.RemoveAt(i); i--; } }

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
            rung = mesh_positions[ladder_i + 2];
            mesh_triangles.Add(ladder_i);
            mesh_triangles.Add(orphans[orphan_i]);
            orphan_i++;
            orphan = mesh_positions[orphans[orphan_i]];
            orung = mesh_positions[orphans[orphan_i + 1]];
            mesh_triangles.Add(orphans[orphan_i]);
            orphan = mesh_positions[orphans[orphan_i]];
            while (ladder_i + 2 < mesh_positions.Count) {
                while (orung.z <= rung.z && orung.y <= rung.y && orphan_i + 1 < orphans.Count) { //increment orphan
                    ai = ladder_i;
                    bi = orphans[orphan_i];
                    ci = orphans[orphan_i + 1];
                    mesh_triangles.Add(ai);
                    mesh_triangles.Add(bi);
                    mesh_triangles.Add(ci);

                    orphan_i++;
                    orphan = mesh_positions[orphans[orphan_i]];
                    if (orphan_i + 1 < orphans.Count) orung = mesh_positions[orphans[orphan_i + 1]]; //yes, both this AND previous line need +1 (+1 for advance, +1 for orung)
                }
                if (ladder_i + 2 < mesh_positions.Count) { //increment ladder
                    ai = ladder_i;
                    bi = orphans[orphan_i];
                    ci = ladder_i + 2;
                    mesh_triangles.Add(ai);
                    mesh_triangles.Add(bi);
                    mesh_triangles.Add(ci);

                    ladder_i += 2;
                    ladder = mesh_positions[ladder_i];
                    if (ladder_i + 2 < mesh_positions.Count) rung = mesh_positions[ladder_i + 2]; //yes, both this AND previous line need +2 (+2 for advance, +2 for rung)
                }
            }
            left_orphan_i = orphan_i;

            //right
            orphans = right_orphans;
            orphan_i = 0;
            orphan = mesh_positions[orphans[orphan_i]];
            orung = mesh_positions[orphans[orphan_i + 1]];
            ladder_i = position_dome_region + 1;
            ladder = mesh_positions[ladder_i];
            rung = mesh_positions[ladder_i + 2];
            mesh_triangles.Add(orphans[orphan_i]);
            mesh_triangles.Add(ladder_i);
            ladder_i += 2;
            ladder = mesh_positions[ladder_i];
            mesh_triangles.Add(ladder_i);
            while (ladder_i + 2 < mesh_positions.Count) {
                while ((ladder.y > orung.y || rung.z > orung.z) && orphan_i + 1 < orphans.Count) { //increment orphan
                    ai = orphans[orphan_i];
                    bi = ladder_i;
                    ci = orphans[orphan_i + 1];
                    mesh_triangles.Add(ai);
                    mesh_triangles.Add(bi);
                    mesh_triangles.Add(ci);

                    orphan_i++;
                    orphan = mesh_positions[orphans[orphan_i]];
                    if (orphan_i + 1 < orphans.Count) orung = mesh_positions[orphans[orphan_i + 1]]; //yes, both this AND previous line need +1 (+1 for advance, +1 for orung)
                }
                if (ladder_i + 2 < mesh_positions.Count) { //increment ladder
                    ai = orphans[orphan_i];
                    bi = ladder_i;
                    ci = ladder_i + 2;
                    mesh_triangles.Add(ai);
                    mesh_triangles.Add(bi);
                    mesh_triangles.Add(ci);

                    ladder_i += 2;
                    ladder = mesh_positions[ladder_i];
                    if (ladder_i + 2 < mesh_positions.Count) rung = mesh_positions[ladder_i + 2]; //yes, both this AND previous line need +2 (+2 for advance, +2 for rung)
                }
            }
            right_orphan_i = orphan_i;
        }

        //fan missing top
        for (int i = left_orphan_i + 1; i < left_orphans.Count; i++) {
            mesh_triangles.Add(left_orphans[i - 1]);
            mesh_triangles.Add(left_orphans[i]);
            mesh_triangles.Add(highest_y_i);
        }
        for (int i = right_orphan_i + 1; i < right_orphans.Count; i++) {
            mesh_triangles.Add(right_orphans[i]);
            mesh_triangles.Add(right_orphans[i - 1]);
            mesh_triangles.Add(highest_y_i);
        }
        mesh_triangles.Add(left_orphans[left_orphans.Count - 1]);
        mesh_triangles.Add(right_orphans[right_orphans.Count - 1]);
        mesh_triangles.Add(highest_y_i);

        //fill in vapor dome
        int triangle_inner_dome_region = mesh_triangles.Count;
        int position_inner_dome_region = mesh_positions.Count;
        for (int i = position_dome_region; i < position_inner_dome_region; i++) //duplicate inner positions so each can have own normal at seam
        {
            mesh_positions.Add(mesh_positions[i]);
        }
        for (int y = 0; y < concentrated_samples - 1; y++) {
            mesh_triangles.Add(position_inner_dome_region + y * 2 + 0);
            mesh_triangles.Add(position_inner_dome_region + y * 2 + 2);
            mesh_triangles.Add(position_inner_dome_region + y * 2 + 1);
            mesh_triangles.Add(position_inner_dome_region + y * 2 + 1);
            mesh_triangles.Add(position_inner_dome_region + y * 2 + 2);
            mesh_triangles.Add(position_inner_dome_region + y * 2 + 3);
        }


        //set normals
        mesh_normals = new List<Vector3>(new Vector3[mesh_positions.Count]);
        for (int i = 0; i < mesh_triangles.Count; i += 3) {
            int ai = mesh_triangles[i + 0];
            int bi = mesh_triangles[i + 1];
            int ci = mesh_triangles[i + 2];
            Vector3 a = mesh_positions[ai];
            Vector3 b = mesh_positions[bi];
            Vector3 c = mesh_positions[ci];
            Vector3 n = Vector3.Cross(Vector3.Normalize(b - a), Vector3.Normalize(c - a));
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
        graphObject.transform.localPosition = new Vector3(0f, 0f, 0f);
        graphObject.transform.localScale = new Vector3(1f, 1f, 1f);
        graphObject.transform.localRotation = Quaternion.identity;
        graphObject.GetComponent<MeshFilter>().mesh = mesh;
        graphObject.GetComponent<MeshRenderer>().material = graph_material;
        graphObject.GetComponent<MeshRenderer>().shadowCastingMode = UnityEngine.Rendering.ShadowCastingMode.Off;
        var lightable = graphObject.GetComponent<Lightable>();
        lightable.use_custom_mats = true;
        lightable.base_mat = graph_material;
        lightable.lit_mat = graph_material_lit;


        /*
            //add tests
            for(int y = 0; y < samples; y++)
            {
              double pt = ((double)y/(samples-1));
              for(int z = 0; z < samples; z++)
              {
                double tt = ((double)z/(samples-1));
                double pst = sample(pt);
                double tst = sample(tt);

                //pvt in Pa, M³/Kg, K
                double p = ThermoMath.p_given_percent(pst);
                double t = ThermoMath.t_given_percent(tst);
                double v = ThermoMath.v_given_pt(p,t, state.region);

                //Debug.LogFormat("p:{0}Pa, v:{1}m³/Kg, t:{2}K",p,v,t);
                int i = samples*y+z;
                pt_positions[i] = plot(p,v,t);
              }
            }
            genTestMeshFromPts(pt_positions, "v_given_pt"); //used to gen main mesh

            for(int y = 0; y < samples; y++)
            {
              double pt = ((double)y/(samples-1));
              for(int z = 0; z < samples; z++)
              {
                double tt = ((double)z/(samples-1));
                double pst = sample(pt);
                double tst = sample(tt);

                //pvt in Pa, M³/Kg, K
                //double p = ThermoMath.p_given_percent(pst);
                //double t = ThermoMath.t_given_percent(tst);
                //double v = ThermoMath.v_given_pt(p,t, state.region);

                double t = ThermoMath.t_given_percent(1.0-tst);
                double v = ThermoMath.v_given_percent(pst);
                double p = ThermoMath.p_given_vt(v,t);

                //Debug.LogFormat("p:{0}Pa, v:{1}m³/Kg, t:{2}K",p,v,t);
                int i = samples*y+z;
                pt_positions[i] = plot(p,v,t);
              }
            }
            genTestMeshFromPts(pt_positions, "p_given_vt");

        //*/
    }

    void findObjects() {
        graph = GameObject.Find("gmodel");
        state_dot = GameObject.Find("gstate");
    }

    public void debug_deltas() {
        state.debug_deltas(debug_write, debug_file);
    }

    //assume starting/ending point consistent for whole API!

    public Vector3 guessPlot(double guess_t, double pcube, double vcube) {
        //attempt to do it "the right way"
        /* //not worth the complexity
            Vector3 plot;
            plot.x = state.Clampf(invplot_dimension(ThermoMath.v_min, ThermoMath.v_max, vcube), (float)ThermoMath.v_min, (float)ThermoMath.v_max);
            plot.y = state.Clampf(invplot_dimension(ThermoMath.p_min, ThermoMath.p_max, pcube), (float)ThermoMath.p_min, (float)ThermoMath.p_max);
            plot.z = (float)ThermoMath.iterate_t_given_pv(guess_t, (double)plot.y, (double)plot.x); //gives off warnings because now these math functions are stateful (bad)
            int region = ThermoMath.region_given_pvt(state.pressure,volume,temperature);
            if (region != 1)
            {
              plot.y = (float)ThermoMath.p_given_vt((double)plot.x, (double)plot.z);
              plot.x = (float)ThermoMath.v_given_pt((double)plot.y, (double)plot.z);
            }
            else
            {
               //TODO: !!
            }
            if (plot.z == 0.0) plot = new Vector3((float)volume, (float)state.pressure, (float)temperature);
            return plot;
        */
        return new Vector3(0, 0, 0);
    }

    public Vector3 guessMeshPlot(double vcube, double pcube, double tcube) {
        Vector3 pt = new Vector3((float)vcube, (float)pcube, (float)tcube);
        float d = Vector3.SqrMagnitude(pt - mesh_positions[0]);
        int closest = 0;
        for (int i = 1; i < mesh_positions.Count; i++) {
            float od = Vector3.SqrMagnitude(pt - mesh_positions[i]);
            if (od < d) {
                d = od;
                closest = i;
            }
        }

        Vector3 plot;
        if (closest < position_dome_region) {
            plot.x = MathUtility.Clampf(invplot_dimension(ThermoMath.v_min, ThermoMath.v_max, mesh_positions[closest].x), (float)ThermoMath.v_min, (float)ThermoMath.v_max);
            plot.y = MathUtility.Clampf(invplot_dimension(ThermoMath.p_min, ThermoMath.p_max, mesh_positions[closest].y), (float)ThermoMath.p_min, (float)ThermoMath.p_max);
            plot.z = MathUtility.Clampf(invplot_dimension(ThermoMath.t_min, ThermoMath.t_max, mesh_positions[closest].z), (float)ThermoMath.t_min, (float)ThermoMath.t_max);
        }
        else {
            plot.x = MathUtility.Clampf(invplot_dimension(ThermoMath.v_min, ThermoMath.v_max, vcube), (float)ThermoMath.v_min, (float)ThermoMath.v_max);
            plot.y = MathUtility.Clampf(invplot_dimension(ThermoMath.p_min, ThermoMath.p_max, mesh_positions[closest].y), (float)ThermoMath.p_min, (float)ThermoMath.p_max);
            plot.z = MathUtility.Clampf(invplot_dimension(ThermoMath.t_min, ThermoMath.t_max, mesh_positions[closest].z), (float)ThermoMath.t_min, (float)ThermoMath.t_max);
        }
        return plot;
    }

    public void warp_pv_partial(double p, double v, double t) {
        state.warp_pv_partial(p, v, t, this);
        visualize_state();
    }

    public void warp_pv(double p, double v, double t) {
        state.warp_pv(p, v, t);
        visualize_state();
    }

    public void add_heat_per_delta_time(double applied_heat, double insulation_coefficient, double delta_time, double p_outside, bool is_internal, double temperature_gradient) {
        state.add_heat_per_delta_time(applied_heat, insulation_coefficient, delta_time, p_outside, is_internal, temperature_gradient);

        visualize_state();
    }

    public void add_pressure_uninsulated_per_delta_time(double p, double delta_time, double insulation_coefficient, double p_outside, double temperature_gradient) {
        state.add_pressure_uninsulated_per_delta_time(p, delta_time, insulation_coefficient, p_outside, temperature_gradient);
        visualize_state();
    }

    public void add_pressure_insulated_per_delta_time(double p, double delta_time, double p_outside, double temperature_gradient) {
        state.add_pressure_insulated_per_delta_time(p, delta_time, p_outside, temperature_gradient);
        visualize_state();
    }

    void visualize_state() {
        state_dot.transform.localPosition = plot(state.pressure, state.volume, state.temperature);

        update_tracker_pos();

        // PISTON HEIGHT

        double q = state.quality;

        // Since the difference between min and max volume is a magnitude of 10^7, linear scaling won't work.
        // Instead, use a log scale.
        // log_offset is the exponent value that produces ThermoMath.v_min * -1 given the piston radius.
        // Applying the offset ensures we never have a volume scale < 0
        double vheight = 0; // vapor height
        double lheight = 0; // liquid height
        double total_height = 0; // combined height

        switch (state.region) {
            case ThermoMath.region_liquid:
                lheight = Math.Log(state.volume / ThermoState.piston_area) + ThermoState.log_offset_volume;
                vheight = 0;
                break;
            case ThermoMath.region_twophase:
                // TODO: under constant volume, applying heat/cooling results in changing overall volume presentation
                double vliq = ThermoMath.vliq_given_p(state.pressure, state.region);
                double vvap = ThermoMath.vvap_given_p(state.pressure, state.region);
                vheight = Math.Log(q * vvap / ThermoState.piston_area) + ThermoState.log_offset_volume;
                lheight = Math.Log((1 - q) * vliq / ThermoState.piston_area) + ThermoState.log_offset_volume;
                if (lheight < 0) {
                    // turns out log offset works for overall volume... but individual pieces like liquid height
                    // can still be much smaller than the minimum overall state volume. But at these values,
                    // we're talking about slivers smaller than the player can even see. So we can ignore them.
                    lheight = 0;
                }
                // round out to current overall volume (take out from vapor)
                double dif = (vheight + lheight) - (Math.Log(state.volume / ThermoState.piston_area) + ThermoState.log_offset_volume);
                vheight -= dif;
                break;
            case ThermoMath.region_vapor:
                lheight = 0;
                vheight = Math.Log(state.volume / ThermoState.piston_area) + ThermoState.log_offset_volume;
                break;
            default:
                break;
        }

        total_height = vheight + lheight;

        // WATER/VAPOR RATIO

        float log_map = (float)(total_height / max_height_log); // map log height to range from 0 to 1

        float piston_span = piston.GetSpan(); // distance between min and max piston positions
        Vector3 new_piston_pos = piston.transform.localPosition;
        new_piston_pos.y = piston.GetMinPos().y + log_map * piston_span;
        piston.transform.localPosition = new_piston_pos;


        float contents_span = contents.GetSpan();
        Vector3 new_contents_scale = contents.transform.localScale;
        new_contents_scale.y = contents.GetMinScale().y + log_map * contents_span;
        contents.transform.localScale = new_contents_scale;

        Vector3 new_liquid_scale = contents.Water.transform.localScale;
        Vector3 new_vapor_scale = contents.Steam.transform.localScale;

        new_liquid_scale.y = (float)(lheight / total_height);
        new_vapor_scale.y = -(1 - new_liquid_scale.y);

        contents.Water.transform.localScale = new_liquid_scale;
        contents.Steam.transform.localScale = new_vapor_scale;
    }

    private void update_tracker_pos() {
        Vector3 state_pos = state_dot.transform.localPosition;

        Vector3 tracker_pos = dot_x_tracker.transform.localPosition;
        dot_x_tracker.transform.localPosition = new Vector3(state_pos.x, tracker_pos.y, tracker_pos.z);

        tracker_pos = dot_y_tracker.transform.localPosition;
        dot_y_tracker.transform.localPosition = new Vector3(tracker_pos.x, state_pos.y, tracker_pos.z);

        tracker_pos = dot_z_tracker.transform.localPosition;
        dot_z_tracker.transform.localPosition = new Vector3(tracker_pos.x, tracker_pos.y, state_pos.z);
    }


    string region_to_name(int region) //0 subcooled liquid, 1 two-phase, 2 superheated vapor
    {
        switch (region) {
            case ThermoMath.region_liquid: return "Subcooled Liquid";
            case ThermoMath.region_twophase: return "Two-Phase";
            case ThermoMath.region_vapor: return "Superheated Vapor";
        }
        return "Undefined";
    }

    private void Start() {
        string update_text = "";
        update_text = "region: " + region_to_name(state.region); DispatchText(update_text, "", state.region / 2.0f, VarID.Region);
        if (state.region == 1) { update_text = string.Format("x: {0:0.000}", (float)(state.quality * 100f)); DispatchText(update_text, Units.Quality, state.quality / quality_range, VarID.Quality); }
        else if (state.region == 0) { update_text = "x: Undefined"; DispatchText(update_text, "", 0, VarID.Quality); }
        else { update_text = "x: Undefined"; DispatchText(update_text, "", 1, VarID.Quality); }
    }

    // Update is called once per frame
    void Update() {
        //detect editor graphgen modifications
        bool modified = false;
        modified = ((plot_lbase != plot_lbase_prev) || (sample_lbase != sample_lbase_prev));
        sample_lbase_prev = sample_lbase;
        plot_lbase_prev = plot_lbase;
        if (modified) genMesh();

        string update_text = "";
        if (Math.Abs(state.pressure - state.prev_pressure) > ThermoMath.p_smallstep) { update_text = string.Format("P: " + DigitFormat.Pressure, (float)state.pressure / 1000f); DispatchText(update_text, Units.Pressure, (state.pressure - ThermoMath.p_min) / pressure_range, VarID.Pressure); }
        if (Math.Abs(state.temperature - state.prev_temperature) > ThermoMath.t_smallstep) { update_text = string.Format("T: " + DigitFormat.TemperatureK + "    " + Units.TemperatureK + "  ({1:0.00} " + Units.TemperatureC + ")", (float)state.temperature, (float)state.temperature - 273.15f); DispatchText(update_text, "", (state.temperature - ThermoMath.t_min) / temperature_range, VarID.Temperature); }
        if (Math.Abs(state.volume - state.prev_volume) > ThermoMath.v_smallstep) { update_text = string.Format("v: " + DigitFormat.Volume, (float)state.volume); DispatchText(update_text, Units.Volume, (state.volume - ThermoMath.v_min) / volume_range, VarID.Volume); }
        if (Math.Abs(state.internalenergy - state.prev_internalenergy) > ThermoMath.u_smallstep) { update_text = string.Format("u: " + DigitFormat.InternalEnergy, (float)state.internalenergy / 1000f); DispatchText(update_text, Units.InternalEnergy, state.internalenergy / internal_energy_range, VarID.InternalEnergy); }
        if (Math.Abs(state.entropy - state.prev_entropy) > ThermoMath.s_smallstep) { update_text = string.Format("s: " + DigitFormat.Entropy, (float)state.entropy / 1000f); DispatchText(update_text, Units.Entropy, state.entropy / entropy_range, VarID.Entropy); }
        if (Math.Abs(state.enthalpy - state.prev_enthalpy) > ThermoMath.h_smallstep) { update_text = string.Format("h: " + DigitFormat.Enthalpy, (float)state.enthalpy / 1000f); DispatchText(update_text, Units.Enthalpy, state.enthalpy / enthalpy_range, VarID.Enthalpy); }
        if (state.region == 1 && Math.Abs(state.quality - state.prev_quality) > ThermoMath.x_smallstep) { update_text = string.Format("x: " + DigitFormat.Quality, (float)(state.quality * 100f)); DispatchText(update_text, Units.Quality, state.quality / quality_range, VarID.Quality); }
        if (true /*state.region != state.prev_region*/) {
            update_text = "Region: " + region_to_name(state.region); DispatchText(update_text, "", state.region / 2.0f, VarID.Region);
            if (state.region == 1) { update_text = string.Format("x: " + DigitFormat.Quality, (float)(state.quality * 100f)); DispatchText(update_text, Units.Quality, state.quality / quality_range, VarID.Quality); }
            else if (state.region == 0) { update_text = "x: Undefined"; DispatchText(update_text, "", 0, VarID.Quality); }
            else { update_text = "x: Undefined"; DispatchText(update_text, "", 1, VarID.Quality); }
        }

        state.stamp_prev();
    }

    private void DispatchText(string update_text, string units, double proportion, VarID varId) {
        GameMgr.Events.Dispatch(GameEvents.UpdateVarText, new VarUpdate(varId, update_text, units, (float)proportion));
    }

}

