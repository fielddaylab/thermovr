using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class World : MonoBehaviour
{

  GameObject room;
  GameObject graph;
  GameObject vessel;
  GameObject lhand;
  GameObject rhand;

  public Material hand_empty;
  public Material hand_intersecting;
  public Material hand_grabbing;

  List<Grabbable> grabbables;
  GameObject lgrabbed = null;
  GameObject rgrabbed = null;

  bool ltrigger = false;
  bool rtrigger = false;

  // Start is called before the first frame update
  void Start()
  {
    room   = GameObject.Find("Room");
    graph  = GameObject.Find("Graph");
    vessel = GameObject.Find("Vessel");
    lhand  = GameObject.Find("LHand");
    rhand  = GameObject.Find("RHand");

    lhand.GetComponent<MeshRenderer>().material = hand_empty;
    rhand.GetComponent<MeshRenderer>().material = hand_empty;

    grabbables = new List<Grabbable>();
    grabbables.Add(graph.GetComponent<Grabbable>());
    grabbables.Add(vessel.GetComponent<Grabbable>());
  }

  // Update is called once per frame
  void Update()
  {
    float trigger_threshhold = 0.1f;

    //left hand
    int ltrigger_delta = 0;
    if(!ltrigger && OVRInput.Get(OVRInput.Axis1D.PrimaryHandTrigger) > trigger_threshhold)
    {
      ltrigger_delta = 1;
      ltrigger = true;
    }
    else if(ltrigger && OVRInput.Get(OVRInput.Axis1D.PrimaryHandTrigger) <= trigger_threshhold)
    {
      ltrigger_delta = -1;
      ltrigger = false;
    }

    if(lgrabbed == null && ltrigger_delta == 1)
    {
      for(int i = 0; lgrabbed == null && i < grabbables.Count; i++)
      {
        if(grabbables[i].lintersect)
        {
          lgrabbed = grabbables[i].gameObject;
          lgrabbed.transform.SetParent(lhand.transform);
          lhand.GetComponent<MeshRenderer>().material = hand_grabbing;
          if(lgrabbed == rgrabbed)
          {
            rgrabbed = null;
            rhand.GetComponent<MeshRenderer>().material = hand_intersecting;
          }
        }
      }
    }
    else if(lgrabbed && ltrigger_delta == -1)
    {
      lgrabbed.transform.SetParent(room.transform);
      lgrabbed = null;
    }
    if(lgrabbed == null)
    {
      lhand.GetComponent<MeshRenderer>().material = hand_empty;
      for(int i = 0; i < grabbables.Count; i++)
      {
        if(grabbables[i].lintersect)
        {
          lhand.GetComponent<MeshRenderer>().material = hand_intersecting;
          break;
        }
      }
    }

    //right hand
    int rtrigger_delta = 0;
    if(!rtrigger && OVRInput.Get(OVRInput.Axis1D.SecondaryHandTrigger) > trigger_threshhold)
    {
      rtrigger_delta = 1;
      rtrigger = true;
    }
    else if(rtrigger && OVRInput.Get(OVRInput.Axis1D.SecondaryHandTrigger) <= trigger_threshhold)
    {
      rtrigger_delta = -1;
      rtrigger = false;
    }
    if(rgrabbed == null && rtrigger_delta == 1)
    {
      for(int i = 0; rgrabbed == null && i < grabbables.Count; i++)
      {
        if(grabbables[i].rintersect)
        {
          rgrabbed = grabbables[i].gameObject;
          rgrabbed.transform.SetParent(rhand.transform);
          rhand.GetComponent<MeshRenderer>().material = hand_grabbing;
          if(rgrabbed == lgrabbed)
          {
            lgrabbed = null;
            lhand.GetComponent<MeshRenderer>().material = hand_intersecting;
          }
        }
      }
    }
    else if(rgrabbed && rtrigger_delta == -1)
    {
      rgrabbed.transform.SetParent(room.transform);
      rgrabbed = null;
    }
    if(rgrabbed == null)
    {
      rhand.GetComponent<MeshRenderer>().material = hand_empty;
      for(int i = 0; i < grabbables.Count; i++)
      {
        if(grabbables[i].rintersect)
        {
          rhand.GetComponent<MeshRenderer>().material = hand_intersecting;
          break;
        }
      }
    }

  }

}

