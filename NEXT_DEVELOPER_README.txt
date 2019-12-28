edit: I started writing this as an attempt to contextualize the subtle ways the project is incomplete, but by the end of writing this I'm realizing it's actually not nearly as bad as I thought. feel free to read on if you want a thorough justification of the situation.

The state of Thermo VR, as of 12/27/19
To anyone next working on this project, it will be EXTREMELY HELPFUL if you read through this whole thing.

The premise of the application was the following:
In pet thermodynamics problems dealing with the states of water (as used for practice, on tests, etc...), there is a common abstract mechanism. This mechanism consists of an open cylindrical container, containing some volume of water, which is sealed by a free-floating piston. The mechanism is "abstract" in that 1. it doesn't actually exist anywhere (it's a totally useless device), and 2. certain physical tolerances can be ignored (the container is "perfectly sealed" by the piston, the piston does NOT experience friction, mechanisms of heat transfer can be assumed to be idealized, etc...).
This application was going to model this mechanism in VR, giving the user the ability to "actually interact with" what was previously only on paper.
The available processes to play with were to be: add heat, remove heat, add pressure, remove pressure.


THE FOLLOWING SOUNDS LIKE REALLY BAD NEWS, BUT IT'S THE BEST WAY I COULD THINK TO CLEARLY EXPLAIN WHERE THINGS ARE AND HOW THEY GOT THERE:

I will now show the faulty reasoning with which we started the project, and afterwards, will contrast these with what is currently understood.

The development of the simulation was predicated on a few (INCORRECT) ideas:
- given any two of P,V,T, the other can be derived (given a specified black-box thermo implementation)
- adding weight increases P while keeping T the same
- adding heat increases T while keeping EITHER P the same [free-floating piston] OR keeping V the same [piston clamped into position]
- we have the ability to easily move from any "correct" PVT to a new "correct" PVT for any of the given interactions.

Pseudocode for the (INCORRECT) processes:
  (meant to be viewed with a fixed-width font) 
  Change Weight:                      new_p = old_p+delta; new_t = old_t; new_v = v_given_pt(new_p,new_t);
  Change Temp (free-floating piston): new_t = old_t+delta; new_p = old_p; new_v = v_given_pt(new_p,new_t);
  Change Temp (clamped       piston): new_t = old_t+delta; new_v = old_v; new_p = p_given_vt(new_v,new_t);


Here are the corrected critiques on this reasoning:

The relevant variables are actually p,v,t,u,s,h,x. "Any two" of them are NOT sufficient to derive the others.
- p and t are insufficient to derive anything else in the case where water is in the "vapor dome" state (part-way between liquid and vapor)
- x (and anything else) is insufficient anywhere EXCEPT the "vapor dome"
- there may be more exceptions I am unaware of?

The specified black-box thermo implementation (IAPWS) actually has multiple specifications (at least a 95 version, and a 97 version), and MANY implementations (across many different languages).
- each specification/implementation combo deals with some subset of the various possible X_given_YZ(Y,Z) [where X,Y,Z are any of p,v,t,u,s,h,x]
- given the two implementations we found/were able to port (IAPWS95.cs and IAPWS97.cs), we have the ability to directly derive MANY of the possible combinations, but NOT all
- the rest (given they are possible) are able to be approximated from our available subset iteratively (see ThermoMath.cs)
- (fortunately we do not NEED every possible X_given_YZ)

The processes specified DO NOT correspond to the previously stated outcomes.
- "adding heat" does NOT mean "increasing temperature" at all ("heat" != "temperature")
- "adding pressure" does NOT imply that "temperature stays the same" (primarily relevant in the vapor dome)

If we translate the processes to human actions ("add pressure" = "place a weight on the piston"; "add heat" = "put a bunson burner under the container"; etc...), they are often only "solvable" in ways that much better fit a word problem than they do a realtime simulation.
(In other words, "given any* two variables, the rest can be derived" is only generally true assuming a timescale that allows the system to reach equilibrium- a simulator performs calculations every 1/60th of a second)
- "adding heat without an insulator" is thought of as "not adding heat" (because over time, that heat would fully dissipate into the environment anyways)
- "adding weight while in the vapor dome" is "unsolvable" unless there is a long-term known resting position (ie, there is a "stop" that you know the piston will compress to and hit, resulting in a known final volume, etc...)


The previous ideas, stated correctly:
- There are many states of water such that knowing two of p,v,t,u,s,h,x, the rest can be derived for its eventual state of equilibrium.
- There are specifications that define many of the "X_given_YZ" derivations for many different states.
- Much of the rest can be determined iteratively from the subset available

As an example of something that CANNOT be calculated:
The piston is free-floating, there is a perfect insulator on the container, and the water is in the vapor dome. Some weight is added, insufficient to get it out of the vapor dome. Finding the new state is not possible. It's _especially_ not possible with the caveat of "exactly 1/60th of a second has passed".


The current state of things:
- Free-floating piston, add heat
  perfectly calculable everywhere
- Clamped piston, add heat
  perfectly calculable everywhere
- Uninsulated, add pressure
  only calculable when vapor
- Insulated, add pressure
  only calculable when vapor

What needs to be done:
- keep weight and balloon static when not a superheated vapor
- apply weight deltas as iterating toward target


