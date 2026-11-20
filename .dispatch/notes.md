# Loading progress part 1 — make Start init yield across frames + instrument

Parent card: `1901208238947304559`. This card is the mechanical enabler; the visible
fill bar is part 2.

**No WebGL build was run and no branch was pushed** (per the card's DO NOT). Everything
below was verified editor-free on this host.

---

## 1. The init chain now yields (AC 1)

### `Assets/Scripts/General/GameMgr.cs`

- **`Start()` — lines 41-48.** Was the entire synchronous chain
  (`m_thermo_present.Init()` → `m_world.Init()` → `m_thermo_present.Reset()` → four
  `EventMgr` dispatches) inside one frame. Now it is just
  `StartCoroutine(InitRoutine())`.
- **`InitRoutine()` — lines 50-108.** Phases separated by `yield return null`:
  - `yield return null` at **line 69** (after `ThermoPresent.InitRoutine()` drains)
  - `yield return null` at **line 76** (after `world.init`)
  - `yield return null` at **line 83** (after `thermo_present.reset`)
  - plus every yield forwarded out of `ThermoPresent.InitRoutine()` at **lines 64-67**
    (`while (thermo_init.MoveNext()) yield return thermo_init.Current;`)
- **`FixedUpdate()` — lines 110-120.** New `m_initComplete` guard (field at line 27, set
  at line 87). Necessary because `ManualFixedUpdate()` used to be unreachable before
  `World.Init()` had run — init finished inside `Start()`. Now that init spans frames,
  `FixedUpdate` can fire first, so it is held off until the chain completes. This is a
  required consequence of yielding, not an unrelated change.

### `Assets/Scripts/General/ThermoPresent.cs`

- **`LoadProgress` / `LoadComplete` / `CountInitWorkUnits()` — lines 102-124.**
- **`Init()` — lines 127-130.** Kept as a synchronous wrapper that drains the coroutine
  (`while (routine.MoveNext()) { }`) so any non-coroutine caller still works.
- **`InitRoutine()` — lines 132-181.** `yield return null` at **lines 158, 171** and the
  forwarded genMesh yields at **lines 160-163**.
- **`genMesh()` — lines 412-422.** Now a synchronous wrapper draining `genMeshRoutine()`.
  This preserves the editor-parameter-change call site in `UpdatePropertyReadouts()`
  (`if (modified) genMesh();`), which is not in a coroutine.
- **`genMeshRoutine()` — lines 424-830.** The chunked version.

### Both genMesh hot loops are chunked (AC 1)

| loop | lines | chunk boundary | chunks at `samples = 350` |
|---|---|---|---|
| surface sweep (350 rows × 350 `v_given_pt` solves) | **446-478** (boundary at 462-473) | every `SURFACE_ROWS_PER_CHUNK = 50` rows (line 389) | 7 |
| spanning-triangle cull (quadratic `RemoveAt`) | **562-630** (boundary at 563-573) | every `CULL_TRIS_PER_CHUNK = 40000` triangles examined (line 390) | ≤ 7 |

**Yield count across the whole init** (within the card's "roughly 10-20" target):
6 mid-sweep + ≤6 mid-cull + 1 after surface_sweep + 1 after dome + 1 after cull +
1 after stitch_normals + 1 after findObjects + 1 after state.reset + 3 in `GameMgr`
= **~22 yields**, i.e. `Update()` runs ~22 times during the load instead of 0.

---

## 2. `LoadProgress` is computed from completed work units (AC 2)

The math lives in a new pure-C# file, `Assets/Scripts/General/LoadProgressTracker.cs`
(no `using UnityEngine`, so it can be reasoned about and exercised off-editor).

**The computing lines, quoted verbatim** (`LoadProgressTracker.cs` lines 76-87):

```csharp
public static float Progress {
    get {
        if (s_completed >= s_total) {
            return 1.0f;
        }
        if (s_completed <= 0) {
            return 0.0f;
        }
        return (float)s_completed / (float)s_total;
    }
}
```

`s_completed` only ever moves up (`Complete()` lines 35-42 early-returns on
`units <= 0` and clamps with `Math.Min(s_total, ...)`; `CompleteTo()` lines 49-54
early-returns unless it is an increase). So `Progress` is **monotonic non-decreasing by
construction**, and it returns **exactly the literal `1.0f`** once
`s_completed >= s_total` — not a float division that happens to round to 1.

### Where the denominator comes from

Fixed up front, before any unit is completed, so the fraction can never move backwards
(`GameMgr.cs` lines 58-62):

```csharp
const int GAMEMGR_UNITS = 3; // world.init, thermo_present.reset, dispatch
LoadProgressTracker.Begin(m_thermo_present.CountInitWorkUnits() + GAMEMGR_UNITS);
```

`ThermoPresent.CountInitWorkUnits()` (line 121) `= 3 + CountGenMeshWorkUnits()`
(findObjects, state.reset, visualize_state), and `CountGenMeshWorkUnits()` (line 409)
`= SurfaceChunkBudget() + CullChunkBudget() + 3` (dome, stitch_normals, mesh_upload).
At `samples = 350` that is **23 units total**.

The cull budget is an *upper bound* (`mesh_triangles` only ever shrinks, so the initial
triangle count bounds the iterations). When a chunked loop finishes in fewer chunks than
its bound, the catch-up calls at **ThermoPresent.cs 475-478 and 632-635** award the
unused reserve, so the denominator stays honest and progress still only rises.

`LoadProgressTracker.Finish()` (`GameMgr.cs` line 105) pins it to exactly `1.0f` after
the final dispatch phase.

**Exposed as `public static float LoadProgress` in three places**, all reading the same
tracker: `ThermoPresent` line 107 (the name the card asks for), `GameMgr` line 37, and
`UILoading` line 130 (what part 2's fill bar should bind to).

### Verified by running it

Compiled `LoadProgressTracker.cs` standalone on net8.0 and replayed the real init
sequence (including a cull loop that finishes early, to exercise the catch-up path):

```
begin                        0/23 = 0.0000
findObjects                  1/23 = 0.0435
surface chunk 1..6           2..7/23
surface catch-up             8/23 = 0.3478
dome                         9/23 = 0.3913
cull chunk 1..4             10..13/23
cull catch-up               16/23 = 0.6957
stitch_normals              17/23 = 0.7391
mesh_upload                 18/23 = 0.7826
state.reset                 19/23 = 0.8261
visualize_state             20/23 = 0.8696
world.init                  21/23 = 0.9130
thermo_present.reset        22/23 = 0.9565
Finish                      23/23 = 1.0000

reaches exactly 1.0f: True
IsComplete: True
after bogus calls still 1.0f: True      (Complete(-5), Complete(0), CompleteTo(1))
over-complete clamps to 1.0f: True      (Begin(10) then Complete(999))
```

A monotonicity assertion ran on every step (it exits non-zero on any decrease) and
never tripped.

---

## 3. Instrumentation (scope item 1)

Every phase is wrapped in a `System.Diagnostics.Stopwatch` and logs one
`[LoadTiming] phase=<name> ms=<n>` line.

| phase | file : line |
|---|---|
| `thermomath.init` | ThermoPresent.cs 94 |
| `awake.state.reset` | ThermoPresent.cs 98 |
| `findObjects` | ThermoPresent.cs 155 |
| `surface_sweep` | ThermoPresent.cs 479 |
| `tri_index` | ThermoPresent.cs 506 |
| `dome` | ThermoPresent.cs 542 |
| `cull` | ThermoPresent.cs 636 |
| `stitch_normals` | ThermoPresent.cs 798 |
| `mesh_upload` | ThermoPresent.cs 826 |
| `genMesh` (whole) | ThermoPresent.cs 828 |
| `state.reset` | ThermoPresent.cs 168 |
| `visualize_state` | ThermoPresent.cs 176 |
| `world.init` | GameMgr.cs 73 |
| `thermo_present.reset` | GameMgr.cs 80 |
| `dispatch` | GameMgr.cs 101 |
| `total ms=<n>` | GameMgr.cs 107 |
| `verts=<n> tris=<n>` | ThermoPresent.cs 827 |

The per-phase stopwatches are `Stop()`ped around each `yield return null` and restarted
afterwards, so a phase's reported ms is *compute* time, not wall-clock including the
frames it waited out. `total ms` is deliberately *not* stopped across yields — it is
the user-visible wall clock, which is the number djgagnon is complaining about.

Two phases beyond the card's list are logged because the instrumentation exposed them:
`thermomath.init` (in `Awake()`, the 17 ms the card had already ruled out — logged so
the next reader does not have to re-derive that) and `tri_index` (the index-buffer build
between the sweep and the dome, previously unattributed).

---

## 4. Mesh output is bit-identical (AC 4) — how I confirmed it by reading the code

Counts cannot be run here: there is no Unity editor or licence on this host, so
`genMesh` cannot execute (`GameObject.Find`, `Mesh`, `Destroy`, `MeshFilter`).
The argument is therefore a code-reading one, plus an executable model of the two loop
rewrites.

**What the edits are.** In both loops the *only* inserted statements are a counter
increment, an `if` on that counter, a `LoadProgressTracker.Complete()`, a stopwatch
`Stop()`/`Start()`, and a `yield return null`. No expression that feeds
`mesh_positions`, `mesh_triangles`, `pt_positions`, the ladder/rung variables, the
orphan lists, `highest_y`, `highest_y_i` or `position_dome_region` is touched. In
particular:

- **Surface sweep.** The boundary is placed *after* the inner `z` loop closes
  (ThermoPresent.cs 460-473), i.e. only on whole-row boundaries. `y` and `z` keep their
  original bounds and step. `pt_positions[i]` is written at `i = samples * y + z` —
  an absolute index, not an append — so even a reordering (which does not happen) could
  not change the array. The guard `&& y + 1 < samples` suppresses a pointless trailing
  yield.
- **Cull.** The boundary is the *first* statement of the loop body (ThermoPresent.cs
  563-574), before `ai/bi/ci` are read. It does not touch `i`, so the
  `RemoveAt(i+2); RemoveAt(i+1); RemoveAt(i+0); i -= 3;` compaction and the
  `i < mesh_triangles.Count` re-evaluation are untouched. The ladder/rung advance and
  the six-way orphan classification are byte-for-byte the original lines.
- Converting a method to `IEnumerator` does not change the order of its statements;
  `yield return null` only suspends. All locals are hoisted into the compiler-generated
  state machine and survive across the suspension, so carried loop state (ladders,
  rungs, `highest_y`, orphan lists) persists exactly as before.
- `genMesh()` is retained as a wrapper that drains the routine to completion inside one
  frame (`while (routine.MoveNext()) { }`), so the `UpdatePropertyReadouts()` call site
  sees identical behaviour, including identical synchronicity.

**Executable corroboration.** I extracted both loop *shapes* — same control flow, same
`RemoveAt`+`i -= 3` compaction, same carried ladder-style state, with a deterministic
stand-in for the geometric predicate — and ran the original against the chunked version
on net8.0 with a deliberately tiny chunk size (7 triangles / 5 rows) to force 138 and 4
boundaries respectively:

```
cull:  tris identical=True (774 tris kept), left identical=True,
       right identical=True, carried-state identical=True, yields=138
sweep: identical=True (529 pts), yields=4
```

Element-by-element equality of the kept triangle list, both orphan lists, the carried
state, and the full point array. This is evidence about the *rewrite*, not about
thermodynamics — which is the right target, since the thermodynamics is unchanged.

**What was not changed:** `samples` is still `350` (ThermoPresent.cs line 47).
`ThermoMath.cs`, `IAPWS97.cs`, `IAPWS95.cs` and `ThermoState`'s solver paths are
untouched — `git diff --stat` lists only `GameMgr.cs`, `ThermoPresent.cs`,
`UILoading.cs` and the new `LoadProgressTracker.cs`. No collision with card
`1901120701390128201`. No prefab, no scene YAML.

---

## 5. `UILoading` no longer hides `m_loadGroup` on Start-interactable (AC 3)

`Assets/Scripts/UI/UILoading.cs`, old lines 97-100 → now **lines 100-125**.

Was:

```csharp
if (ModeMgr.Instance.IsDesktop && m_startButton.interactable)
{
    m_loadGroup.gameObject.SetActive(false);
}
```

The condition `m_startButton.interactable` stays true forever once
`HandleAsyncLoadComplete()` fires, so the loading visuals were destroyed at precisely
the moment the freeze began. Now the hide condition is split in two, keyed on a new
`m_startClicked` flag (field line 42, set in `HandleStartGameClicked()` line 159):

```csharp
bool awaiting_start_click = ModeMgr.Instance.IsDesktop && m_startButton.interactable && !m_startClicked;
bool init_finished = ModeMgr.Instance.IsDesktop && m_startClicked && ThermoVR.LoadProgressTracker.IsComplete;
```

- *Before* the click: hidden, as before — the title screen still looks the same while
  idle at the Start button.
- *From the click until `LoadProgress == 1f`*: shown, and the `else` branch keeps
  calling `m_loadIcon.transform.Rotate(...)` (line 123), which now actually runs because
  `Update()` is reached between phases.
- *After init*: hidden.

`HandleStartGameClicked()` also calls `m_loadGroup.gameObject.SetActive(true)` directly
(line 163) rather than waiting for the next `Update()`, since the scene activates
immediately on the click.

`IsComplete` is used rather than `Progress >= 1f` so the group is never hidden before a
load has even been started (`HasStarted` is false → `IsComplete` false).

---

## 6. Verification actually performed on this host (AC 5, AC 6)

### `tools/thermo-tests` (AC 5) — passes 15/15, but read the caveat

**`tools/thermo-tests` does not exist on this branch, and it cannot be made to compile
against this branch.** Both facts matter, so here they are precisely.

It was added by commit `36578c7`, which lives only on
`card-1901120701390128201-fix-entropy-drifts-and-t-is-pinned-in-superhe` — the staging
card this card is explicitly told not to collide with. It has not reached `desktop`, and
`desktop` is what this branch is cut from. So the registered command as literally
written fails at `cd tools/thermo-tests`: no such directory.

My first attempt was to check that branch's harness into a scratch tree and point it at
**this** worktree's `Assets/` (the csproj pulls production sources in via
`<Compile Include="../../Assets/Scripts/General/...">`). That does not build, and the
failure is informative rather than my fault:

```
ThermoMathPsTests.cs(58,29): error CS0117: 'ThermoMath' does not contain a definition for 't_given_ps'
ThermoMathPsTests.cs(93,29): error CS0117: 'ThermoMath' does not contain a definition for 'v_given_ps'
... (9 such errors)
```

`t_given_ps` / `v_given_ps` are added by *that* card's implementation commit
(`9c0b324`). The harness was written to fail to compile against pre-fix `ThermoMath`
— that is deliberate, it is how it serves as that card's regression gate. This branch
has pre-fix `ThermoMath`, so the harness cannot bind. Nothing I changed caused this.

So I ran the suite the only way it is meaningful: against that branch's complete tree,
with exactly the registered environment.

```
export DOTNET_ROOT=/opt/data/tools/dotnet; export PATH=$DOTNET_ROOT:$PATH
export DOTNET_CLI_TELEMETRY_OPTOUT=1; export DOTNET_SYSTEM_GLOBALIZATION_INVARIANT=1
export HOME=/opt/data
cd <scratch>/tools/thermo-tests && dotnet test
```

```
Passed!  - Failed: 0, Passed: 15, Skipped: 0, Total: 15, Duration: 26 ms - ThermoMathTests.dll (net8.0)
```

`DOTNET_SYSTEM_GLOBALIZATION_INVARIANT=1` is genuinely mandatory — without it the SDK
aborts with a core dump (no libicu, no root), as the card warns.

**What that result does and does not establish.** It is a green run of the thermo suite,
but it is a run of *the other branch's* sources, so on its own it says nothing about my
change. The claim it supports is made valid by a separate check: the three production
files that suite compiles are byte-identical between `HEAD` and my worktree, i.e. I did
not touch any of them.

```
Assets/Scripts/General/IAPWS97.cs     HEAD=02404bd536f9222f  worktree=02404bd536f9222f  IDENTICAL
Assets/Scripts/General/IAPWS95.cs     HEAD=db1ce45e4154710a  worktree=db1ce45e4154710a  IDENTICAL
Assets/Scripts/General/ThermoMath.cs  HEAD=48da5c9782959d99  worktree=48da5c9782959d99  IDENTICAL
```

Since the suite's entire compiled surface is unchanged by this commit, it is impossible
for this commit to alter its outcome. That — not the 15/15 line by itself — is the real
AC 5 evidence here.

And as the card already states: **these tests cover `ThermoMath` only and verify none of
criteria 1-4.**

### The changed C# parses — but there is no editor-free full compile (AC 6)

Stated plainly: **I could not compile the Unity-dependent sources, and I am not implying
a build passed.** There is no editor-free compile path for anything that references
`UnityEngine`, and the host has no Unity.

What I did instead, per file, with Roslyn (`csc.dll` from the 8.0.427 SDK):

```
Assets/Scripts/General/LoadProgressTracker.cs: no syntax/semantic errors other than missing refs
Assets/Scripts/General/GameMgr.cs:             ditto
Assets/Scripts/General/ThermoPresent.cs:       ditto
Assets/Scripts/UI/UILoading.cs:                ditto
```

Every error reported was `CS0246`/`CS0518`/`CS0234`/`CS0012` — "type or namespace not
found" for `UnityEngine`, `TMPro`, `BeauUtil` etc. **Zero** errors outside that set,
which is what rules out the failure modes this change could plausibly introduce:
unbalanced braces, a `yield return` in a non-iterator, `yield` inside a `try`/`finally`,
a duplicate or out-of-scope local from the hoisted chunk counters, or a bad member
declaration.

`LoadProgressTracker.cs` has no Unity dependency, so it **fully compiled** against
`net8.0` reference assemblies (`tracker.dll built OK`) and then actually ran — that is
where the progress table in section 2 comes from.

### Not verifiable without Unity — explicitly not claimed

- **That frames actually render during init.** Needs the WebGL build, which is part 2's
  job. I claim only that the chain now `yield return null`s between phases and that
  `Update()` is therefore reachable ~22 times during the load. **I am not claiming the
  freeze is fixed.**
- **Real phase timings.** The `[LoadTiming]` lines are the instrument; nobody has yet
  seen their output, because that needs a Unity run. The card's baseline (`genMesh`
  dominant, `surface_sweep` ~344 ms, `cull` ~369 ms, `ThermoMath.Init()` ~17 ms) was
  measured standalone-native and is unverified under IL2CPP/WebGL. Trust the
  distribution, not the absolutes, until the logs come back.
- **Actual vertex/triangle counts.** The `verts=`/`tris=` line will print them on first
  run; `genMesh` cannot execute here. Part 2 (or anyone with an editor) should capture
  one `[LoadTiming]` log before and after this commit and diff the `verts=`/`tris=`
  line — that is the empirical form of AC 4, and it is cheap once a build exists.
- **That the spinner is visually where the user expects it — this is the biggest open
  risk, and it is scene-level, so it belongs to part 2.** Two facts from reading the
  scenes (not editable here — the card forbids scene YAML):
  1. `Boot.LoadAsync()` calls `SceneManager.LoadSceneAsync(m_firstScene)` with the
     default `LoadSceneMode.Single`. When `allowSceneActivation` flips to true on the
     Start click, Unity unloads `LaunchScreen` — including the `UILoading` instance that
     owns the spinning `m_loadIcon` — before `GameMgr.Start()` in `Desktop_Layout` runs.
     The `UILoading.cs` fix in this card is correct and necessary (the old 97-100 logic
     would have hidden the group even if the object survived), but on its own it cannot
     make the LaunchScreen spinner outlive the scene swap.
  2. `Desktop_Layout.unity` has its *own* `Loading Canvas` prefab instance
     (`&339961287`) and `GameMgr.m_UILoading` is already wired to it (fileID
     `1180044090`), but that instance's root carries an `m_IsActive: 0` override, and the
     `m_UILoading` field is unused in code. Enabling it during init would re-run
     `UILoading.Awake()` (re-dispatching `TitleScreenDisplayed`, which is an analytics
     event), so I did not do it blind.
  Part 2 should pick one of: mark the LaunchScreen canvas `DontDestroyOnLoad` /
  `LoadSceneMode.Additive`, or activate the `Desktop_Layout` instance with its
  title-screen side effects suppressed. Either way the progress source
  (`LoadProgressTracker.Progress`, surfaced as `GameMgr.LoadProgress`,
  `ThermoPresent.LoadProgress`, `UILoading.LoadProgress`) is static, so it is reachable
  from whichever `UILoading` instance survives.
- **`FixedUpdate` gating side effects.** The `m_initComplete` guard is logically
  required, but whether any system depends on `ManualFixedUpdate()` running in those
  first ~22 frames is not observable without a run. Reading `ManualFixedUpdate()`, it
  only stabilises hands, applies tools, enforces limits and dispatches readouts — all of
  which were already unreachable before init finished, since init used to complete
  inside `Start()`.
