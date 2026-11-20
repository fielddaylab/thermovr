# ThermoMath (p,s) test harness

Standalone xUnit tests for the isentropic-path fix in the 100%-insulated
superheated-vapor pressure solver. No Unity install required: the csproj
compiles `IAPWS97.cs`, `IAPWS95.cs` and `ThermoMath.cs` straight out of
`Assets/` and `UnityShim.cs` supplies the handful of UnityEngine /
Unity.IL2CPP symbols they reference.

## Run

    dotnet test

On a host without the .NET SDK, install it without root:

    curl -sSL https://dot.net/v1/dotnet-install.sh | bash -s -- --channel 8.0 --install-dir ~/dotnet --no-path
    export DOTNET_ROOT=~/dotnet PATH=~/dotnet:$PATH

If the host has no ICU libraries (`Couldn't find a valid ICU package`) and you
cannot install libicu, run in invariant mode:

    export DOTNET_SYSTEM_GLOBALIZATION_INVARIANT=1

## What is covered

- `t_given_ps` / `v_given_ps` reproduce IAPWS-IF97 backward-function reference
  values in region 2 (vapor) and region 1 (liquid).
- Unit conventions: the wrappers take Pa and J/(kg K) and return K and m^3/kg,
  matching the rest of `ThermoMath`.
- **AC2 regression for the reported bug** (`IsentropicVaporPathTests`): from
  p = 0.5 MPa, T = 500 K, an insulated pressure change conserves entropy
  (<0.1% drift) while temperature tracks the IF97 isentropic path to within
  1 K (511.30 / 549.73 / 587.65 K at 0.55 / 0.75 / 1.00 MPa). Also asserts
  compression heats and expansion cools.

`if97_region2.py` is an independent IF97 oracle written from the release
equations (not ported from `IAPWS97.cs`), so a bug in the shipped backward
functions cannot be masked by the same bug in the oracle. Its values were
additionally cross-checked against the `iapws` Python package.

## Scope limits

These tests exercise `ThermoMath`, not `ThermoState.helper_add_p_insulated_vapor`
itself — that type pulls in the Unity object graph. The in-app acceptance
criteria (dial at 100%, liquid/two-phase non-regression, vapor-dome crossing)
still need a real Unity build or a human in the headset.
