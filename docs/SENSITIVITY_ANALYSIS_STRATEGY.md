# Strategy: Differentiable Sensitivity Analysis of the Emanuel Scheme

This document outlines the architectural and mathematical strategy for leveraging the JAX-compiled, fully differentiable Emanuel Convection Scheme to conduct exact sensitivity analyses. 

Instead of relying on computationally expensive and mathematically approximate finite-difference methods (e.g., perturbing a parameter by $\epsilon$ and re-running the model), this strategy uses Automatic Differentiation (AutoDiff) to compute exact gradients and Jacobians. This allows us to rigorously quantify how different physical parameters control convective tendencies and rainfall across diverse climate states.

---

## 1. The Mathematical Framework

The Emanuel convection scheme can be viewed as a non-linear vector-valued function $F$:
$$ Y = F(X; \theta) $$
Where:
*   $X$ is the input atmospheric column state (Temperature $T$, Moisture $Q$, Winds $U, V$, Pressure $P$).
*   $\theta$ is the vector of internal scheme parameters (e.g., Entrainment `ENTP`, Fall speeds `OMTRAIN`, Downdraft fraction `SIGD`).
*   $Y$ is the output convective tendencies ($\Delta T, \Delta Q, \Delta U, \Delta V$) and diagnostic rainfall.

### Exact Sensitivities via AutoDiff
Because $F$ is written in pure JAX, we can analytically compute the derivative of any output with respect to any input parameter.

1.  **Scalar Targets (Gradients):** If our target is a scalar metric, such as total column rainfall $R$, we compute the gradient vector $\nabla_\theta R$:
    $$ \nabla_\theta R = \left[ \frac{\partial R}{\partial \text{ENTP}}, \frac{\partial R}{\partial \text{OMTRAIN}}, \dots \right] $$
    *Implementation:* `jax.grad(rainfall_fn, argnums=...)`

2.  **Vector Targets (Jacobians):** If our target is a vertical profile, such as the heating tendency $\Delta T$ at all $N$ model levels, we compute the $N \times K$ Jacobian matrix $J$:
    $$ J_{i,j} = \frac{\partial (\Delta T_i)}{\partial \theta_j} $$
    *Implementation:* `jax.jacrev(tendency_fn, argnums=...)`

---

## 2. Defining Representative Climate States

Sensitivities $\frac{\partial Y}{\partial \theta}$ are strongly state-dependent. A parameter that heavily controls tropical rainfall might be entirely inert in the subtropics. We must define distinct regimes $X_{regime}$ to evaluate the scheme.

### 1. Deep Tropics (e.g., ITCZ, Warm Pool)
*   **Characteristics:** High CAPE, high lower-tropospheric moisture, weak wind shear, deep conditionally unstable layers.
*   **Expected Behavior:** Deep, undiluted plumes. High sensitivity to entrainment (`ENTP`) and autoconversion (`ELCRIT`), as these dictate whether plumes survive to the tropopause and how much liquid water is unloaded.

### 2. Subtropics (e.g., Trade-Wind Regimes, Stratocumulus)
*   **Characteristics:** Strong trade inversion, dry free troposphere, subsidence, moderate boundary layer moisture.
*   **Expected Behavior:** Shallow convection. High sensitivity to mixing parameters and evaporation coefficients (`COEFFR`), as rain often evaporates in the dry sub-cloud layer before reaching the surface.

### 3. Midlatitudes (e.g., Storm Tracks)
*   **Characteristics:** Strong vertical wind shear, strong synoptic forcing, moderate CAPE, distinct baroclinic zones.
*   **Expected Behavior:** Organized convection. Higher sensitivity to downdraft parameters (`SIGD`, `BETA`) which drive cold-pool dynamics and convective momentum transport (`CU`).

*Data Source:* These profiles can be extracted from reanalysis datasets (ERA5), idealized Radiative-Convective Equilibrium (RCE) models, or Cloud-Resolving Models (CRMs).

---

## 3. Software Architecture and Implementation

The implementation strategy requires wrapping the JAX Emanuel kernel in target-specific closures that isolate the parameters of interest.

### Step 1: Parameter Injection Wrapper
Create a pure function that reconstructs the `EmanuelParams` struct from a subset of floating-point arguments, allowing JAX to trace them.

```python
def target_physics_wrapper(entp, omtrain, sigd, state_X):
    # Rebuild parameter struct with tracing
    params = get_default_params(ENTP=entp, OMTRAIN=omtrain, SIGD=sigd)
    # Run scheme
    tendencies = _jax_vectorized_convect(state_X.T, state_X.Q, ..., params)
    return tendencies
```

### Step 2: Metric Extraction
Define the specific outputs to differentiate.
*   **For Rainfall:** Extract $\Delta Q$, integrate over pressure levels, and return a scalar.
*   **For Heating Profiles:** Extract $\Delta T$ and return the vector.

### Step 3: Differentiation and Vectorization
Use JAX transformations to compute the sensitivities across thousands of profiles simultaneously.

```python
# 1. Define the Jacobian function with respect to parameters (args 0, 1, 2)
compute_jacobian = jax.jacrev(heating_profile_fn, argnums=(0, 1, 2))

# 2. Vectorize over a batch of diverse climate profiles (batch dimension 3)
batch_jacobian = jax.vmap(compute_jacobian, in_axes=(None, None, None, 0))

# 3. Execute on GPU/TPU
J_entp, J_omtrain, J_sigd = batch_jacobian(1.5, 50.0, 0.05, batch_of_states)
```

---

## 4. Scientific Applications and Outcomes

By computing these exact Jacobians across global climate states, several critical analyses become possible:

### 1. Parameter Identifiability (Ill-posedness)
If $\frac{\partial Y}{\partial \theta} \approx 0$ for a given parameter across a specific climate regime, that parameter is "unidentifiable" in that regime. 
*   *Application:* This tells ML practitioners that a Neural Network cannot safely predict that parameter from data in that regime, preventing catastrophic overfitting to noise.

### 2. Regime Shift Detection
By plotting the magnitude of a Jacobian term (e.g., $|\frac{\partial R}{\partial \text{ENTP}}|$) along a continuous geographic transect (e.g., from the Equator to the Subtropics), sudden spikes or drops mathematically define the boundary where the scheme transitions from deep to shallow convection regimes.

### 3. Targeted Tuning and Data Assimilation
When tuning GCMs to match observations, Data Assimilation (DA) systems (like Ensemble Kalman Filters) require the covariance between parameters and outputs.
*   *Application:* Instead of estimating this covariance statistically, the JAX scheme provides the exact mathematical tangent-linear model, drastically accelerating GCM tuning and reducing parameter uncertainty.

### 4. Climate Change Sensitivity
How will precipitation respond to warming? By evaluating the mixed second-derivative $\frac{\partial}{\partial T_{surf}} \left( \frac{\partial R}{\partial \theta} \right)$, we can quantify whether a parameter's influence amplifies or dampens in a warmer climate ($+4K$), exposing structural structural biases in the parameterization's formulation.