# -*- coding: utf-8 -*-
"""
Script 2: Predict internal parameters of the Emanuel scheme using a NN.
The NN takes the state as input and predicts parameters like 'entrainment_coeff'.
Loss is computed based on the difference between the scheme's output and observations.
"""
import jax
import jax.numpy as jnp
import flax.linen as nn
import optax
from datetime import timedelta
from climt import get_grid, get_default_state, EmanuelConvectionPythonV3
from climt._components.emanuel.pure_python_v3 import _jax_vectorized_convect, EmanuelParams, _bolton_q_sat_jax

# 1. Setup JAX
jax.config.update("jax_enable_x64", True)

# 2. Define the Parameter Predictor (MLP)
class ParameterPredictor(nn.Module):
    @nn.compact
    def __call__(self, x):
        # x is the concatenated state vector per column
        # shape: (nlev*5 + nlev+1 + 1, ncol) -> Flatten
        x = x.reshape(-1)
        x = nn.Dense(64)(x)
        x = nn.relu(x)
        x = nn.Dense(64)(x)
        x = nn.relu(x)
        # Predict 2 parameters: ENTP and OMTRAIN
        x = nn.Dense(2)(x)
        # Use softplus for a smooth, physically bounded output
        entp = jax.nn.softplus(x[0]) * 1.5  # scale to typical range ~1.5
        omtrain = jax.nn.softplus(x[1]) * 50.0 # scale to typical range ~50.0
        return entp, omtrain

# 3. Training Logic
def train():
    nlev = 30
    ncol = 1
    grid = get_grid(nx=ncol, ny=1, nz=nlev)
    conv = EmanuelConvectionPythonV3()
    timestep = timedelta(minutes=15)
    delt = timestep.total_seconds()
    
    # Initialize Predictor
    model = ParameterPredictor()
    key = jax.random.PRNGKey(1)
    
    # Inputs: t, q, u, v, p (all nlev), ph (nlev+1), cbmf (1)
    input_dim = 5 * nlev + (nlev + 1) + 1
    dummy_input = jnp.zeros((input_dim, ncol))
    params = model.init(key, dummy_input)
    
    # Optimizer
    optimizer = optax.adam(1e-3)
    opt_state = optimizer.init(params)
    
    # Helper to get default params
    def get_default_params(entp, omtrain):
        # We replace the two we are predicting
        return EmanuelParams(
            IPBL=0, MINORIG=1, ELCRIT=0.0011, TLCRIT=-55.0,
            ENTP=entp, SIGD=0.05, SIGS=0.12, OMTRAIN=omtrain,
            OMTSNOW=5.5, COEFFR=1.0, COEFFS=0.8, CU=0.7,
            BETA=10.0, DTMAX=0.9, ALPHA=0.1, DAMP=0.1,
            CPD=1005.7, CPV=1870.0, CL=2500.0, RV=461.5,
            RD=287.04, LV0=2.501e6, G=9.8, ROWL=1000.0, DELT0=300.0
        )

    def get_input_vec(t_vals, q_vals, u_vals, v_vals, p_vals, ph_vals, cbmf_vals):
        cbmf_2d = cbmf_vals.reshape(1, ncol)
        return jnp.concatenate([t_vals, q_vals, u_vals, v_vals, p_vals, ph_vals, cbmf_2d], axis=0)

    # 4. Mock "Observed" Tendency
    true_entp = 1.5
    true_omtrain = 50.0
    
    def generate_data(t_surf, q_surf):
        temp = jnp.linspace(t_surf, 250, nlev).reshape(nlev, ncol)
        q = jnp.zeros((nlev, ncol)).at[0:10].set(q_surf)
        u = jnp.zeros((nlev, ncol))
        v = jnp.zeros((nlev, ncol))
        cbmf = jnp.zeros(ncol)
        state_sympl = get_default_state([conv], grid_state=grid)
        p = jnp.array(state_sympl['air_pressure'].values / 100.0).reshape(nlev, ncol)
        ph = jnp.array(state_sympl['air_pressure_on_interface_levels'].values / 100.0).reshape(nlev+1, ncol)
        
        qs = _bolton_q_sat_jax(temp, p * 100, 287.04, 461.5)
        
        true_p = get_default_params(true_entp, true_omtrain)
        results = _jax_vectorized_convect(
            temp, q, qs, u, v,
            p, ph, nlev, nlev - 3, 0, delt, cbmf, true_p
        )
        
        # Suppose our observations include tendencies for t, q, u, v
        target_tendencies = jnp.concatenate([results[0], results[1], results[2], results[3]], axis=0)
        
        return temp, q, u, v, p, ph, cbmf, qs, target_tendencies

    # Loss Function
    @jax.jit
    def loss_fn(params, t, q, u, v, p, ph, cbmf, qs, target):
        input_vec = get_input_vec(t, q, u, v, p, ph, cbmf)
        
        # 1. Predict parameters
        pred_entp, pred_omtrain = model.apply(params, input_vec)
        
        # 2. Run Emanuel with predicted parameters
        p_struct = get_default_params(pred_entp, pred_omtrain)
        results = _jax_vectorized_convect(
            t, q, qs, u, v,
            p, ph, nlev, nlev - 3, 0, delt, cbmf, p_struct
        )
        
        pred_tendencies = jnp.concatenate([results[0], results[1], results[2], results[3]], axis=0)
        
        # Loss: Difference in output tendencies
        return jnp.mean((pred_tendencies - target)**2)

    # Training step
    @jax.jit
    def step(params, opt_state, t, q, u, v, p, ph, cbmf, qs, target):
        loss, grads = jax.value_and_grad(loss_fn)(params, t, q, u, v, p, ph, cbmf, qs, target)
        updates, opt_state = optimizer.update(grads, opt_state)
        params = optax.apply_updates(params, updates)
        return params, opt_state, loss

    print("Starting Parameter Prediction demonstration...")
    current_params = params
    current_opt_state = opt_state
    
    for i in range(5):
        t_surf = 305.0 + i
        q_surf = 0.018 - 0.001 * i
        t, q, u, v, p, ph, cbmf, qs, target = generate_data(t_surf, q_surf)
        
        current_params, current_opt_state, loss = step(
            current_params, current_opt_state, t, q, u, v, p, ph, cbmf, qs, target
        )
        
        # Check what parameters are being predicted
        input_vec = get_input_vec(t, q, u, v, p, ph, cbmf)
        pred_entp, pred_omtrain = model.apply(current_params, input_vec)
        print(f"Epoch {i}, Loss: {loss:.6f}, Pred ENTP: {pred_entp:.4f}, Pred OMTRAIN: {pred_omtrain:.4f}")

if __name__ == "__main__":
    train()
