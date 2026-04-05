# -*- coding: utf-8 -*-
"""
Script 2 (Extremes): Predict internal parameters of the Emanuel scheme using a NN.
Incorporates Magnitude-Based Importance Weighting to focus on rare extreme convective events.
"""
import jax
import jax.numpy as jnp
import flax.linen as nn
import optax
from datetime import timedelta
from climt import get_grid, get_default_state, EmanuelConvectionPythonV3
from climt._components.emanuel.pure_python_v3 import _jax_vectorized_convect, EmanuelParams, _bolton_q_sat_jax

jax.config.update("jax_enable_x64", True)

class ParameterPredictor(nn.Module):
    @nn.compact
    def __call__(self, x):
        x = x.reshape(-1)
        x = nn.Dense(64)(x)
        x = nn.relu(x)
        x = nn.Dense(64)(x)
        x = nn.relu(x)
        x = nn.Dense(2)(x)
        entp = jax.nn.softplus(x[0]) * 1.5
        omtrain = jax.nn.softplus(x[1]) * 50.0
        return entp, omtrain

def train():
    nlev = 30
    ncol = 1
    grid = get_grid(nx=ncol, ny=1, nz=nlev)
    conv = EmanuelConvectionPythonV3()
    delt = timedelta(minutes=15).total_seconds()
    
    model = ParameterPredictor()
    key = jax.random.PRNGKey(1)
    input_dim = 5 * nlev + (nlev + 1) + 1
    params = model.init(key, jnp.zeros((input_dim, ncol)))
    optimizer = optax.adam(1e-3)
    opt_state = optimizer.init(params)
    
    def get_default_params(entp, omtrain):
        return EmanuelParams(
            IPBL=0, MINORIG=1, ELCRIT=0.0011, TLCRIT=-55.0, ENTP=entp, SIGD=0.05, 
            SIGS=0.12, OMTRAIN=omtrain, OMTSNOW=5.5, COEFFR=1.0, COEFFS=0.8, CU=0.7,
            BETA=10.0, DTMAX=0.9, ALPHA=0.1, DAMP=0.1, CPD=1005.7, CPV=1870.0, 
            CL=2500.0, RV=461.5, RD=287.04, LV0=2.501e6, G=9.8, ROWL=1000.0, DELT0=300.0
        )

    def get_input_vec(t, q, u, v, p, ph, cbmf):
        return jnp.concatenate([t, q, u, v, p, ph, cbmf.reshape(1, ncol)], axis=0)
    
    def generate_data(t_surf, q_surf):
        temp = jnp.linspace(t_surf, 250, nlev).reshape(nlev, ncol)
        q = jnp.zeros((nlev, ncol)).at[0:10].set(q_surf)
        u, v, cbmf = jnp.zeros((nlev, ncol)), jnp.zeros((nlev, ncol)), jnp.zeros(ncol)
        state_sympl = get_default_state([conv], grid_state=grid)
        p = jnp.array(state_sympl['air_pressure'].values / 100.0).reshape(nlev, ncol)
        ph = jnp.array(state_sympl['air_pressure_on_interface_levels'].values / 100.0).reshape(nlev+1, ncol)
        qs = _bolton_q_sat_jax(temp, p * 100, 287.04, 461.5)
        
        # Target with "true" parameters
        true_p = get_default_params(1.5, 50.0)
        results = _jax_vectorized_convect(temp, q, qs, u, v, p, ph, nlev, nlev - 3, 0, delt, cbmf, true_p)
        target = jnp.concatenate([results[0], results[1], results[2], results[3]], axis=0)
        return temp, q, u, v, p, ph, cbmf, qs, target

    # --- MAGNITUDE-WEIGHTED LOSS FUNCTION ---
    @jax.jit
    def loss_fn(params, t, q, u, v, p, ph, cbmf, qs, target):
        input_vec = get_input_vec(t, q, u, v, p, ph, cbmf)
        pred_entp, pred_omtrain = model.apply(params, input_vec)
        
        p_struct = get_default_params(pred_entp, pred_omtrain)
        results = _jax_vectorized_convect(t, q, qs, u, v, p, ph, nlev, nlev - 3, 0, delt, cbmf, p_struct)
        pred_tendencies = jnp.concatenate([results[0], results[1], results[2], results[3]], axis=0)
        
        # 1. Identify Extremes: Calculate the magnitude of the target tendency
        # Tendencies are usually very small (units per second), so we scale them up to evaluate the exponent
        tendency_magnitude = jnp.mean(jnp.abs(target))
        
        # 2. Importance Weighting: Normal events get weight ~1.0, extremes get exponentially higher weights.
        # This forces the optimizer to prioritize matching the parameters of extreme convective events.
        weight = 1.0 + jnp.exp(5000.0 * tendency_magnitude)
        
        mse = jnp.mean((pred_tendencies - target)**2)
        return weight * mse

    @jax.jit
    def step(params, opt_state, t, q, u, v, p, ph, cbmf, qs, target):
        loss, grads = jax.value_and_grad(loss_fn)(params, t, q, u, v, p, ph, cbmf, qs, target)
        updates, opt_state = optimizer.update(grads, opt_state)
        return optax.apply_updates(params, updates), opt_state, loss

    current_params, current_opt_state = params, opt_state
    
    print("Starting Extremes Parameter Prediction demonstration...")
    for i in range(5):
        # Ramp up surface temperature and humidity to simulate increasingly extreme profiles
        t_surf = 300.0 + (i * 3.0) 
        q_surf = 0.015 + (i * 0.002)
        
        t, q, u, v, p, ph, cbmf, qs, target = generate_data(t_surf, q_surf)
        current_params, current_opt_state, loss = step(current_params, current_opt_state, t, q, u, v, p, ph, cbmf, qs, target)
        
        pred_entp, pred_omtrain = model.apply(current_params, get_input_vec(t, q, u, v, p, ph, cbmf))
        
        tendency_magnitude = jnp.mean(jnp.abs(target))
        print(f"Profile {i} | Target Mag: {tendency_magnitude:.6f} | Loss: {loss:.6f} | ENTP: {pred_entp:.4f}")

if __name__ == "__main__":
    train()