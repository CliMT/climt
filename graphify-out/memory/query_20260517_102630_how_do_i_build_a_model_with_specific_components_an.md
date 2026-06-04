---
type: "query"
date: "2026-05-17T10:26:30.838416+00:00"
question: "How do I build a model with specific components and backends in CliMT?"
contributor: "graphify"
source_nodes: ["RRTMGLongwave", "RRTMGShortwave", "EmanuelConvection", "SlabSurface", "ComponentBase3D", "ComponentBaseColumn", "load_atmospheric_properties()", "UnytBackend", "get_grid()", "get_default_state()"]
---

# Q: How do I build a model with specific components and backends in CliMT?

## Answer

To build a model in CliMT, follow the plug-and-play architecture defined in docs/PLUG_AND_PLAY_ARCHITECTURE.md.
1. Set the backend using sympl.set_backend(climt.UnytBackend()) (or DataArrayBackend).
2. Initialize independent components from climt._components (e.g., RRTMGLongwave, RRTMGShortwave, EmanuelConvection, SlabSurface).
3. Create the grid using climt.get_grid(nx, ny, nz).
4. Create the model state using climt.get_default_state(components, grid_state=grid).
5. Wrap the tendency components in a time stepper like sympl.AdamsBashforth(components).
6. In the simulation loop, pass the state to diagnostic components (like SimplePhysics) and then pass the state and timestep to the time stepper. Update the state with the returned diagnostics and new_state.

## Source Nodes

- RRTMGLongwave
- RRTMGShortwave
- EmanuelConvection
- SlabSurface
- ComponentBase3D
- ComponentBaseColumn
- load_atmospheric_properties()
- UnytBackend
- get_grid()
- get_default_state()