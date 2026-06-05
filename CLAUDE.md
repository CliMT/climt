## graphify

This project has a knowledge graph at graphify-out/ with god nodes, community structure, and cross-file relationships.

Rules:
- For codebase questions, first run `graphify query "<question>"` when graphify-out/graph.json exists. Use `graphify path "<A>" "<B>"` for relationships and `graphify explain "<concept>"` for focused concepts. These return a scoped subgraph, usually much smaller than GRAPH_REPORT.md or raw grep output.
- For constructor signatures, kwargs, docstrings, and import paths of a specific symbol (including symbols imported from external packages like `sympl.set_backend`), run `python scripts/augment_graph.py show "<Symbol>"`. The augmenter populates these on every code node and resolves aliases through `__init__.py` re-exports.
- If graphify-out/wiki/index.md exists, use it for broad navigation instead of raw source browsing.
- Read graphify-out/GRAPH_REPORT.md only for broad architecture review or when query/path/explain do not surface enough context.
- After modifying code, run `graphify update .` followed by `python scripts/augment_graph.py` to keep the graph current. The augment step is required: graphify's AST extractor skips `.pyx` files and doesn't link the `docs/PLUG_AND_PLAY_ARCHITECTURE.md` concept nodes to their code counterparts. `scripts/augment_graph.py` adds Cython def nodes, Python-wrapper → `.pyx` `imports`/`calls` edges, `.pyx` → Fortran/C extern `calls` edges, and `concept → code` `aliases` edges. It is idempotent.
