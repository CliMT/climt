"""Post-process graphify-out/graph.json so LLMs can trace the wrapper boundary.

Graphify's AST extractor skips .pyx files entirely and never links
docs/PLUG_AND_PLAY_ARCHITECTURE.md concept nodes to the matching code nodes.
This script adds the missing nodes and edges in-place. Idempotent: re-runs
detect existing additions via the ``augmented`` marker.

Run after every ``graphify update``.
"""

from __future__ import annotations

import ast
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GRAPH = ROOT / "graphify-out" / "graph.json"
PYX_GLOB = "climt/_components/**/*.pyx"
ARCH_DOC = "docs/PLUG_AND_PLAY_ARCHITECTURE.md"

DEF_RE = re.compile(r"^def\s+([A-Za-z_][A-Za-z0-9_]*)\s*\(", re.MULTILINE)
EXTERN_BLOCK_RE = re.compile(
    r"^cdef\s+extern[^\n]*:\n((?:[ \t]+.*\n)+)", re.MULTILINE
)
EXTERN_SYMBOL_RE = re.compile(
    r"^\s+(?:void|int|double|float|long|char)\s+\*?\s*([A-Za-z_][A-Za-z0-9_]*)\s*\(",
    re.MULTILINE,
)


def slug(path: str) -> str:
    return re.sub(r"[^a-z0-9_]", "_", path.lower())


def pyx_file_id(rel_path: str) -> str:
    return slug(rel_path)


def pyx_def_id(rel_path: str, name: str) -> str:
    return f"{slug(rel_path)}_{name.lower()}"


def python_module_id(rel_path: str) -> str:
    # Mirror graphify AST convention: <parent_dir>_<stem>
    p = Path(rel_path)
    return slug(f"{p.parent.name}_{p.stem}")


def load_graph() -> dict:
    return json.loads(GRAPH.read_text())


def save_graph(g: dict) -> None:
    GRAPH.write_text(json.dumps(g, indent=2, ensure_ascii=False))


def edges_list(g: dict):
    return g.get("links") if "links" in g else g.setdefault("edges", [])


def add_node(g, nodes_by_id, node):
    if node["id"] in nodes_by_id:
        return False
    node["augmented"] = True
    g["nodes"].append(node)
    nodes_by_id[node["id"]] = node
    return True


def add_edge(edges, seen, source, target, relation, confidence="EXTRACTED", score=1.0, source_file=""):
    key = (source, target, relation)
    if key in seen:
        return False
    edges.append({
        "source": source,
        "target": target,
        "relation": relation,
        "confidence": confidence,
        "confidence_score": score,
        "source_file": source_file,
        "weight": 1.0,
        "augmented": True,
    })
    seen.add(key)
    return True


def discover_pyx_files():
    return sorted(ROOT.glob(PYX_GLOB))


def parse_pyx(pyx_path: Path):
    text = pyx_path.read_text()
    defs = DEF_RE.findall(text)
    externs = []
    for block in EXTERN_BLOCK_RE.findall(text):
        externs.extend(EXTERN_SYMBOL_RE.findall(block))
    return defs, externs


def python_importers_of(pyx_path: Path) -> list[Path]:
    """Return .py siblings that do `from . import <pyx_stem>` (or `as ...`)."""
    stem = pyx_path.stem  # e.g. _emanuel_convection
    pat = re.compile(rf"from\s+\.\s+import\s+{re.escape(stem)}(?:\s+as\s+\w+)?\b")
    out = []
    for py in pyx_path.parent.glob("*.py"):
        try:
            if pat.search(py.read_text()):
                out.append(py)
        except OSError:
            continue
    return out


def fortran_node_by_label(nodes_by_label_lower, name: str):
    return nodes_by_label_lower.get(name.lower())


def _fmt_arg(a: ast.arg, default=None) -> str:
    s = a.arg
    if a.annotation is not None:
        try:
            s += f": {ast.unparse(a.annotation)}"
        except Exception:
            pass
    if default is not None:
        try:
            s += f"={ast.unparse(default)}"
        except Exception:
            s += "=..."
    return s


def signature_of(fn: ast.FunctionDef | ast.AsyncFunctionDef) -> str:
    a = fn.args
    parts: list[str] = []
    pos = list(a.posonlyargs) + list(a.args)
    pos_defaults = [None] * (len(pos) - len(a.defaults)) + list(a.defaults)
    for arg, d in zip(pos, pos_defaults):
        if arg.arg == "self":
            continue
        parts.append(_fmt_arg(arg, d))
    if a.vararg:
        parts.append("*" + _fmt_arg(a.vararg))
    elif a.kwonlyargs:
        parts.append("*")
    for arg, d in zip(a.kwonlyargs, a.kw_defaults):
        parts.append(_fmt_arg(arg, d))
    if a.kwarg:
        parts.append("**" + _fmt_arg(a.kwarg))
    return f"{fn.name}({', '.join(parts)})"


def first_doc_line(node) -> str | None:
    doc = ast.get_docstring(node)
    if not doc:
        return None
    for line in doc.splitlines():
        line = line.strip()
        if line:
            return line[:280]
    return None


def full_doc(node, max_chars: int = 4000) -> str | None:
    doc = ast.get_docstring(node)
    if not doc:
        return None
    doc = doc.strip()
    return doc[:max_chars] + ("..." if len(doc) > max_chars else "")


def base_names(cls_node: ast.ClassDef) -> list[str]:
    out = []
    for b in cls_node.bases:
        try:
            out.append(ast.unparse(b))
        except Exception:
            continue
    return out


SYMPL_PROPERTY_ATTRS = (
    "input_properties",
    "output_properties",
    "diagnostic_properties",
    "tendency_properties",
)


def extract_sympl_properties(cls_node: ast.ClassDef) -> dict:
    """Pull out sympl property dicts as a flat schema: {attr: {field: {dims, units}}}."""
    out: dict[str, dict] = {}
    for stmt in cls_node.body:
        if not isinstance(stmt, ast.Assign):
            continue
        if len(stmt.targets) != 1 or not isinstance(stmt.targets[0], ast.Name):
            continue
        name = stmt.targets[0].id
        if name not in SYMPL_PROPERTY_ATTRS:
            continue
        if not isinstance(stmt.value, ast.Dict):
            continue
        attr: dict[str, dict] = {}
        for k, v in zip(stmt.value.keys, stmt.value.values):
            if not isinstance(k, ast.Constant) or not isinstance(v, ast.Dict):
                continue
            inner: dict[str, str] = {}
            for ik, iv in zip(v.keys, v.values):
                if not isinstance(ik, ast.Constant):
                    continue
                try:
                    inner[str(ik.value)] = ast.literal_eval(iv)
                except Exception:
                    try:
                        inner[str(ik.value)] = ast.unparse(iv)
                    except Exception:
                        continue
            attr[str(k.value)] = inner
        out[name] = attr
    return out


def collect_python_signatures():
    """Walk every .py under climt/, scripts/, examples/, tests/ and
    return {(source_file_rel, symbol_label) -> {signature, doc, imports, ...}}.
    """
    out: dict[tuple[str, str], dict] = {}
    import_map: dict[str, dict[str, str]] = {}  # rel_file -> {alias: source_module}
    roots = [ROOT / d for d in ("climt", "scripts", "examples", "tests")]
    for r in roots:
        for py in r.rglob("*.py"):
            try:
                src = py.read_text()
                tree = ast.parse(src)
            except (SyntaxError, OSError, UnicodeDecodeError):
                continue
            rel = str(py.relative_to(ROOT))

            imports: dict[str, str] = {}
            for node in tree.body:
                if isinstance(node, ast.ImportFrom) and node.module:
                    for n in node.names:
                        imports[n.asname or n.name] = f"from {node.module} import {n.name}"
                elif isinstance(node, ast.Import):
                    for n in node.names:
                        imports[n.asname or n.name.split(".")[0]] = f"import {n.name}"
            if imports:
                import_map[rel] = imports

            for node in ast.walk(tree):
                if isinstance(node, ast.ClassDef):
                    init = next(
                        (b for b in node.body
                         if isinstance(b, (ast.FunctionDef, ast.AsyncFunctionDef))
                         and b.name == "__init__"),
                        None,
                    )
                    sig = f"{node.name}{signature_of(init)[len('__init__'):]}" if init else f"{node.name}()"
                    doc = first_doc_line(node) or (first_doc_line(init) if init else None)
                    props = extract_sympl_properties(node)
                    out[(rel, node.name)] = {
                        "signature": sig,
                        "docstring": doc,
                        "full_docstring": full_doc(node) or (full_doc(init) if init else None),
                        "bases": base_names(node) or None,
                        "kind": "class",
                        "sympl_properties": props or None,
                    }
                elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    out[(rel, node.name)] = {
                        "signature": signature_of(node),
                        "docstring": first_doc_line(node),
                        "full_docstring": full_doc(node),
                        "kind": "function",
                    }
    return out, import_map


def main():
    g = load_graph()
    # Strip previously-augmented additions so logic changes take effect.
    g["nodes"] = [n for n in g["nodes"] if not n.get("augmented")]
    edges_key = "links" if "links" in g else "edges"
    g[edges_key] = [e for e in g[edges_key] if not e.get("augmented")]
    nodes_by_id = {n["id"]: n for n in g["nodes"]}
    # Label -> node, prefer fortran/.f90, then any code
    label_to_code = {}
    for n in g["nodes"]:
        if n.get("file_type") != "code":
            continue
        lbl = (n.get("label") or "").rstrip("()").lower()
        if not lbl:
            continue
        sf = n.get("source_file") or ""
        # Prefer fortran sources for extern symbol matching
        existing = label_to_code.get(lbl)
        if existing is None or (sf.endswith((".f90", ".F90")) and not (existing.get("source_file") or "").endswith((".f90", ".F90"))):
            label_to_code[lbl] = n

    edges = edges_list(g)
    seen_edges = {(e["source"], e["target"], e.get("relation")) for e in edges}

    n_added_nodes = 0
    n_added_edges = 0

    # ---- 1. Pyx file + def nodes, wrapper imports, extern bridges ----
    for pyx in discover_pyx_files():
        rel = str(pyx.relative_to(ROOT))
        defs, externs = parse_pyx(pyx)

        file_id = pyx_file_id(rel)
        if add_node(g, nodes_by_id, {
            "id": file_id,
            "label": pyx.name,
            "file_type": "code",
            "source_file": rel,
            "source_location": "L1",
        }):
            n_added_nodes += 1

        # Def nodes
        def_ids = []
        for name in defs:
            did = pyx_def_id(rel, name)
            def_ids.append(did)
            if add_node(g, nodes_by_id, {
                "id": did,
                "label": f"{name}()",
                "file_type": "code",
                "source_file": rel,
                "source_location": None,
            }):
                n_added_nodes += 1
            if add_edge(edges, seen_edges, did, file_id, "contains", source_file=rel):
                n_added_edges += 1

        # Wrapper Python module -> pyx file (+ each def).
        # Graphify keeps multiple duplicate IDs for the same .py file under
        # different normalizations; link the bridge to every variant we find.
        for py in python_importers_of(pyx):
            py_rel = str(py.relative_to(ROOT))
            full = slug(py_rel)
            collapsed = re.sub(r"_+", "_", full)
            candidates = {full, collapsed, python_module_id(py_rel)}
            targets = [c for c in candidates if c in nodes_by_id]
            if not targets:
                stub = full
                if add_node(g, nodes_by_id, {
                    "id": stub,
                    "label": py.name,
                    "file_type": "code",
                    "source_file": py_rel,
                    "source_location": "L1",
                }):
                    n_added_nodes += 1
                targets = [stub]
            for target_id in targets:
                if add_edge(edges, seen_edges, target_id, file_id, "imports", source_file=py_rel):
                    n_added_edges += 1
                for did in def_ids:
                    if add_edge(edges, seen_edges, target_id, did, "calls",
                                confidence="INFERRED", score=0.85, source_file=py_rel):
                        n_added_edges += 1

        # Pyx def -> fortran/C extern symbol (best-effort)
        # Heuristic: any def named `init_*`, `set_*`, `do_*`, `convect`, `get_new_state`
        # routes to all externs declared in the file.
        if externs:
            for did in def_ids:
                for sym in externs:
                    fnode = label_to_code.get(sym.lower())
                    if fnode is None:
                        # Create a stub extern node so the boundary is visible
                        stub_id = f"{file_id}_extern_{sym.lower()}"
                        if add_node(g, nodes_by_id, {
                            "id": stub_id,
                            "label": f"{sym}()",
                            "file_type": "code",
                            "source_file": rel,
                            "source_location": None,
                            "extern": True,
                        }):
                            n_added_nodes += 1
                        target = stub_id
                    else:
                        target = fnode["id"]
                    if add_edge(edges, seen_edges, did, target, "calls",
                                confidence="INFERRED", score=0.75, source_file=rel):
                        n_added_edges += 1

    # ---- 2a. Build canonical class -> module-file map from public API ----
    api_init = ROOT / "climt" / "_components" / "__init__.py"
    class_to_module_file: dict[str, str] = {}
    if api_init.exists():
        text = api_init.read_text()
        # Normalize parenthesised imports onto a single logical line.
        text = re.sub(
            r"from\s+\.(\w+)\s+import\s*\(([^)]*)\)",
            lambda m: f"from .{m.group(1)} import " + ",".join(
                p.strip() for p in m.group(2).split(",") if p.strip()
            ),
            text,
        )
        line_re = re.compile(r"^from\s+\.(\w+)\s+import\s+(.+)$", re.MULTILINE)
        for mod, names in line_re.findall(text):
            for raw in names.split(","):
                cls = raw.strip().rstrip("\\").strip()
                if not cls or cls.startswith("#"):
                    continue
                # Module may be a file (mod.py) or a package (mod/__init__.py)
                candidates = [
                    f"climt/_components/{mod}.py",
                    f"climt/_components/{mod}/__init__.py",
                    f"climt/_components/{mod}/component.py",
                ]
                for c in candidates:
                    if (ROOT / c).exists():
                        class_to_module_file[cls] = c
                        break

    file_to_node_ids: dict[str, list[str]] = {}
    file_any_node: dict[str, list[str]] = {}
    for n in g["nodes"]:
        sf = n.get("source_file") or ""
        if not sf:
            continue
        file_any_node.setdefault(sf, []).append(n["id"])
        if n.get("source_location") in (None, "L1", ""):
            file_to_node_ids.setdefault(sf, []).append(n["id"])

    # ---- 2b. Doc-concept -> code-node alias edges ----
    for n in g["nodes"]:
        if n.get("file_type") != "concept":
            continue
        sf = n.get("source_file") or ""
        if ARCH_DOC not in sf and "graphify-out/memory/" not in sf:
            continue
        raw_label = (n.get("label") or "").rstrip("()")
        lbl = raw_label.lower()
        targets: list[str] = []
        code = label_to_code.get(lbl)
        if code is not None and code["id"] != n["id"]:
            targets.append(code["id"])
        # Fallback: resolve via public API registry to the module file
        mod_file = class_to_module_file.get(raw_label)
        if mod_file:
            mod_targets = file_to_node_ids.get(mod_file) or file_any_node.get(mod_file, [])
            for tid in mod_targets:
                if tid != n["id"]:
                    targets.append(tid)
        for tid in dict.fromkeys(targets):  # dedupe, preserve order
            if add_edge(edges, seen_edges, n["id"], tid, "aliases",
                        confidence="INFERRED", score=0.95, source_file=sf):
                n_added_edges += 1

    # ---- 3. Enrich Python code nodes with signature + docstring ----
    sigs, import_map = collect_python_signatures()

    # 3a. Create missing class/function nodes that graphify's AST step dropped.
    # Index existing (file, label) pairs to detect gaps.
    existing_by_file_label = {
        (n.get("source_file") or "", (n.get("label") or "").rstrip("()"))
        for n in g["nodes"]
    }
    for (sf, name), info in sigs.items():
        if info.get("kind") != "class":
            continue
        if (sf, name) in existing_by_file_label:
            continue
        node_id = f"{slug(sf)}_{name.lower()}"
        if node_id in nodes_by_id:
            continue
        if add_node(g, nodes_by_id, {
            "id": node_id,
            "label": name,
            "file_type": "code",
            "source_file": sf,
            "source_location": None,
            "signature": info.get("signature"),
            "docstring": info.get("docstring"),
            "full_docstring": info.get("full_docstring"),
            "bases": info.get("bases"),
            "symbol_kind": "class",
            "sympl_properties": info.get("sympl_properties"),
        }):
            n_added_nodes += 1
            # Link to the file-level node if present so traversal can find it
            file_node = next(
                (nid for nid, nn in nodes_by_id.items()
                 if nn.get("source_file") == sf
                 and nn.get("source_location") in (None, "L1", "")
                 and (nn.get("label") or "").endswith(".py")),
                None,
            )
            if file_node:
                if add_edge(edges, seen_edges, node_id, file_node, "contains", source_file=sf):
                    n_added_edges += 1

    n_enriched = 0
    for n in g["nodes"]:
        sf = n.get("source_file") or ""
        if not sf.endswith(".py"):
            continue
        # Strip trailing () that graphify adds to function labels.
        label = (n.get("label") or "").rstrip("()")
        if not label or label.endswith(".py"):
            continue
        info = sigs.get((sf, label))
        if info is None:
            continue
        changed = False
        if info.get("signature") and n.get("signature") != info["signature"]:
            n["signature"] = info["signature"]; changed = True
        if info.get("docstring") and n.get("docstring") != info["docstring"]:
            n["docstring"] = info["docstring"]; changed = True
        if info.get("kind") and n.get("symbol_kind") != info["kind"]:
            n["symbol_kind"] = info["kind"]; changed = True
        if info.get("sympl_properties") and n.get("sympl_properties") != info["sympl_properties"]:
            n["sympl_properties"] = info["sympl_properties"]; changed = True
        if info.get("full_docstring") and n.get("full_docstring") != info["full_docstring"]:
            n["full_docstring"] = info["full_docstring"]; changed = True
        if info.get("bases") and n.get("bases") != info["bases"]:
            n["bases"] = info["bases"]; changed = True
        if changed:
            n["augmented_meta"] = True
            n_enriched += 1

    # Attach imports table to module-file nodes (so set_backend's origin is findable).
    for n in g["nodes"]:
        sf = n.get("source_file") or ""
        if not sf.endswith(".py"):
            continue
        if n.get("source_location") not in (None, "L1", ""):
            continue
        imports = import_map.get(sf)
        if imports and n.get("imports") != imports:
            n["imports"] = imports
            n["augmented_meta"] = True
            n_enriched += 1

    save_graph(g)
    print(f"augment_graph: +{n_added_nodes} nodes, +{n_added_edges} edges, "
          f"{n_enriched} python nodes enriched with signature/docstring/imports")


def show(symbol: str) -> int:
    """Print signature + docstring + import path for a symbol.

    Handles three cases:
      1. Symbol has a code node -> print its signature/docstring.
      2. Symbol is a concept node aliased to code -> follow the alias.
      3. Symbol has no node but is imported elsewhere -> resolve from import map +
         re-parse the source file for the actual signature.
    """
    g = load_graph()
    edges_key = "links" if "links" in g else "edges"
    target = symbol.rstrip("()").lower()
    by_id = {n["id"]: n for n in g["nodes"]}

    hits = [n for n in g["nodes"]
            if (n.get("label") or "").rstrip("()").lower() == target
            and (n.get("label") or "") not in ("__init__.py", "component.py")
            and not (n.get("label") or "").endswith(".py")]
    # If a hit is itself a re-exporting file (label.endswith(.py)), follow its
    # imports table to the real definition file.

    # Follow aliases: if any hit is a concept with an aliases edge to a code node,
    # add the code node to the display set. If that code node is itself a
    # re-exporting __init__.py / component.py, walk its imports table to find
    # where the symbol is actually defined and look up that node instead.
    sigs_cache, _ = collect_python_signatures()
    for e in g[edges_key]:
        if e.get("relation") != "aliases":
            continue
        src = by_id.get(e["source"]); tgt = by_id.get(e["target"])
        if not (src and (src.get("label") or "").rstrip("()").lower() == target):
            continue
        if tgt is None:
            continue
        label = tgt.get("label") or ""
        sf = tgt.get("source_file") or ""
        if label.endswith(".py"):
            # Walk re-exports: __init__.py → component.py → class def
            visited = set()
            current_file = sf
            while current_file and current_file not in visited:
                visited.add(current_file)
                imp = next((n.get("imports") for n in g["nodes"]
                            if n.get("source_file") == current_file and n.get("imports")), None)
                if not imp or symbol.rstrip("()") not in imp:
                    break
                stmt = imp[symbol.rstrip("()")]
                # `from .component import Foo` -> sibling component.py
                m = re.match(r"from\s+(\.?[\w\.]+)\s+import\s+", stmt)
                if not m:
                    break
                mod = m.group(1).lstrip(".")
                base = str(Path(current_file).parent)
                next_file = f"{base}/{mod.replace('.', '/')}.py"
                if not (ROOT / next_file).exists():
                    next_file = f"{base}/{mod.replace('.', '/')}/__init__.py"
                if not (ROOT / next_file).exists():
                    break
                # Found the file — look for a class/function node there
                resolved = next(
                    (n for n in g["nodes"]
                     if n.get("source_file") == next_file
                     and (n.get("label") or "").rstrip("()").lower() == target),
                    None,
                )
                if resolved is not None and resolved not in hits:
                    hits.append(resolved)
                    # Re-enrich on the fly in case the node lacks signature
                    info = sigs_cache.get((next_file, symbol.rstrip("()")))
                    if info:
                        resolved.setdefault("signature", info.get("signature"))
                        resolved.setdefault("docstring", info.get("docstring"))
                    break
                current_file = next_file
        elif tgt not in hits:
            hits.append(tgt)

    # Fallback: no node at all -> search import map and re-parse the source.
    if not hits:
        importers: dict[str, list[str]] = {}
        for n in g["nodes"]:
            imp = n.get("imports") or {}
            for alias, stmt in imp.items():
                if alias == symbol or alias == symbol.rstrip("()"):
                    importers.setdefault(stmt, []).append(n.get("source_file"))
        if not importers:
            print(f"No node, alias, or import found for {symbol!r}")
            return 1
        print(f"\n{symbol}  (no node — resolved via imports)")
        for stmt, files in importers.items():
            print(f"  {stmt}   ({len(files)} file{'s' if len(files)!=1 else ''}, e.g. {files[0]})")
        return 0

    for n in hits:
        sf = n.get("source_file") or "?"
        loc = n.get("source_location") or ""
        print(f"\n{n.get('label')}  [{sf}{(' ' + loc) if loc else ''}]")
        if n.get("signature"):
            print(f"  signature: {n['signature']}")
        if n.get("docstring"):
            print(f"  doc:       {n['docstring']}")
        if n.get("symbol_kind"):
            print(f"  kind:      {n['symbol_kind']}")
        if n.get("bases"):
            print(f"  bases:     {', '.join(n['bases'])}")
        if n.get("full_docstring") and n["full_docstring"] != n.get("docstring"):
            print("  docstring:")
            for line in n["full_docstring"].splitlines():
                print(f"      {line}")
        props = n.get("sympl_properties") or {}
        for attr in ("input_properties", "output_properties",
                     "diagnostic_properties", "tendency_properties"):
            fields = props.get(attr)
            if not fields:
                continue
            print(f"  {attr}:")
            for fname, meta in fields.items():
                units = meta.get("units", "?")
                dims = meta.get("dims", "?")
                print(f"      {fname:50s} units={units!s:18s} dims={dims}")
    # Show how this symbol is imported elsewhere
    edges_key = "links" if "links" in g else "edges"
    importers = {}
    for n in g["nodes"]:
        imp = n.get("imports") or {}
        if symbol in imp or symbol.rstrip("()") in imp:
            key = symbol if symbol in imp else symbol.rstrip("()")
            importers.setdefault(imp[key], []).append(n.get("source_file"))
    if importers:
        print("\n  imported as:")
        for stmt, files in importers.items():
            print(f"    {stmt}   ({len(files)} file{'s' if len(files)!=1 else ''}, e.g. {files[0]})")
    return 0


if __name__ == "__main__":
    import sys
    if len(sys.argv) >= 3 and sys.argv[1] == "show":
        raise SystemExit(show(sys.argv[2]))
    main()
