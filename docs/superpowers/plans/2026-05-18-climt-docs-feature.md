# climt.docs Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Ship a `python -m climt.docs show <Symbol>` CLI inside the climt package that surfaces constructor signatures, full docstrings, base classes, sympl `input/output/diagnostic/tendency_properties` (units + dims), and import paths for any public climt symbol — without requiring graphify or any external tooling.

**Architecture:** Pure-AST extraction over the `climt/` source tree (no import of climt at runtime — avoids Fortran/Cython compile dependency). A small `climt/_docs/` package holds the extractor, an in-memory registry built lazily on first call, a plain-text formatter, and a `__main__` CLI. Resolution follows three paths in order: direct name match → alias chain through `climt/_components/__init__.py` re-exports → imports-table fallback for symbols defined outside climt (e.g. `sympl.set_backend`). Results cache to `climt/_docs/_cache.json` keyed by source-tree mtime hash so repeat calls are O(load JSON).

**Tech Stack:** Python `ast` stdlib only. No new runtime dependencies. Tests use pytest (already in dev deps).

---

## File Structure

```
climt/_docs/
├── __init__.py      # Public API: show(symbol) -> str, list_symbols() -> list[str]
├── __main__.py      # `python -m climt.docs ...` dispatch
├── extract.py       # AST walker: source file -> list[SymbolInfo]
├── registry.py      # Build + cache full registry, resolve a query string
├── format.py        # SymbolInfo -> plain-text rendering
└── _cache.json      # Generated; not in git

tests/_docs/
├── __init__.py
├── test_extract.py
├── test_registry.py
├── test_format.py
└── test_cli.py
```

Each file has one responsibility. `extract.py` knows AST and nothing else. `registry.py` knows how to walk the climt tree and resolve aliases. `format.py` knows plain-text rendering. `__main__.py` is a thin CLI dispatch. They communicate via a single `SymbolInfo` dataclass defined in `extract.py`.

---

### Task 1: SymbolInfo dataclass + failing test

**Files:**
- Create: `climt/_docs/__init__.py` (empty placeholder so import works)
- Create: `climt/_docs/extract.py`
- Create: `tests/_docs/__init__.py` (empty)
- Create: `tests/_docs/test_extract.py`

- [ ] **Step 1: Write the failing test**

`tests/_docs/test_extract.py`:
```python
from climt._docs.extract import SymbolInfo


def test_symbolinfo_holds_required_fields():
    s = SymbolInfo(
        name="Foo",
        kind="class",
        source_file="climt/x.py",
        source_line=10,
        signature="Foo(a=1)",
        docstring="One liner.",
        full_docstring="One liner.\n\nDetails.",
        bases=["Bar"],
        sympl_properties={"input_properties": {"x": {"units": "K", "dims": ["*"]}}},
    )
    assert s.name == "Foo"
    assert s.kind == "class"
    assert s.sympl_properties["input_properties"]["x"]["units"] == "K"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/_docs/test_extract.py::test_symbolinfo_holds_required_fields -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'climt._docs'`.

- [ ] **Step 3: Write minimal implementation**

`climt/_docs/__init__.py`:
```python
```

`climt/_docs/extract.py`:
```python
"""AST-based metadata extractor for climt symbols. Imports nothing from climt."""
from __future__ import annotations

import ast
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class SymbolInfo:
    name: str
    kind: str  # "class" or "function"
    source_file: str
    source_line: int
    signature: str
    docstring: Optional[str] = None
    full_docstring: Optional[str] = None
    bases: list[str] = field(default_factory=list)
    sympl_properties: dict = field(default_factory=dict)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/_docs/test_extract.py::test_symbolinfo_holds_required_fields -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add climt/_docs/__init__.py climt/_docs/extract.py tests/_docs/__init__.py tests/_docs/test_extract.py
git commit -m "feat(docs): add SymbolInfo dataclass for climt.docs metadata"
```

---

### Task 2: Extract class signature, docstring, bases from a file

**Files:**
- Modify: `climt/_docs/extract.py`
- Modify: `tests/_docs/test_extract.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/_docs/test_extract.py`:
```python
from pathlib import Path
from climt._docs.extract import extract_file


def test_extract_class_with_init_signature_docstring_bases(tmp_path):
    src = tmp_path / "demo.py"
    src.write_text(
        '''
class Foo(Base):
    """Short summary.

    Detailed paragraph that goes on.
    """
    def __init__(self, a=1, b="x", **kwargs):
        pass
'''
    )
    [sym] = extract_file(src, rel_path="demo.py")
    assert sym.name == "Foo"
    assert sym.kind == "class"
    assert sym.signature == "Foo(a=1, b='x', **kwargs)"
    assert sym.docstring == "Short summary."
    assert "Detailed paragraph" in sym.full_docstring
    assert sym.bases == ["Base"]
    assert sym.source_line == 2
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/_docs/test_extract.py::test_extract_class_with_init_signature_docstring_bases -v`
Expected: FAIL with `ImportError: cannot import name 'extract_file'`.

- [ ] **Step 3: Write minimal implementation**

Append to `climt/_docs/extract.py`:
```python
from pathlib import Path


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


def _signature_of(fn: ast.FunctionDef | ast.AsyncFunctionDef) -> str:
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


def _first_doc_line(node) -> Optional[str]:
    doc = ast.get_docstring(node)
    if not doc:
        return None
    for line in doc.splitlines():
        line = line.strip()
        if line:
            return line[:280]
    return None


def _full_doc(node, max_chars: int = 4000) -> Optional[str]:
    doc = ast.get_docstring(node)
    if not doc:
        return None
    doc = doc.strip()
    return doc[:max_chars] + ("..." if len(doc) > max_chars else "")


def _base_names(cls_node: ast.ClassDef) -> list[str]:
    out = []
    for b in cls_node.bases:
        try:
            out.append(ast.unparse(b))
        except Exception:
            continue
    return out


def extract_file(path: Path, rel_path: str) -> list[SymbolInfo]:
    """Parse one Python file and return SymbolInfo entries for every top-level
    class and function. Returns empty list if the file fails to parse."""
    try:
        tree = ast.parse(path.read_text())
    except (SyntaxError, OSError, UnicodeDecodeError):
        return []

    out: list[SymbolInfo] = []
    for node in tree.body:
        if isinstance(node, ast.ClassDef):
            init = next(
                (b for b in node.body
                 if isinstance(b, (ast.FunctionDef, ast.AsyncFunctionDef))
                 and b.name == "__init__"),
                None,
            )
            sig = (
                f"{node.name}{_signature_of(init)[len('__init__'):]}"
                if init else f"{node.name}()"
            )
            out.append(SymbolInfo(
                name=node.name,
                kind="class",
                source_file=rel_path,
                source_line=node.lineno,
                signature=sig,
                docstring=_first_doc_line(node)
                          or (_first_doc_line(init) if init else None),
                full_docstring=_full_doc(node)
                               or (_full_doc(init) if init else None),
                bases=_base_names(node),
                sympl_properties={},
            ))
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            out.append(SymbolInfo(
                name=node.name,
                kind="function",
                source_file=rel_path,
                source_line=node.lineno,
                signature=_signature_of(node),
                docstring=_first_doc_line(node),
                full_docstring=_full_doc(node),
            ))
    return out
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/_docs/test_extract.py -v`
Expected: both tests PASS.

- [ ] **Step 5: Commit**

```bash
git add climt/_docs/extract.py tests/_docs/test_extract.py
git commit -m "feat(docs): extract class signature, docstring, bases from source"
```

---

### Task 3: Extract sympl property dicts (input/output/diagnostic/tendency)

**Files:**
- Modify: `climt/_docs/extract.py`
- Modify: `tests/_docs/test_extract.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/_docs/test_extract.py`:
```python
def test_extract_sympl_properties(tmp_path):
    src = tmp_path / "comp.py"
    src.write_text(
        '''
class C:
    input_properties = {
        "air_temperature": {"units": "degK", "dims": ["*", "mid_levels"]},
    }
    diagnostic_properties = {
        "olr": {"units": "W m^-2", "dims": ["*"]},
    }
    def __init__(self):
        pass
'''
    )
    [sym] = extract_file(src, rel_path="comp.py")
    assert sym.sympl_properties["input_properties"]["air_temperature"]["units"] == "degK"
    assert sym.sympl_properties["input_properties"]["air_temperature"]["dims"] == ["*", "mid_levels"]
    assert sym.sympl_properties["diagnostic_properties"]["olr"]["units"] == "W m^-2"
    assert "output_properties" not in sym.sympl_properties
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/_docs/test_extract.py::test_extract_sympl_properties -v`
Expected: FAIL — `sympl_properties` is currently always `{}`.

- [ ] **Step 3: Write minimal implementation**

Add to `climt/_docs/extract.py` (above `extract_file`):
```python
SYMPL_PROPERTY_ATTRS = (
    "input_properties",
    "output_properties",
    "diagnostic_properties",
    "tendency_properties",
)


def _extract_sympl_properties(cls_node: ast.ClassDef) -> dict:
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
            inner: dict = {}
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
```

In `extract_file`, modify the class branch to populate `sympl_properties`:
```python
            out.append(SymbolInfo(
                name=node.name,
                kind="class",
                source_file=rel_path,
                source_line=node.lineno,
                signature=sig,
                docstring=_first_doc_line(node)
                          or (_first_doc_line(init) if init else None),
                full_docstring=_full_doc(node)
                               or (_full_doc(init) if init else None),
                bases=_base_names(node),
                sympl_properties=_extract_sympl_properties(node),
            ))
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/_docs/test_extract.py -v`
Expected: all three tests PASS.

- [ ] **Step 5: Commit**

```bash
git add climt/_docs/extract.py tests/_docs/test_extract.py
git commit -m "feat(docs): extract sympl input/output/diagnostic/tendency property dicts"
```

---

### Task 4: Extract per-file imports table

**Files:**
- Modify: `climt/_docs/extract.py`
- Modify: `tests/_docs/test_extract.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/_docs/test_extract.py`:
```python
from climt._docs.extract import extract_imports


def test_extract_imports_records_aliases(tmp_path):
    src = tmp_path / "m.py"
    src.write_text(
        '''
from sympl import set_backend, AdamsBashforth as AB
from . import component as comp
import numpy as np
'''
    )
    imports = extract_imports(src)
    assert imports["set_backend"] == "from sympl import set_backend"
    assert imports["AB"] == "from sympl import AdamsBashforth"
    assert imports["comp"] == "from . import component"
    assert imports["np"] == "import numpy"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/_docs/test_extract.py::test_extract_imports_records_aliases -v`
Expected: FAIL — `extract_imports` doesn't exist.

- [ ] **Step 3: Write minimal implementation**

Append to `climt/_docs/extract.py`:
```python
def extract_imports(path: Path) -> dict[str, str]:
    """Return {alias_name: 'from MOD import NAME' or 'import NAME'} for one file."""
    try:
        tree = ast.parse(path.read_text())
    except (SyntaxError, OSError, UnicodeDecodeError):
        return {}
    out: dict[str, str] = {}
    for node in tree.body:
        if isinstance(node, ast.ImportFrom):
            mod_prefix = "." * node.level + (node.module or "")
            for n in node.names:
                key = n.asname or n.name
                out[key] = f"from {mod_prefix} import {n.name}"
        elif isinstance(node, ast.Import):
            for n in node.names:
                key = n.asname or n.name.split(".")[0]
                out[key] = f"import {n.name}"
    return out
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/_docs/test_extract.py -v`
Expected: all four tests PASS.

- [ ] **Step 5: Commit**

```bash
git add climt/_docs/extract.py tests/_docs/test_extract.py
git commit -m "feat(docs): extract per-file imports table with alias resolution"
```

---

### Task 5: Build registry — walk climt tree, index by symbol name

**Files:**
- Create: `climt/_docs/registry.py`
- Create: `tests/_docs/test_registry.py`

- [ ] **Step 1: Write the failing test**

`tests/_docs/test_registry.py`:
```python
from climt._docs.registry import build_registry


def test_build_registry_finds_known_climt_classes():
    reg = build_registry()
    names = {s.name for s in reg.symbols}
    # These three are stable public components.
    assert "RRTMGLongwave" in names
    assert "EmanuelConvection" in names
    assert "get_grid" in names


def test_registry_resolve_direct_hit():
    reg = build_registry()
    hits = reg.resolve("RRTMGLongwave")
    assert len(hits) >= 1
    assert any(h.kind == "class" and h.source_file.endswith("rrtmg/lw/component.py")
               for h in hits)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/_docs/test_registry.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'climt._docs.registry'`.

- [ ] **Step 3: Write minimal implementation**

`climt/_docs/registry.py`:
```python
"""Walk the climt source tree and index all extractable symbols."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from climt._docs.extract import SymbolInfo, extract_file, extract_imports

CLIMT_ROOT = Path(__file__).resolve().parent.parent  # .../climt/


@dataclass
class Registry:
    symbols: list[SymbolInfo] = field(default_factory=list)
    imports_by_file: dict[str, dict[str, str]] = field(default_factory=dict)

    def resolve(self, query: str) -> list[SymbolInfo]:
        target = query.rstrip("()").lower()
        return [s for s in self.symbols if s.name.lower() == target]


def build_registry() -> Registry:
    reg = Registry()
    for py in CLIMT_ROOT.rglob("*.py"):
        if any(part.startswith(".") or part == "__pycache__" for part in py.parts):
            continue
        rel = str(py.relative_to(CLIMT_ROOT.parent))  # e.g. "climt/_core/x.py"
        reg.symbols.extend(extract_file(py, rel_path=rel))
        imports = extract_imports(py)
        if imports:
            reg.imports_by_file[rel] = imports
    return reg
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/_docs/test_registry.py -v`
Expected: both tests PASS.

- [ ] **Step 5: Commit**

```bash
git add climt/_docs/registry.py tests/_docs/test_registry.py
git commit -m "feat(docs): build symbol registry by walking climt tree"
```

---

### Task 6: Resolve aliased symbols through `_components/__init__.py` re-exports

Context: graphify's earlier failure mode was that `EmanuelConvection` is re-exported from `climt/_components/__init__.py` via `from .emanuel import EmanuelConvection`, but the *class* `EmanuelConvection` lives in `climt/_components/emanuel/component.py`. A registry that only indexes by direct file location misses the chain. This task adds chain-following.

**Files:**
- Modify: `climt/_docs/registry.py`
- Modify: `tests/_docs/test_registry.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/_docs/test_registry.py`:
```python
def test_registry_resolve_via_components_reexport():
    reg = build_registry()
    # SlabSurface is re-exported from climt/_components/__init__.py
    # and from climt/__init__.py — both paths should find the class def.
    hits = reg.resolve("SlabSurface")
    sources = {h.source_file for h in hits if h.kind == "class"}
    assert "climt/_components/slab_surface.py" in sources
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/_docs/test_registry.py::test_registry_resolve_via_components_reexport -v`
Expected: PASS if direct lookup already works (since `SlabSurface` class is a top-level class in `slab_surface.py` and the walker picks it up). If it FAILS for a re-exported symbol that's actually defined in a sub-package (`EmanuelConvection` → `emanuel/component.py`), the registry needs the chain-follow logic. Run a quick `python -c "from climt._docs.registry import build_registry; r=build_registry(); print([s.source_file for s in r.resolve('EmanuelConvection')])"` to confirm. If the direct walk finds it, this task is a no-op and the test still validates correctness.

- [ ] **Step 3: If a resolve gap exists, add re-export following**

Only execute this step if Step 2 revealed a missing symbol. Append to `Registry.resolve`:
```python
    def resolve(self, query: str) -> list[SymbolInfo]:
        target = query.rstrip("()").lower()
        hits = [s for s in self.symbols if s.name.lower() == target]
        if hits:
            return hits
        # Follow re-exports: find any file that imports `query` and look in
        # the source module the import names.
        for file, imports in self.imports_by_file.items():
            if query in imports:
                stmt = imports[query]
                # "from <mod> import <name>" -> resolve relative to `file`'s dir
                parts = stmt.split()
                if parts[0] != "from":
                    continue
                mod = parts[1].lstrip(".")
                base = Path(file).parent
                for candidate in (
                    base / f"{mod}.py",
                    base / mod / "__init__.py",
                    base / mod / "component.py",
                ):
                    if (CLIMT_ROOT.parent / candidate).exists():
                        rel = str(candidate)
                        sub_hits = [s for s in self.symbols
                                    if s.source_file == rel and s.name == query]
                        if sub_hits:
                            return sub_hits
        return []
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/_docs/test_registry.py -v`
Expected: all registry tests PASS.

- [ ] **Step 5: Commit**

```bash
git add climt/_docs/registry.py tests/_docs/test_registry.py
git commit -m "feat(docs): follow component re-export chain when resolving names"
```

---

### Task 7: Imports-only fallback for external symbols (e.g. `sympl.set_backend`)

**Files:**
- Modify: `climt/_docs/registry.py`
- Modify: `tests/_docs/test_registry.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/_docs/test_registry.py`:
```python
def test_registry_resolve_external_symbol_via_imports():
    reg = build_registry()
    info = reg.resolve_external("set_backend")
    assert info is not None
    assert info["import_statement"].startswith("from sympl")
    assert info["files"]  # at least one climt file imports it
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/_docs/test_registry.py::test_registry_resolve_external_symbol_via_imports -v`
Expected: FAIL — `resolve_external` doesn't exist.

- [ ] **Step 3: Write minimal implementation**

Append to `Registry`:
```python
    def resolve_external(self, query: str) -> dict | None:
        """For symbols not defined in climt but imported from elsewhere
        (e.g. sympl.set_backend), return where they come from."""
        statements: dict[str, list[str]] = {}
        for file, imports in self.imports_by_file.items():
            if query in imports:
                statements.setdefault(imports[query], []).append(file)
        if not statements:
            return None
        # Pick the most-frequent import statement; surface all files
        top_stmt = max(statements, key=lambda s: len(statements[s]))
        return {
            "import_statement": top_stmt,
            "files": statements[top_stmt],
            "all_statements": statements,
        }
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/_docs/test_registry.py -v`
Expected: all tests PASS.

- [ ] **Step 5: Commit**

```bash
git add climt/_docs/registry.py tests/_docs/test_registry.py
git commit -m "feat(docs): fall back to imports table for external symbols"
```

---

### Task 8: Plain-text formatter

**Files:**
- Create: `climt/_docs/format.py`
- Create: `tests/_docs/test_format.py`

- [ ] **Step 1: Write the failing test**

`tests/_docs/test_format.py`:
```python
from climt._docs.extract import SymbolInfo
from climt._docs.format import format_symbol, format_external


def test_format_symbol_includes_signature_bases_properties():
    s = SymbolInfo(
        name="Foo",
        kind="class",
        source_file="climt/x.py",
        source_line=10,
        signature="Foo(a=1)",
        docstring="Short.",
        full_docstring="Short.\n\nDetail.",
        bases=["Bar"],
        sympl_properties={
            "input_properties": {
                "air_temperature": {"units": "degK", "dims": ["*", "mid_levels"]},
            },
        },
    )
    out = format_symbol(s, importers={"from climt import Foo": ["tests/x.py"]})
    assert "Foo  [climt/x.py L10]" in out
    assert "signature: Foo(a=1)" in out
    assert "bases:     Bar" in out
    assert "air_temperature" in out
    assert "units=degK" in out
    assert "from climt import Foo" in out


def test_format_external():
    out = format_external("set_backend", {
        "import_statement": "from sympl import set_backend",
        "files": ["tests/test_components.py", "examples/x.py"],
    })
    assert "set_backend" in out
    assert "from sympl import set_backend" in out
    assert "2 file" in out
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/_docs/test_format.py -v`
Expected: FAIL — `format.py` doesn't exist.

- [ ] **Step 3: Write minimal implementation**

`climt/_docs/format.py`:
```python
"""Plain-text rendering for SymbolInfo."""
from __future__ import annotations

from climt._docs.extract import SymbolInfo


def format_symbol(s: SymbolInfo, importers: dict[str, list[str]] | None = None) -> str:
    lines: list[str] = []
    header = f"\n{s.name}  [{s.source_file}"
    if s.source_line:
        header += f" L{s.source_line}"
    header += "]"
    lines.append(header)
    if s.signature:
        lines.append(f"  signature: {s.signature}")
    if s.docstring:
        lines.append(f"  doc:       {s.docstring}")
    if s.kind:
        lines.append(f"  kind:      {s.kind}")
    if s.bases:
        lines.append(f"  bases:     {', '.join(s.bases)}")
    if s.full_docstring and s.full_docstring != s.docstring:
        lines.append("  docstring:")
        for line in s.full_docstring.splitlines():
            lines.append(f"      {line}")
    for attr in ("input_properties", "output_properties",
                 "diagnostic_properties", "tendency_properties"):
        fields = s.sympl_properties.get(attr) if s.sympl_properties else None
        if not fields:
            continue
        lines.append(f"  {attr}:")
        for fname, meta in fields.items():
            units = meta.get("units", "?")
            dims = meta.get("dims", "?")
            lines.append(f"      {fname:50s} units={units!s:18s} dims={dims}")
    if importers:
        lines.append("")
        lines.append("  imported as:")
        for stmt, files in importers.items():
            n = len(files)
            lines.append(f"    {stmt}   ({n} file{'s' if n != 1 else ''}, e.g. {files[0]})")
    return "\n".join(lines)


def format_external(name: str, info: dict) -> str:
    files = info.get("files", [])
    n = len(files)
    out = [f"\n{name}  (no climt definition — resolved via imports)"]
    out.append(f"  {info['import_statement']}   "
               f"({n} file{'s' if n != 1 else ''}, e.g. {files[0] if files else '?'})")
    return "\n".join(out)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/_docs/test_format.py -v`
Expected: both tests PASS.

- [ ] **Step 5: Commit**

```bash
git add climt/_docs/format.py tests/_docs/test_format.py
git commit -m "feat(docs): plain-text formatter for SymbolInfo and external imports"
```

---

### Task 9: Public API in `climt/_docs/__init__.py`

**Files:**
- Modify: `climt/_docs/__init__.py`
- Create: `tests/_docs/test_public_api.py`

- [ ] **Step 1: Write the failing test**

`tests/_docs/test_public_api.py`:
```python
from climt._docs import show, list_symbols


def test_show_real_climt_class():
    out = show("RRTMGLongwave")
    assert "RRTMGLongwave" in out
    assert "signature:" in out


def test_show_external_symbol():
    out = show("set_backend")
    assert "from sympl import set_backend" in out


def test_show_unknown_symbol_returns_message():
    out = show("ThisDoesNotExist")
    assert "not found" in out.lower() or "no climt" in out.lower()


def test_list_symbols_returns_known_names():
    names = list_symbols()
    assert "RRTMGLongwave" in names
    assert "EmanuelConvection" in names
    assert "get_grid" in names
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/_docs/test_public_api.py -v`
Expected: FAIL — `show` and `list_symbols` don't exist.

- [ ] **Step 3: Write minimal implementation**

`climt/_docs/__init__.py`:
```python
"""Introspect climt's component constructors, docstrings, and sympl property dicts.

Usage:
    >>> from climt._docs import show
    >>> print(show("RRTMGLongwave"))

CLI:
    $ python -m climt.docs show RRTMGLongwave
"""
from __future__ import annotations

from climt._docs.registry import build_registry
from climt._docs.format import format_symbol, format_external

_REGISTRY = None


def _get_registry():
    global _REGISTRY
    if _REGISTRY is None:
        _REGISTRY = build_registry()
    return _REGISTRY


def show(symbol: str) -> str:
    """Return plain-text metadata for `symbol` — signature, docstring, bases,
    sympl input/output/diagnostic/tendency properties, and where it's imported
    from. Falls back to import-trail resolution for symbols defined outside
    climt (e.g. ``sympl.set_backend``).
    """
    reg = _get_registry()
    hits = reg.resolve(symbol)
    if hits:
        # Collect importers across the whole tree
        importers: dict[str, list[str]] = {}
        for file, imports in reg.imports_by_file.items():
            if symbol in imports:
                importers.setdefault(imports[symbol], []).append(file)
        return "\n".join(format_symbol(h, importers=importers) for h in hits)
    ext = reg.resolve_external(symbol)
    if ext:
        return format_external(symbol, ext)
    return f"No climt symbol found for {symbol!r}."


def list_symbols() -> list[str]:
    """Return the sorted list of all symbol names known to the registry."""
    reg = _get_registry()
    return sorted({s.name for s in reg.symbols})
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/_docs/test_public_api.py -v`
Expected: all four tests PASS.

- [ ] **Step 5: Commit**

```bash
git add climt/_docs/__init__.py tests/_docs/test_public_api.py
git commit -m "feat(docs): public show() and list_symbols() API"
```

---

### Task 10: `python -m climt.docs` CLI

**Files:**
- Create: `climt/_docs/__main__.py`
- Create: `tests/_docs/test_cli.py`

- [ ] **Step 1: Write the failing test**

`tests/_docs/test_cli.py`:
```python
import subprocess
import sys


def _run(*args) -> tuple[int, str]:
    res = subprocess.run(
        [sys.executable, "-m", "climt.docs", *args],
        capture_output=True, text=True,
    )
    return res.returncode, res.stdout + res.stderr


def test_cli_show_known_symbol():
    code, out = _run("show", "RRTMGLongwave")
    assert code == 0
    assert "RRTMGLongwave" in out
    assert "signature:" in out


def test_cli_show_unknown_symbol_nonzero_exit():
    code, out = _run("show", "DefinitelyNotASymbol")
    assert code != 0
    assert "not found" in out.lower() or "no climt" in out.lower()


def test_cli_list_includes_known_names():
    code, out = _run("list")
    assert code == 0
    assert "RRTMGLongwave" in out
    assert "EmanuelConvection" in out


def test_cli_list_with_filter():
    code, out = _run("list", "--filter", "RRTMG")
    assert code == 0
    assert "RRTMGLongwave" in out
    assert "EmanuelConvection" not in out


def test_cli_no_args_prints_usage():
    code, out = _run()
    assert code != 0
    assert "usage" in out.lower() or "show" in out.lower()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/_docs/test_cli.py -v`
Expected: FAIL — `__main__.py` doesn't exist.

- [ ] **Step 3: Write minimal implementation**

`climt/_docs/__main__.py`:
```python
"""CLI entry: `python -m climt.docs <command> [args]`."""
from __future__ import annotations

import argparse
import sys

from climt._docs import show, list_symbols


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="python -m climt.docs",
        description="Introspect climt component constructors and state-field contracts.",
    )
    sub = p.add_subparsers(dest="command", required=True)

    p_show = sub.add_parser("show", help="Show metadata for one symbol.")
    p_show.add_argument("symbol", help="Symbol name (e.g. 'RRTMGLongwave', 'set_backend').")

    p_list = sub.add_parser("list", help="List all known symbol names.")
    p_list.add_argument(
        "--filter", default=None,
        help="Case-insensitive substring filter.",
    )

    args = p.parse_args(argv)
    if args.command == "show":
        out = show(args.symbol)
        print(out)
        return 0 if "no climt symbol found" not in out.lower() else 1
    if args.command == "list":
        names = list_symbols()
        if args.filter:
            f = args.filter.lower()
            names = [n for n in names if f in n.lower()]
        print("\n".join(names))
        return 0
    return 2  # unreachable thanks to required=True


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/_docs/test_cli.py -v`
Expected: all five tests PASS.

- [ ] **Step 5: Commit**

```bash
git add climt/_docs/__main__.py tests/_docs/test_cli.py
git commit -m "feat(docs): `python -m climt.docs show|list` CLI"
```

---

### Task 11: On-disk cache keyed by source mtime hash

Rationale: rebuilding the registry walks ~80 .py files each call. That's ~30 ms cold but adds up if an LLM tool calls `show` repeatedly. Cache the registry to JSON next to `_docs/`, invalidated by a hash of all source mtimes.

**Files:**
- Modify: `climt/_docs/registry.py`
- Create: `tests/_docs/test_cache.py`

- [ ] **Step 1: Write the failing test**

`tests/_docs/test_cache.py`:
```python
import time
from pathlib import Path

from climt._docs.registry import build_registry, _cache_path


def test_second_build_uses_cache(tmp_path, monkeypatch):
    # Point the cache at a tmp location so we don't clobber the real one.
    cache = tmp_path / "_cache.json"
    monkeypatch.setattr("climt._docs.registry._cache_path", lambda: cache)
    reg1 = build_registry()
    assert cache.exists()
    mtime1 = cache.stat().st_mtime
    time.sleep(0.01)
    reg2 = build_registry()
    assert cache.stat().st_mtime == mtime1  # cache reused, not rewritten
    assert len(reg2.symbols) == len(reg1.symbols)


def test_cache_invalidated_when_source_changes(tmp_path, monkeypatch):
    cache = tmp_path / "_cache.json"
    monkeypatch.setattr("climt._docs.registry._cache_path", lambda: cache)
    build_registry()
    first = cache.stat().st_mtime
    # Touch any climt .py file to bump its mtime.
    import climt
    sample = Path(climt.__file__).parent / "__init__.py"
    sample.touch()
    time.sleep(0.01)
    build_registry()
    assert cache.stat().st_mtime > first
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/_docs/test_cache.py -v`
Expected: FAIL — `_cache_path` and cache logic don't exist.

- [ ] **Step 3: Write minimal implementation**

Modify `climt/_docs/registry.py`:
```python
import json
from dataclasses import asdict


def _cache_path() -> Path:
    return Path(__file__).resolve().parent / "_cache.json"


def _source_fingerprint() -> str:
    """Cheap hash of every climt .py mtime. Order-stable."""
    parts: list[str] = []
    for py in sorted(CLIMT_ROOT.rglob("*.py")):
        if any(p.startswith(".") or p == "__pycache__" for p in py.parts):
            continue
        parts.append(f"{py}:{py.stat().st_mtime_ns}")
    import hashlib
    return hashlib.sha1("\n".join(parts).encode()).hexdigest()
```

Then refactor `build_registry`:
```python
def build_registry() -> Registry:
    cache = _cache_path()
    fp = _source_fingerprint()
    if cache.exists():
        try:
            data = json.loads(cache.read_text())
            if data.get("fingerprint") == fp:
                return Registry(
                    symbols=[SymbolInfo(**s) for s in data["symbols"]],
                    imports_by_file=data["imports_by_file"],
                )
        except (json.JSONDecodeError, KeyError, TypeError):
            pass  # fall through and rebuild
    reg = Registry()
    for py in CLIMT_ROOT.rglob("*.py"):
        if any(part.startswith(".") or part == "__pycache__" for part in py.parts):
            continue
        rel = str(py.relative_to(CLIMT_ROOT.parent))
        reg.symbols.extend(extract_file(py, rel_path=rel))
        imports = extract_imports(py)
        if imports:
            reg.imports_by_file[rel] = imports
    cache.write_text(json.dumps({
        "fingerprint": fp,
        "symbols": [asdict(s) for s in reg.symbols],
        "imports_by_file": reg.imports_by_file,
    }))
    return reg
```

Also add `climt/_docs/_cache.json` to `.gitignore`:
```bash
echo "climt/_docs/_cache.json" >> .gitignore
```

- [ ] **Step 4: Run all tests to verify everything still passes**

Run: `pytest tests/_docs/ -v`
Expected: every test passes (including the new cache tests). The previous `test_registry.py` tests still work because the cache is transparent.

- [ ] **Step 5: Commit**

```bash
git add climt/_docs/registry.py tests/_docs/test_cache.py .gitignore
git commit -m "feat(docs): cache registry on disk keyed by source mtime hash"
```

---

### Task 12: Document the feature

**Files:**
- Modify: `README.md`
- Modify: `CLAUDE.md`
- Modify: `GEMINI.md`

- [ ] **Step 1: Add a README section**

Append the following section to `README.md` (find an appropriate placement near other usage sections; if there's a "Documentation" or "Development" section, put it there; otherwise append before any final "License" section):

```markdown
## Introspecting components from the command line

`climt` ships with a built-in metadata browser. From any environment that has
climt installed:

```bash
python -m climt.docs show RRTMGLongwave
python -m climt.docs show set_backend
python -m climt.docs list --filter Cork
```

`show` prints the constructor signature, full docstring, base class, the
sympl `input_properties` / `output_properties` / `diagnostic_properties` /
`tendency_properties` dicts (with units and dimensions for every state
field), and every import statement used to bring the symbol into other
files. It works for any climt class or function and for symbols imported
from sympl/numpy/etc.

This is especially useful when pairing with an LLM coding assistant — it
gives the assistant a reliable source of truth for kwargs and state-field
names without having to grep the source.

Programmatic use: `from climt._docs import show, list_symbols`.
```

- [ ] **Step 2: Update CLAUDE.md and GEMINI.md to point at the new CLI**

In both `CLAUDE.md` and `GEMINI.md`, replace the existing line about `python scripts/augment_graph.py show "<Symbol>"` with:

```markdown
- For constructor signatures, kwargs, docstrings, state-field contracts (units + dims), and import paths of any climt symbol — or any symbol imported by climt code (e.g. `sympl.set_backend`) — run `python -m climt.docs show "<Symbol>"`. This ships inside climt itself; no extra tooling needed. The augmenter in `scripts/augment_graph.py` was the prototype; `climt.docs` is its successor and now the canonical path.
```

(Leave the rest of the graphify-related rules in place — `climt.docs` complements graphify rather than replacing it.)

- [ ] **Step 3: Verify the CLI output matches what README claims**

Run: `python -m climt.docs show RRTMGLongwave`
Expected: output contains `signature:`, `bases:`, `input_properties:`, `imported as:` sections. If anything is missing, the earlier tasks have a bug — go fix.

Run: `python -m climt.docs show set_backend`
Expected: `(no climt definition — resolved via imports)` and `from sympl import set_backend`.

- [ ] **Step 4: Commit**

```bash
git add README.md CLAUDE.md GEMINI.md
git commit -m "docs: announce climt.docs CLI in README and AGENT instructions"
```

---

## Self-Review checklist (executed during plan write)

- **Spec coverage:** All four answered design points covered — plain text only (Task 8), no Sphinx integration, AST-only with no climt import at runtime (Task 2-4), sympl_properties extraction (Task 3).
- **Placeholder scan:** No TBDs, no "handle edge cases" lines, every code block is complete.
- **Type consistency:** `SymbolInfo` field names are stable from Task 1 through Task 11. `Registry.resolve`, `Registry.resolve_external`, `build_registry`, `show`, `list_symbols`, `format_symbol`, `format_external` all keep their signatures across tasks.
- **Bite-size:** Every task is 5 steps, each step is one action. No multi-action steps.
