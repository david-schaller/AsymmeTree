"""Generate an interactive dependency diagram for ``example_simulations.py``.

The script performs a lightweight static analysis of the Python sources under the local
``development notes/Example code`` tree and writes a self-contained HTML file. The diagram is
organized by module lane so that the relationships between ``example_simulations.py``,
``asymmetree.hologenome``, and ``asymmetree.holoevolve`` remain readable.
"""

from __future__ import annotations

import argparse
import ast
import json
from dataclasses import dataclass, field
from pathlib import Path


DEFAULT_ENTRY = "example_simulations.py"
DEFAULT_OUTPUT = "dependency_diagram.html"

PREFERRED_LANES = (
    "example_simulations",
    "asymmetree.hologenome",
    "asymmetree.holoevolve",
    "asymmetree.treeevolve",
    "asymmetree.utils",
    "asymmetree.other",
    "external",
)

LANE_STYLES = {
    "example_simulations": {"color": "#d4a72c", "band": "#fff6d5"},
    "asymmetree.hologenome": {"color": "#1f78b4", "band": "#e8f3fb"},
    "asymmetree.holoevolve": {"color": "#2f9e44", "band": "#eaf7ed"},
    "asymmetree.treeevolve": {"color": "#d97706", "band": "#fff2e2"},
    "asymmetree.utils": {"color": "#0f766e", "band": "#e6f6f4"},
    "asymmetree.other": {"color": "#6b7280", "band": "#f3f4f6"},
    "external": {"color": "#7c7c7c", "band": "#f7f7f7"},
}

KIND_RANK = {
    "module": 0,
    "class": 1,
    "main": 2,
    "function": 3,
    "method": 4,
    "external-class": 5,
    "external-function": 6,
}


@dataclass
class ImportBinding:
    """Representation of one import alias."""

    alias: str
    module: str
    symbol: str | None
    lineno: int


@dataclass
class DefinitionInfo:
    """Representation of a top-level function, class, method, or ``__main__`` block."""

    node_id: str
    module: str
    name: str
    kind: str
    path: Path
    lineno: int
    end_lineno: int
    doc: str
    class_name: str | None = None
    ast_body: list[ast.stmt] = field(default_factory=list)
    private: bool = False

    @property
    def label(self) -> str:
        """Short label used in the HTML diagram."""
        if self.kind == "module":
            return self.module
        if self.kind == "class":
            return self.name
        if self.kind == "external-class":
            return self.name
        if self.kind == "method":
            return f"{self.class_name}.{self.name}()"
        if self.kind == "main":
            return "__main__"
        return f"{self.name}()"

    @property
    def full_name(self) -> str:
        """Qualified name of the definition."""
        if self.kind == "module":
            return self.module
        if self.class_name:
            return f"{self.module}.{self.class_name}.{self.name}"
        return f"{self.module}.{self.name}"


@dataclass
class ModuleInfo:
    """Static information extracted from one Python module."""

    name: str
    path: Path
    doc: str
    is_package: bool
    imports: dict[str, ImportBinding] = field(default_factory=dict)
    top_level_symbols: dict[str, str] = field(default_factory=dict)
    class_symbols: dict[str, str] = field(default_factory=dict)
    method_symbols: dict[tuple[str, str], str] = field(default_factory=dict)
    definitions: dict[str, DefinitionInfo] = field(default_factory=dict)


@dataclass
class EdgeInfo:
    """Graph edge between nodes."""

    source: str
    target: str
    edge_type: str
    calls: set[str] = field(default_factory=set)
    lines: list[int] = field(default_factory=list)
    count: int = 0

    def register(self, raw_call: str, lineno: int | None) -> None:
        """Record one additional call or relation occurrence."""
        if raw_call:
            self.calls.add(raw_call)
        if lineno is not None:
            self.lines.append(lineno)
        self.count += 1


class RepositoryModel:
    """Model of the local source tree plus the derived dependency graph."""

    def __init__(self, root_dir: Path) -> None:
        self.root_dir = root_dir.resolve()
        self.module_paths: dict[str, Path] = {}
        self.modules: dict[str, ModuleInfo] = {}
        self.nodes: dict[str, DefinitionInfo] = {}
        self.module_nodes: dict[str, DefinitionInfo] = {}
        self.edges: dict[tuple[str, str, str], EdgeInfo] = {}
        self._scan_module_paths()
        self._parse_modules()
        self._analyze_calls()
        self._add_module_relation_edges()

    def build_diagram(self, entry_module: str) -> dict[str, object]:
        """Return the graph data needed by the HTML viewer."""
        roots = [
            f"{entry_module}.run_example_simulations",
            f"{entry_module}.__main__",
        ]
        reachable = self._reachable_nodes(roots)
        included = self._expand_with_structure(reachable)
        graph_nodes = self._serialize_nodes(included, roots)
        graph_edges = self._serialize_edges(included)
        lanes = self._layout_nodes(graph_nodes)

        return {
            "entry_module": entry_module,
            "nodes": graph_nodes,
            "edges": graph_edges,
            "lanes": lanes,
            "defaultSelection": f"{entry_module}.run_example_simulations",
        }

    def _scan_module_paths(self) -> None:
        """Index all Python files below the local example-code root."""
        for path in sorted(self.root_dir.rglob("*.py")):
            if "__pycache__" in path.parts:
                continue
            self.module_paths[self._module_name(path)] = path

    def _parse_modules(self) -> None:
        """Parse all indexed modules and collect top-level definitions."""
        for module_name, path in self.module_paths.items():
            source = path.read_text(encoding="utf-8")
            tree = ast.parse(source, filename=str(path))
            module = ModuleInfo(
                name=module_name,
                path=path,
                doc=ast.get_docstring(tree) or "",
                is_package=path.name == "__init__.py",
            )
            self.modules[module_name] = module
            self.module_nodes[module_name] = DefinitionInfo(
                node_id=self._module_node_id(module_name),
                module=module_name,
                name=module_name,
                kind="module",
                path=path,
                lineno=1,
                end_lineno=len(source.splitlines()) or 1,
                doc=module.doc,
            )
            collector = ModuleCollector(self, module)
            collector.visit(tree)

    def _analyze_calls(self) -> None:
        """Analyze bodies of functions and methods to collect call edges."""
        for definition in list(self.nodes.values()):
            if definition.kind in {"function", "method", "main"}:
                analyzer = BodyCallAnalyzer(self, definition)
                for statement in definition.ast_body:
                    analyzer.visit(statement)

    def _add_module_relation_edges(self) -> None:
        """Add containment, import, and re-export edges between modules and definitions."""
        for module in self.modules.values():
            module_id = self._module_node_id(module.name)

            for definition in module.definitions.values():
                if definition.kind == "method":
                    class_id = module.class_symbols[definition.class_name or ""]
                    self._register_edge(
                        class_id,
                        definition.node_id,
                        "contains",
                        raw_call="contains",
                        lineno=definition.lineno,
                    )
                else:
                    self._register_edge(
                        module_id,
                        definition.node_id,
                        "contains",
                        raw_call="contains",
                        lineno=definition.lineno,
                    )

            for binding in module.imports.values():
                if binding.symbol is None:
                    if binding.module in self.modules:
                        self._register_edge(
                            module_id,
                            self._module_node_id(binding.module),
                            "imports",
                            raw_call=binding.alias,
                            lineno=binding.lineno,
                        )
                    continue

                target = self.resolve_module_member(binding.module, binding.symbol)
                edge_type = "reexports" if module.is_package else "imports"
                self._register_edge(
                    module_id,
                    target,
                    edge_type,
                    raw_call=binding.alias,
                    lineno=binding.lineno,
                )

    def resolve_name(self, module_name: str, name: str) -> str | None:
        """Resolve a bare name inside one module."""
        module = self.modules[module_name]

        if name in module.top_level_symbols:
            return module.top_level_symbols[name]

        if name in module.imports:
            return self._resolve_import_binding(module.imports[name])

        return None

    def resolve_module_member(self, module_name: str, member: str) -> str:
        """Resolve one named symbol exported by a module."""
        if module_name in self.modules:
            module = self.modules[module_name]
            if member in module.top_level_symbols:
                return module.top_level_symbols[member]
            if member in module.imports:
                return self._resolve_import_binding(module.imports[member])

        return self._ensure_external_node(
            f"{module_name}.{member}",
            kind=self._guess_external_kind(member),
        )

    def resolve_method(self, class_id: str, method_name: str) -> str:
        """Resolve a method on a class node."""
        class_info = self.nodes.get(class_id)
        if not class_info:
            return self._ensure_external_node(
                f"{class_id}.{method_name}",
                kind="external-function",
            )

        if class_info.module not in self.modules:
            return self._ensure_external_node(
                f"{class_info.full_name}.{method_name}",
                kind="external-function",
            )

        module = self.modules[class_info.module]
        key = (class_info.name, method_name)
        if key in module.method_symbols:
            return module.method_symbols[key]

        return self._ensure_external_node(
            f"{class_info.full_name}.{method_name}",
            kind="external-function",
        )

    def is_class_node(self, node_id: str) -> bool:
        """Return whether the given node id denotes a class."""
        node = self.nodes.get(node_id)
        if node and node.kind == "class":
            return True
        return node_id.startswith("external:") and self.nodes[node_id].kind == "external-class"

    def add_definition(self, definition: DefinitionInfo) -> None:
        """Register a new definition node."""
        self.nodes[definition.node_id] = definition
        self.modules[definition.module].definitions[definition.node_id] = definition

    def add_top_level_symbol(self, module_name: str, name: str, node_id: str) -> None:
        """Register a top-level function or class name within a module."""
        module = self.modules[module_name]
        module.top_level_symbols[name] = node_id
        if self.nodes[node_id].kind == "class":
            module.class_symbols[name] = node_id

    def add_method_symbol(
        self,
        module_name: str,
        class_name: str,
        method_name: str,
        node_id: str,
    ) -> None:
        """Register a method name within a class."""
        self.modules[module_name].method_symbols[(class_name, method_name)] = node_id

    def add_import(self, module_name: str, binding: ImportBinding) -> None:
        """Register one import binding within a module."""
        self.modules[module_name].imports[binding.alias] = binding

    def _resolve_import_binding(self, binding: ImportBinding) -> str:
        """Resolve an import binding to a local node whenever possible."""
        if binding.symbol is None:
            return self._module_node_id(binding.module)
        return self.resolve_module_member(binding.module, binding.symbol)

    def _ensure_external_node(self, full_name: str, kind: str) -> str:
        """Create an external node if it does not exist yet."""
        node_id = f"external:{full_name}"
        if node_id in self.nodes:
            return node_id

        module = full_name.rsplit(".", 1)[0] if "." in full_name else full_name
        name = full_name.split(".")[-1]
        self.nodes[node_id] = DefinitionInfo(
            node_id=node_id,
            module=module,
            name=name,
            kind=kind,
            path=self.root_dir,
            lineno=0,
            end_lineno=0,
            doc="External dependency resolved from static analysis.",
            private=name.startswith("_"),
        )
        return node_id

    def _reachable_nodes(self, roots: list[str]) -> set[str]:
        """Return all nodes reachable from the selected entry points via calls/constructs."""
        adjacency = {}
        for edge in self.edges.values():
            if edge.edge_type not in {"calls", "constructs"}:
                continue
            adjacency.setdefault(edge.source, set()).add(edge.target)

        queue = [root for root in roots if root in self.nodes]
        seen = set(queue)

        while queue:
            current = queue.pop(0)
            for target in adjacency.get(current, set()):
                if target in seen:
                    continue
                seen.add(target)
                queue.append(target)

        return seen

    def _expand_with_structure(self, reachable: set[str]) -> set[str]:
        """Add modules, classes, and package re-export nodes needed for context."""
        included = set(reachable)

        for node_id in list(reachable):
            node = self.nodes[node_id]
            if node.kind == "method" and node.class_name:
                class_id = self.modules[node.module].class_symbols[node.class_name]
                included.add(class_id)
            if node.module in self.modules:
                included.add(self._module_node_id(node.module))

        changed = True
        while changed:
            changed = False
            for edge in self.edges.values():
                if edge.edge_type not in {"contains", "imports", "reexports"}:
                    continue
                if edge.target in included and edge.source not in included:
                    included.add(edge.source)
                    changed = True

        return included

    def _serialize_nodes(self, included: set[str], roots: list[str]) -> list[dict[str, object]]:
        """Create the node payload for the HTML viewer."""
        distances = self._call_distances(roots)
        serialized = []

        for node_id in sorted(included):
            if node_id.startswith("module:"):
                module_name = node_id.removeprefix("module:")
                module = self.module_nodes[module_name]
                lane = self._lane_key(module_name)
                color = LANE_STYLES[lane]["color"]
                serialized.append(
                    {
                        "id": node_id,
                        "label": module.name,
                        "fullName": module.name,
                        "module": module.name,
                        "kind": "module",
                        "lane": lane,
                        "color": color,
                        "path": self._relative_path(module.path),
                        "lineno": module.lineno,
                        "endLineno": module.end_lineno,
                        "doc": module.doc or "Module node.",
                        "private": False,
                        "external": False,
                        "distance": min(
                            (
                                distances[child.node_id]
                                for child in self.modules[module_name].definitions.values()
                                if child.node_id in distances
                            ),
                            default=0,
                        ),
                    }
                )
                continue

            node = self.nodes[node_id]
            lane = self._lane_key(node.module)
            color = LANE_STYLES[lane]["color"]
            serialized.append(
                {
                    "id": node_id,
                    "label": node.label,
                    "fullName": node.full_name,
                    "module": node.module,
                    "kind": node.kind,
                    "lane": lane,
                    "color": color,
                    "path": self._relative_path(node.path) if node.lineno else "",
                    "lineno": node.lineno,
                    "endLineno": node.end_lineno,
                    "doc": node.doc or "No docstring available.",
                    "private": node.private,
                    "external": node_id.startswith("external:"),
                    "distance": distances.get(node_id, 99),
                }
            )

        return serialized

    def _serialize_edges(self, included: set[str]) -> list[dict[str, object]]:
        """Create the edge payload for the HTML viewer."""
        serialized = []

        for edge in self.edges.values():
            if edge.source not in included or edge.target not in included:
                continue
            lines = sorted(set(edge.lines))
            calls = sorted(edge.calls)
            serialized.append(
                {
                    "source": edge.source,
                    "target": edge.target,
                    "type": edge.edge_type,
                    "label": ", ".join(calls[:4]),
                    "count": edge.count,
                    "lines": lines,
                }
            )

        return serialized

    def _layout_nodes(self, nodes: list[dict[str, object]]) -> list[dict[str, object]]:
        """Assign deterministic lane-based coordinates to the serialized nodes."""
        grouped: dict[str, list[dict[str, object]]] = {}
        for node in nodes:
            grouped.setdefault(node["lane"], []).append(node)

        lane_names = [lane for lane in PREFERRED_LANES if lane in grouped]
        lane_names.extend(sorted(set(grouped) - set(lane_names)))

        lane_width = 340
        lane_gap = 60
        top = 120
        row_height = 78
        margin = 60
        lanes = []

        for lane_index, lane_name in enumerate(lane_names):
            lane_nodes = grouped[lane_name]
            lane_nodes.sort(
                key=lambda node: (
                    str(node["module"]),
                    KIND_RANK.get(str(node["kind"]), 99),
                    int(node["distance"]),
                    str(node["label"]).lower(),
                )
            )

            x = margin + lane_index * (lane_width + lane_gap) + lane_width / 2
            for row_index, node in enumerate(lane_nodes):
                width = max(154, len(str(node["label"])) * 8 + 28)
                height = 52 if node["kind"] == "module" else 40
                node["x"] = x
                node["y"] = top + row_index * row_height
                node["width"] = width
                node["height"] = height

            lanes.append(
                {
                    "name": lane_name,
                    "label": lane_name,
                    "x": x - lane_width / 2,
                    "width": lane_width,
                    "color": LANE_STYLES[lane_name]["color"],
                    "band": LANE_STYLES[lane_name]["band"],
                }
            )

        return lanes

    def _call_distances(self, roots: list[str]) -> dict[str, int]:
        """Shortest path distances from the chosen roots using call/construct edges."""
        adjacency = {}
        for edge in self.edges.values():
            if edge.edge_type not in {"calls", "constructs"}:
                continue
            adjacency.setdefault(edge.source, set()).add(edge.target)

        queue = [(root, 0) for root in roots if root in self.nodes]
        distances = {root: 0 for root in roots if root in self.nodes}

        while queue:
            current, depth = queue.pop(0)
            for target in adjacency.get(current, set()):
                if target in distances and distances[target] <= depth + 1:
                    continue
                distances[target] = depth + 1
                queue.append((target, depth + 1))

        return distances

    def _register_edge(
        self,
        source: str,
        target: str,
        edge_type: str,
        raw_call: str,
        lineno: int | None,
    ) -> None:
        """Register or update one graph edge."""
        if source == target and edge_type == "contains":
            return

        key = (source, target, edge_type)
        if key not in self.edges:
            self.edges[key] = EdgeInfo(source=source, target=target, edge_type=edge_type)
        self.edges[key].register(raw_call, lineno)

    def _lane_key(self, module_name: str) -> str:
        """Return the lane name used for layout and coloring."""
        if module_name == "example_simulations" or module_name.startswith("example_simulations."):
            return "example_simulations"
        if module_name.startswith("asymmetree.hologenome"):
            return "asymmetree.hologenome"
        if module_name.startswith("asymmetree.holoevolve"):
            return "asymmetree.holoevolve"
        if module_name.startswith("asymmetree.treeevolve"):
            return "asymmetree.treeevolve"
        if module_name.startswith("asymmetree.utils"):
            return "asymmetree.utils"
        if module_name.startswith("asymmetree."):
            return "asymmetree.other"
        return "external"

    def _module_name(self, path: Path) -> str:
        """Convert a file path below the root directory into a module name."""
        rel = path.resolve().relative_to(self.root_dir)
        parts = list(rel.parts)
        if parts[-1] == "__init__.py":
            parts = parts[:-1]
        else:
            parts[-1] = Path(parts[-1]).stem
        return ".".join(parts)

    def _module_node_id(self, module_name: str) -> str:
        """Return the node id used for module nodes."""
        return f"module:{module_name}"

    def _relative_path(self, path: Path) -> str:
        """Return a path relative to the local example-code root."""
        try:
            return str(path.resolve().relative_to(self.root_dir))
        except ValueError:
            return str(path)

    def _guess_external_kind(self, symbol: str) -> str:
        """Best-effort distinction between imported classes and functions."""
        return "external-class" if symbol[:1].isupper() else "external-function"


class ModuleCollector(ast.NodeVisitor):
    """Collect imports and top-level definitions from a module."""

    def __init__(self, repository: RepositoryModel, module: ModuleInfo) -> None:
        self.repository = repository
        self.module = module

    def visit_Import(self, node: ast.Import) -> None:
        for alias in node.names:
            local_name = alias.asname or alias.name.split(".")[0]
            self.repository.add_import(
                self.module.name,
                ImportBinding(
                    alias=local_name,
                    module=alias.name,
                    symbol=None,
                    lineno=node.lineno,
                ),
            )

    def visit_ImportFrom(self, node: ast.ImportFrom) -> None:
        source_module = self._resolve_from_module(node.module, node.level)
        for alias in node.names:
            self.repository.add_import(
                self.module.name,
                ImportBinding(
                    alias=alias.asname or alias.name,
                    module=source_module,
                    symbol=alias.name,
                    lineno=node.lineno,
                ),
            )

    def visit_FunctionDef(self, node: ast.FunctionDef) -> None:
        definition = DefinitionInfo(
            node_id=f"{self.module.name}.{node.name}",
            module=self.module.name,
            name=node.name,
            kind="function",
            path=self.module.path,
            lineno=node.lineno,
            end_lineno=node.end_lineno or node.lineno,
            doc=ast.get_docstring(node) or "",
            ast_body=node.body,
            private=node.name.startswith("_"),
        )
        self.repository.add_definition(definition)
        self.repository.add_top_level_symbol(self.module.name, node.name, definition.node_id)

    def visit_ClassDef(self, node: ast.ClassDef) -> None:
        class_info = DefinitionInfo(
            node_id=f"{self.module.name}.{node.name}",
            module=self.module.name,
            name=node.name,
            kind="class",
            path=self.module.path,
            lineno=node.lineno,
            end_lineno=node.end_lineno or node.lineno,
            doc=ast.get_docstring(node) or "",
            private=node.name.startswith("_"),
        )
        self.repository.add_definition(class_info)
        self.repository.add_top_level_symbol(self.module.name, node.name, class_info.node_id)

        for statement in node.body:
            if not isinstance(statement, ast.FunctionDef):
                continue
            method_info = DefinitionInfo(
                node_id=f"{self.module.name}.{node.name}.{statement.name}",
                module=self.module.name,
                name=statement.name,
                kind="method",
                path=self.module.path,
                lineno=statement.lineno,
                end_lineno=statement.end_lineno or statement.lineno,
                doc=ast.get_docstring(statement) or "",
                class_name=node.name,
                ast_body=statement.body,
                private=statement.name.startswith("_"),
            )
            self.repository.add_definition(method_info)
            self.repository.add_method_symbol(
                self.module.name,
                node.name,
                statement.name,
                method_info.node_id,
            )

    def visit_If(self, node: ast.If) -> None:
        if self._is_main_guard(node):
            definition = DefinitionInfo(
                node_id=f"{self.module.name}.__main__",
                module=self.module.name,
                name="__main__",
                kind="main",
                path=self.module.path,
                lineno=node.lineno,
                end_lineno=node.end_lineno or node.lineno,
                doc="Module entry point guarded by ``if __name__ == '__main__'``.",
                ast_body=node.body,
            )
            self.repository.add_definition(definition)
            self.repository.add_top_level_symbol(self.module.name, "__main__", definition.node_id)
            return

        self.generic_visit(node)

    def _is_main_guard(self, node: ast.If) -> bool:
        """Return whether the ``if`` statement is the module's main guard."""
        test = node.test
        if not isinstance(test, ast.Compare) or len(test.ops) != 1 or len(test.comparators) != 1:
            return False
        if not isinstance(test.ops[0], ast.Eq):
            return False
        if not isinstance(test.left, ast.Name) or test.left.id != "__name__":
            return False
        comparator = test.comparators[0]
        if not isinstance(comparator, ast.Constant) or comparator.value != "__main__":
            return False
        return True

    def _resolve_from_module(self, module: str | None, level: int) -> str:
        """Resolve the absolute module name of a ``from ... import ...`` statement."""
        if level == 0:
            return module or ""

        if self.module.is_package:
            package = self.module.name
        elif "." in self.module.name:
            package = self.module.name.rsplit(".", 1)[0]
        else:
            package = ""

        package_parts = package.split(".") if package else []
        keep = max(0, len(package_parts) - (level - 1))
        prefix = package_parts[:keep]
        if module:
            prefix.extend(module.split("."))
        return ".".join(part for part in prefix if part)


class BodyCallAnalyzer(ast.NodeVisitor):
    """Collect call edges within one function or method body."""

    def __init__(self, repository: RepositoryModel, definition: DefinitionInfo) -> None:
        self.repository = repository
        self.definition = definition
        self.local_bindings: dict[str, str] = {}

    def visit_FunctionDef(self, node: ast.FunctionDef) -> None:
        return

    def visit_ClassDef(self, node: ast.ClassDef) -> None:
        return

    def visit_Assign(self, node: ast.Assign) -> None:
        self.visit(node.value)
        bound_class = self._resolve_receiver(node.value)
        if bound_class and bound_class[0] == "class":
            for target in node.targets:
                for name in self._assigned_names(target):
                    self.local_bindings[name] = bound_class[1]

    def visit_AnnAssign(self, node: ast.AnnAssign) -> None:
        if node.value is not None:
            self.visit(node.value)
            bound_class = self._resolve_receiver(node.value)
            if bound_class and bound_class[0] == "class":
                for name in self._assigned_names(node.target):
                    self.local_bindings[name] = bound_class[1]

    def visit_Call(self, node: ast.Call) -> None:
        target = self._resolve_call(node.func)
        if target:
            edge_type = "constructs" if self.repository.is_class_node(target) else "calls"
            self.repository._register_edge(
                self.definition.node_id,
                target,
                edge_type,
                raw_call=self._call_text(node.func),
                lineno=node.lineno,
            )
        self.generic_visit(node)

    def _resolve_call(self, func: ast.expr) -> str | None:
        """Resolve the callee of one call expression."""
        if isinstance(func, ast.Name):
            return self._resolve_name(func.id)
        if isinstance(func, ast.Attribute):
            return self._resolve_attribute(func)
        return None

    def _resolve_name(self, name: str) -> str | None:
        """Resolve a bare name in the current scope."""
        if name in self.local_bindings:
            return self.local_bindings[name]
        return self.repository.resolve_name(self.definition.module, name)

    def _resolve_attribute(self, node: ast.Attribute) -> str | None:
        """Resolve an attribute access used as a callee."""
        receiver = self._resolve_receiver(node.value)
        if not receiver:
            return None

        receiver_kind, receiver_id = receiver
        if receiver_kind == "module":
            module_name = receiver_id.removeprefix("module:")
            return self.repository.resolve_module_member(module_name, node.attr)
        if receiver_kind == "class":
            return self.repository.resolve_method(receiver_id, node.attr)
        return None

    def _resolve_receiver(self, expr: ast.expr) -> tuple[str, str] | None:
        """Resolve a call receiver to either a module or a class."""
        if isinstance(expr, ast.Name):
            if expr.id == "self" and self.definition.class_name:
                class_id = self.repository.modules[self.definition.module].class_symbols[
                    self.definition.class_name
                ]
                return ("class", class_id)
            if expr.id in self.local_bindings:
                return ("class", self.local_bindings[expr.id])
            module = self.repository.modules[self.definition.module]
            if expr.id in module.imports and module.imports[expr.id].symbol is None:
                target_module = module.imports[expr.id].module
                if target_module in self.repository.modules:
                    return ("module", self.repository._module_node_id(target_module))
                return ("module", self.repository._module_node_id(target_module))
            return None

        if isinstance(expr, ast.Call):
            target = self._resolve_call(expr.func)
            if target and self.repository.is_class_node(target):
                return ("class", target)
            return None

        if isinstance(expr, ast.Attribute):
            target = self._resolve_attribute(expr)
            if not target:
                return None
            if target.startswith("module:"):
                return ("module", target)
            if self.repository.is_class_node(target):
                return ("class", target)
            return None

        return None

    def _assigned_names(self, target: ast.expr) -> list[str]:
        """Collect variable names assigned by a target expression."""
        if isinstance(target, ast.Name):
            return [target.id]
        if isinstance(target, (ast.Tuple, ast.List)):
            names = []
            for item in target.elts:
                names.extend(self._assigned_names(item))
            return names
        return []

    def _call_text(self, expr: ast.expr) -> str:
        """Best-effort textual representation of a callee expression."""
        if isinstance(expr, ast.Name):
            return expr.id
        if isinstance(expr, ast.Attribute):
            return f"{self._call_text(expr.value)}.{expr.attr}"
        return "<call>"


def render_html(graph: dict[str, object]) -> str:
    """Return the self-contained HTML viewer."""
    payload = json.dumps(graph, indent=2)
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Dependency Diagram</title>
  <style>
    :root {{
      --bg: #f7f4ea;
      --panel: #fffdf8;
      --ink: #172033;
      --muted: #5f6674;
      --line: #d8d4c5;
      --accent: #d4a72c;
      --shadow: 0 18px 40px rgba(23, 32, 51, 0.12);
    }}
    * {{
      box-sizing: border-box;
    }}
    body {{
      margin: 0;
      font-family: "IBM Plex Sans", "Avenir Next", "Segoe UI", sans-serif;
      color: var(--ink);
      background:
        radial-gradient(circle at top left, #fff7db 0, transparent 32%),
        linear-gradient(180deg, #faf7ef 0%, #f2ede0 100%);
    }}
    .app {{
      display: grid;
      grid-template-columns: 360px 1fr;
      min-height: 100vh;
    }}
    .sidebar {{
      padding: 24px 22px 18px;
      background: rgba(255, 253, 248, 0.94);
      border-right: 1px solid var(--line);
      box-shadow: var(--shadow);
      position: relative;
      z-index: 2;
      overflow: auto;
    }}
    .sidebar h1 {{
      margin: 0;
      font-size: 1.55rem;
      letter-spacing: -0.03em;
    }}
    .lead {{
      margin: 8px 0 0;
      color: var(--muted);
      line-height: 1.45;
      font-size: 0.95rem;
    }}
    .controls {{
      margin-top: 20px;
      display: grid;
      gap: 10px;
    }}
    .controls input[type="search"] {{
      width: 100%;
      padding: 11px 12px;
      border: 1px solid var(--line);
      border-radius: 10px;
      background: #fff;
      font: inherit;
    }}
    .controls label {{
      display: flex;
      align-items: center;
      gap: 10px;
      font-size: 0.94rem;
      color: var(--ink);
    }}
    .button-row {{
      display: flex;
      gap: 10px;
      flex-wrap: wrap;
    }}
    button {{
      border: 1px solid #c9bc92;
      background: linear-gradient(180deg, #f7e7b0, #efd27e);
      color: #493300;
      padding: 9px 12px;
      border-radius: 999px;
      cursor: pointer;
      font: inherit;
      font-weight: 600;
    }}
    .legend {{
      margin-top: 20px;
      display: grid;
      gap: 8px;
    }}
    .legend-row {{
      display: flex;
      align-items: center;
      gap: 10px;
      font-size: 0.93rem;
      color: var(--muted);
    }}
    .swatch {{
      width: 13px;
      height: 13px;
      border-radius: 999px;
      border: 1px solid rgba(23, 32, 51, 0.15);
      flex: 0 0 auto;
    }}
    .stats, .details {{
      margin-top: 20px;
      padding: 14px;
      border: 1px solid var(--line);
      border-radius: 14px;
      background: #fff;
      box-shadow: 0 10px 24px rgba(23, 32, 51, 0.06);
    }}
    .details h2 {{
      margin: 0 0 8px;
      font-size: 1rem;
    }}
    .meta {{
      margin: 0;
      padding: 0;
      list-style: none;
      display: grid;
      gap: 6px;
      color: var(--muted);
      font-size: 0.92rem;
    }}
    .details p {{
      margin: 10px 0 0;
      line-height: 1.5;
      font-size: 0.94rem;
    }}
    .list {{
      margin: 12px 0 0;
      padding-left: 18px;
      color: var(--muted);
      font-size: 0.92rem;
    }}
    .stage {{
      position: relative;
      overflow: hidden;
    }}
    svg {{
      width: 100%;
      height: 100vh;
      display: block;
      cursor: grab;
    }}
    svg.dragging {{
      cursor: grabbing;
    }}
    .lane-band {{
      rx: 28px;
    }}
    .lane-label {{
      font-size: 16px;
      font-weight: 700;
      fill: #41506a;
    }}
    .edge {{
      fill: none;
      stroke-width: 2.1;
      opacity: 0.72;
    }}
    .edge.calls {{
      stroke: #394354;
    }}
    .edge.constructs {{
      stroke: #8c4f00;
      stroke-dasharray: 8 6;
    }}
    .edge.contains,
    .edge.imports,
    .edge.reexports {{
      stroke: #b8b4a6;
      stroke-dasharray: 5 7;
      opacity: 0.55;
    }}
    .node-rect {{
      stroke-width: 1.5;
      filter: drop-shadow(0 10px 18px rgba(23, 32, 51, 0.08));
    }}
    .node-label {{
      font-size: 13px;
      text-anchor: middle;
      dominant-baseline: central;
      pointer-events: none;
      fill: #172033;
      font-weight: 600;
    }}
    .node.module .node-label {{
      font-size: 14px;
      font-weight: 700;
    }}
    .dimmed {{
      opacity: 0.13 !important;
    }}
    .hidden {{
      display: none;
    }}
    .selected .node-rect {{
      stroke: #172033;
      stroke-width: 2.6;
    }}
    .selected-edge {{
      stroke-width: 3.4;
      opacity: 0.95;
    }}
    .hint {{
      color: var(--muted);
      font-size: 0.9rem;
      line-height: 1.45;
    }}
    @media (max-width: 980px) {{
      .app {{
        grid-template-columns: 1fr;
      }}
      .sidebar {{
        border-right: none;
        border-bottom: 1px solid var(--line);
      }}
      svg {{
        height: 72vh;
      }}
    }}
  </style>
</head>
<body>
  <div class="app">
    <aside class="sidebar">
      <h1>Dependency Diagram</h1>
      <p class="lead">
        Static call/dependency view rooted in <code>{graph["entry_module"]}.py</code>. The default
        view emphasizes the local script plus the <code>asymmetree.hologenome</code> and
        <code>asymmetree.holoevolve</code> paths it triggers.
      </p>

      <div class="controls">
        <input id="search" type="search" placeholder="Search modules, functions, methods...">
        <label><input id="showExternal" type="checkbox"> Show external dependencies</label>
        <label><input id="showPrivate" type="checkbox" checked> Show private helpers</label>
        <label><input id="showStructure" type="checkbox" checked> Show contains / import edges</label>
        <div class="button-row">
          <button id="resetView" type="button">Reset View</button>
          <button id="clearSelection" type="button">Clear Selection</button>
        </div>
      </div>

      <div class="legend" id="legend"></div>
      <div class="stats" id="stats"></div>
      <div class="details" id="details">
        <h2>Selection</h2>
        <p class="hint">Click a node to inspect its callers, callees, file location, and docstring.</p>
      </div>
    </aside>

    <main class="stage">
      <svg id="graph" viewBox="0 0 1800 1200" aria-label="Dependency diagram">
        <g id="viewport"></g>
      </svg>
    </main>
  </div>

  <script>
    const graph = {payload};

    const svg = document.getElementById("graph");
    const viewport = document.getElementById("viewport");
    const searchInput = document.getElementById("search");
    const showExternal = document.getElementById("showExternal");
    const showPrivate = document.getElementById("showPrivate");
    const showStructure = document.getElementById("showStructure");
    const resetViewButton = document.getElementById("resetView");
    const clearSelectionButton = document.getElementById("clearSelection");
    const statsPanel = document.getElementById("stats");
    const detailsPanel = document.getElementById("details");
    const legend = document.getElementById("legend");

    const laneMap = new Map(graph.lanes.map((lane) => [lane.name, lane]));
    const nodeMap = new Map(graph.nodes.map((node) => [node.id, node]));
    const outgoing = new Map();
    const incoming = new Map();

    for (const edge of graph.edges) {{
      if (!outgoing.has(edge.source)) {{
        outgoing.set(edge.source, []);
      }}
      if (!incoming.has(edge.target)) {{
        incoming.set(edge.target, []);
      }}
      outgoing.get(edge.source).push(edge);
      incoming.get(edge.target).push(edge);
    }}

    let selectedId = graph.defaultSelection || null;
    let transform = {{
      x: 24,
      y: 24,
      scale: 1
    }};
    let isDragging = false;
    let dragOrigin = null;

    function updateLegend() {{
      legend.innerHTML = "";
      for (const lane of graph.lanes) {{
        const row = document.createElement("div");
        row.className = "legend-row";
        row.innerHTML = `
          <span class="swatch" style="background:${{lane.color}}"></span>
          <span>${{lane.label}}</span>
        `;
        legend.appendChild(row);
      }}
    }}

    function passesFilters(node, query) {{
      if (!showExternal.checked && node.external) {{
        return false;
      }}
      if (!showPrivate.checked && node.private) {{
        return false;
      }}
      if (!query) {{
        return true;
      }}
      const haystack = `${{node.label}} ${{node.fullName}} ${{node.module}} ${{node.doc}}`.toLowerCase();
      return haystack.includes(query);
    }}

    function render() {{
      const query = searchInput.value.trim().toLowerCase();
      const visibleNodes = new Map();
      for (const node of graph.nodes) {{
        if (passesFilters(node, query)) {{
          visibleNodes.set(node.id, node);
        }}
      }}

      const visibleEdges = graph.edges.filter((edge) => {{
        if (!visibleNodes.has(edge.source) || !visibleNodes.has(edge.target)) {{
          return false;
        }}
        if (!showStructure.checked && ["contains", "imports", "reexports"].includes(edge.type)) {{
          return false;
        }}
        return true;
      }});

      const highlight = new Set();
      const highlightEdges = new Set();
      if (selectedId && visibleNodes.has(selectedId)) {{
        highlight.add(selectedId);
        for (const edge of outgoing.get(selectedId) || []) {{
          if (visibleNodes.has(edge.target)) {{
            highlight.add(edge.target);
            highlightEdges.add(edge);
          }}
        }}
        for (const edge of incoming.get(selectedId) || []) {{
          if (visibleNodes.has(edge.source)) {{
            highlight.add(edge.source);
            highlightEdges.add(edge);
          }}
        }}
      }}

      const maxX = Math.max(...graph.lanes.map((lane) => lane.x + lane.width + 36), 1600);
      const maxY = Math.max(...graph.nodes.map((node) => node.y + node.height + 96), 1000);
      svg.setAttribute("viewBox", `0 0 ${{maxX}} ${{maxY}}`);
      viewport.innerHTML = "";

      for (const lane of graph.lanes) {{
        const rect = document.createElementNS("http://www.w3.org/2000/svg", "rect");
        rect.setAttribute("x", lane.x);
        rect.setAttribute("y", 36);
        rect.setAttribute("width", lane.width);
        rect.setAttribute("height", maxY - 72);
        rect.setAttribute("fill", lane.band);
        rect.setAttribute("class", "lane-band");
        viewport.appendChild(rect);

        const text = document.createElementNS("http://www.w3.org/2000/svg", "text");
        text.setAttribute("x", lane.x + 18);
        text.setAttribute("y", 62);
        text.setAttribute("class", "lane-label");
        text.textContent = lane.label;
        viewport.appendChild(text);
      }}

      for (const edge of visibleEdges) {{
        const source = visibleNodes.get(edge.source);
        const target = visibleNodes.get(edge.target);
        const sameLane = source.lane === target.lane;
        const startX = source.x + source.width / 2;
        const endX = target.x - target.width / 2;
        const startY = source.y;
        const endY = target.y;
        const bend = sameLane ? 120 : Math.max(80, Math.abs(endX - startX) * 0.24);
        const controlX = sameLane ? startX + bend : (startX + endX) / 2;
        const controlY = sameLane ? (startY + endY) / 2 : (startY + endY) / 2;

        const path = document.createElementNS("http://www.w3.org/2000/svg", "path");
        path.setAttribute(
          "d",
          `M ${{startX}} ${{startY}} Q ${{controlX}} ${{controlY}} ${{endX}} ${{endY}}`
        );
        path.setAttribute("class", `edge ${{edge.type}}`);
        if (selectedId && highlight.size > 0 && !highlightEdges.has(edge)) {{
          path.classList.add("dimmed");
        }}
        if (highlightEdges.has(edge)) {{
          path.classList.add("selected-edge");
        }}
        const title = document.createElementNS("http://www.w3.org/2000/svg", "title");
        const linePart = edge.lines.length ? ` line${{edge.lines.length > 1 ? "s" : ""}} ${{edge.lines.join(", ")}}` : "";
        title.textContent = `${{edge.type}}: ${{edge.label || "relation"}}${{linePart}}`;
        path.appendChild(title);
        viewport.appendChild(path);
      }}

      for (const node of visibleNodes.values()) {{
        const group = document.createElementNS("http://www.w3.org/2000/svg", "g");
        group.setAttribute("class", `node ${{node.kind}}`);
        group.dataset.id = node.id;
        if (selectedId === node.id) {{
          group.classList.add("selected");
        }} else if (selectedId && highlight.size > 0 && !highlight.has(node.id)) {{
          group.classList.add("dimmed");
        }}

        const rect = document.createElementNS("http://www.w3.org/2000/svg", "rect");
        rect.setAttribute("x", node.x - node.width / 2);
        rect.setAttribute("y", node.y - node.height / 2);
        rect.setAttribute("width", node.width);
        rect.setAttribute("height", node.height);
        rect.setAttribute("rx", node.kind === "module" ? 18 : 14);
        rect.setAttribute("fill", node.kind === "module" ? node.color : "white");
        rect.setAttribute("stroke", node.color);
        rect.setAttribute("class", "node-rect");
        group.appendChild(rect);

        const label = document.createElementNS("http://www.w3.org/2000/svg", "text");
        label.setAttribute("x", node.x);
        label.setAttribute("y", node.y);
        label.setAttribute("class", "node-label");
        label.textContent = node.label;
        group.appendChild(label);

        const title = document.createElementNS("http://www.w3.org/2000/svg", "title");
        title.textContent = `${{node.fullName}}${{node.lineno ? ` (line ${{node.lineno}})` : ""}}`;
        group.appendChild(title);

        group.addEventListener("click", (event) => {{
          event.stopPropagation();
          selectedId = node.id;
          updateDetails(node, visibleNodes);
          render();
        }});
        viewport.appendChild(group);
      }}

      updateStats(visibleNodes.size, visibleEdges.length);
      if (selectedId && visibleNodes.has(selectedId)) {{
        updateDetails(visibleNodes.get(selectedId), visibleNodes);
      }} else {{
        renderHint();
      }}

      applyTransform();
    }}

    function updateStats(nodeCount, edgeCount) {{
      statsPanel.innerHTML = `
        <strong>${{nodeCount}}</strong> visible nodes<br>
        <strong>${{edgeCount}}</strong> visible edges<br>
        <span class="hint">Calls and constructors define reachability. Import and containment edges provide structure.</span>
      `;
    }}

    function updateDetails(node, visibleNodes) {{
      const callers = (incoming.get(node.id) || [])
        .filter((edge) => visibleNodes.has(edge.source))
        .map((edge) => `${{nodeMap.get(edge.source).label}} [${{edge.type}}]`);
      const callees = (outgoing.get(node.id) || [])
        .filter((edge) => visibleNodes.has(edge.target))
        .map((edge) => `${{nodeMap.get(edge.target).label}} [${{edge.type}}]`);
      const fileLine = node.path && node.lineno ? `${{node.path}}:${{node.lineno}}` : "external / inferred";

      detailsPanel.innerHTML = `
        <h2>${{node.label}}</h2>
        <ul class="meta">
          <li><strong>Kind:</strong> ${{node.kind}}</li>
          <li><strong>Module:</strong> ${{node.module}}</li>
          <li><strong>Source:</strong> ${{fileLine}}</li>
          <li><strong>Full name:</strong> ${{node.fullName}}</li>
        </ul>
        <p>${{node.doc || "No docstring available."}}</p>
        <p><strong>Callers / parents</strong></p>
        <ul class="list">
          ${{callers.length ? callers.map((entry) => `<li>${{entry}}</li>`).join("") : "<li>None visible</li>"}}
        </ul>
        <p><strong>Callees / children</strong></p>
        <ul class="list">
          ${{callees.length ? callees.map((entry) => `<li>${{entry}}</li>`).join("") : "<li>None visible</li>"}}
        </ul>
      `;
    }}

    function renderHint() {{
      detailsPanel.innerHTML = `
        <h2>Selection</h2>
        <p class="hint">
          Click a node to inspect its neighborhood. Search narrows the graph; hiding external nodes
          is useful when you want to concentrate on the script, <code>asymmetree.hologenome</code>,
          and <code>asymmetree.holoevolve</code>.
        </p>
      `;
    }}

    function applyTransform() {{
      viewport.setAttribute(
        "transform",
        `translate(${{transform.x}}, ${{transform.y}}) scale(${{transform.scale}})`
      );
    }}

    svg.addEventListener("click", () => {{
      selectedId = null;
      render();
    }});

    svg.addEventListener("pointerdown", (event) => {{
      if (event.target.closest(".node")) {{
        return;
      }}
      isDragging = true;
      dragOrigin = {{ x: event.clientX, y: event.clientY, tx: transform.x, ty: transform.y }};
      svg.classList.add("dragging");
    }});

    window.addEventListener("pointermove", (event) => {{
      if (!isDragging || !dragOrigin) {{
        return;
      }}
      transform.x = dragOrigin.tx + (event.clientX - dragOrigin.x);
      transform.y = dragOrigin.ty + (event.clientY - dragOrigin.y);
      applyTransform();
    }});

    window.addEventListener("pointerup", () => {{
      isDragging = false;
      dragOrigin = null;
      svg.classList.remove("dragging");
    }});

    svg.addEventListener("wheel", (event) => {{
      event.preventDefault();
      const direction = event.deltaY < 0 ? 1.08 : 0.92;
      transform.scale = Math.max(0.4, Math.min(2.8, transform.scale * direction));
      applyTransform();
    }}, {{ passive: false }});

    resetViewButton.addEventListener("click", () => {{
      transform = {{ x: 24, y: 24, scale: 1 }};
      applyTransform();
    }});

    clearSelectionButton.addEventListener("click", () => {{
      selectedId = null;
      render();
    }});

    searchInput.addEventListener("input", render);
    showExternal.addEventListener("change", render);
    showPrivate.addEventListener("change", render);
    showStructure.addEventListener("change", render);

    updateLegend();
    render();
  </script>
</body>
</html>
"""


def module_name_from_entry(root_dir: Path, entry: Path) -> str:
    """Return the module name corresponding to the selected entry file."""
    rel = entry.resolve().relative_to(root_dir.resolve())
    if rel.name == "__init__.py":
        parts = rel.parts[:-1]
    else:
        parts = rel.parts[:-1] + (rel.stem,)
    return ".".join(parts)


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--entry",
        type=Path,
        default=Path(DEFAULT_ENTRY),
        help="Entry script or module file to analyze.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(DEFAULT_OUTPUT),
        help="Output HTML file.",
    )
    return parser.parse_args()


def main() -> None:
    """Generate the HTML dependency diagram."""
    args = parse_args()
    root_dir = Path(__file__).resolve().parent
    entry_path = args.entry if args.entry.is_absolute() else root_dir / args.entry
    output_path = args.output if args.output.is_absolute() else root_dir / args.output

    repository = RepositoryModel(root_dir)
    entry_module = module_name_from_entry(root_dir, entry_path)
    html = render_html(repository.build_diagram(entry_module))
    output_path.write_text(html, encoding="utf-8")

    print(f"Wrote interactive diagram to {output_path}")


if __name__ == "__main__":
    main()
