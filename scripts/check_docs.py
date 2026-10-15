"""Check that the docs cover the code (plan Task 16).

- Every function and class in physics.py, equations.py and boundary.py appears in
  docs/physics.md on a line that gives a paper equation ("eq.") or section ("§").
- Every module in src/outflow/ and app/ appears in docs/architecture.md.

Run from the repository root: ``python scripts/check_docs.py``. Exits with 1 on a gap.
"""

import ast
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def defined_names(path: Path) -> list[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    names = []
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef | ast.ClassDef) and not node.name.startswith("__"):
            names.append(node.name)
    return names


def main() -> int:
    problems = []
    physics_lines = (ROOT / "docs" / "physics.md").read_text(encoding="utf-8").splitlines()
    for module in ("physics.py", "equations.py", "boundary.py"):
        for name in defined_names(ROOT / "src" / "outflow" / module):
            if name == "derivative":  # local closure inside a trial integration
                continue
            rows = [line for line in physics_lines if f"`{name}`" in line or f".{name}`" in line]
            if not any("eq" in row or "§" in row for row in rows):
                problems.append(f"physics.md: {module}:{name} missing or without an equation")
    architecture = (ROOT / "docs" / "architecture.md").read_text(encoding="utf-8")
    modules = sorted((ROOT / "src" / "outflow").glob("*.py")) + sorted((ROOT / "app").glob("*.py"))
    for path in modules:
        if f"`{path.name}`" not in architecture:
            problems.append(f"architecture.md: module {path.name} missing")
    for problem in problems:
        print(problem)
    print("docs check:", "OK" if not problems else f"{len(problems)} gaps")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
