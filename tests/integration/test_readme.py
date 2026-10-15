"""The README's quickstart code runs as written (its first ```python block)."""

import re
from pathlib import Path

README = Path(__file__).resolve().parents[2] / "README.md"


def first_python_block(text: str) -> str:
    match = re.search(r"```python\n(.*?)```", text, re.S)
    assert match, "README.md has no ```python block"
    return match.group(1)


def test_readme_quickstart_runs(tmp_path, monkeypatch):
    code = first_python_block(README.read_text(encoding="utf-8"))
    assert "outer_radius=5 * params.starburst_radius" in code  # keeps the test fast
    monkeypatch.chdir(tmp_path)  # the quickstart writes wind.csv
    namespace: dict = {}
    exec(compile(code, "README.md", "exec"), namespace)
    profile = namespace["profile"]
    assert profile.r[-1] >= 5 * namespace["params"].starburst_radius
    assert (tmp_path / "wind.csv").exists()
