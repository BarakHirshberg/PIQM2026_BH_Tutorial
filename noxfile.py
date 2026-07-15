"""Task runner for the bosons/fermions PIMD tutorial.

Because every dependency is available from PyPI (i-PI ships its own driver, so
there is no compiled code and no conda-only package), nox can build a plain
virtual environment -- no conda required for the automated path.

    nox -e bosons-fermions-pimd   # install deps in a venv and run the recipe
    nox -e lint                   # check formatting/style
    nox -e format                 # auto-fix formatting/style
"""

import nox


nox.options.reuse_existing_virtualenvs = True
nox.options.sessions = ["lint"]

EXAMPLE = "examples/bosons-fermions-pimd"
RECIPE_DEPS = ["ipi>=3.2.0", "numpy", "matplotlib", "ase", "chemiscope>=1.0"]
LINT_TARGETS = [f"{EXAMPLE}/bosons-fermions-pimd.py", f"{EXAMPLE}/analysis.py", "noxfile.py"]


@nox.session(name="bosons-fermions-pimd")
def bosons_fermions_pimd(session):
    """Install dependencies and execute the recipe end-to-end."""
    session.install(*RECIPE_DEPS)
    session.chdir(EXAMPLE)
    session.run("python", "bosons-fermions-pimd.py", env={"MPLBACKEND": "Agg"})


@nox.session
def lint(session):
    """Check code style with ruff."""
    session.install("ruff")
    session.run("ruff", "check", *LINT_TARGETS)
    session.run("ruff", "format", "--check", *LINT_TARGETS)


@nox.session
def format(session):
    """Auto-format the code with ruff."""
    session.install("ruff")
    session.run("ruff", "check", "--fix", *LINT_TARGETS)
    session.run("ruff", "format", *LINT_TARGETS)
