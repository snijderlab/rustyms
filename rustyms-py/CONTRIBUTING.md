# Contributing

## Before you begin

If you have an idea for a feature, use case to add or an approach for a bugfix, you are welcome to communicate it with the community by opening a thread in GitHub Discussions or in GitHub Issues.

Not sure where to start? Check out the rustyms Rust code and see which parts are missing in the Python bindings. Also check out the open issues that carry the `good first issue` or `help wanted` labels.

## Development setup

Both Rust and Python are required:

- [Rust](https://rustup.rs/)
- [Python](https://github.com/pyenv/pyenv?tab=readme-ov-file#getting-pyenv)

Ideally, setup a virtual environment for Python:

```bash
python -m venv .venv
source .venv/bin/activate
```

PyO3 and Maturin are used to build the Rust extensions for Python.

Clone the repository:

```bash
git clone https://github.com/snijderlab/rustyms
cd rustyms
```

Either install Maturin and use it to build and install the Python bindings (recommended):

```bash
pip install maturin
maturin develop -m ./rustyms-py/Cargo.toml
```

Or use pip:

```bash
pip install --editable ./rustyms-py
```

## Documentation

The Python documentation is built with Sphinx and hosted on Read the Docs. To build the documentation locally, install the `docs` dependencies and run `sphinx-autobuild`:

```bash
pip install -e ./rustyms-py[docs]
sphinx-autobuild docs/python/source/ docs/python/_build
```
