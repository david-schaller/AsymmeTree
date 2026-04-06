# Installation

AsymmeTree requires **Python 3.10 or higher**.

## From PyPI

```bash
pip install asymmetree
```

## From source

```bash
git clone https://github.com/david-schaller/AsymmeTree.git
cd AsymmeTree
pip install .
```

## Contributors

### uv

If you want to contribute to `AsymmeTree`, please use the package and project manager
[uv](https://docs.astral.sh/uv/).
See [this page](https://docs.astral.sh/uv/getting-started/installation/) for installation
instructions.

To set up `uv` for your local `AsymmeTree` repository, navigate to the root directory of your
local `AsymmeTree` repository create a new virtual environment that is managed by `uv`:

```bash
cd <MY_PATH_TO>/AsymmeTree
uv sync
```

A new virtual environment will be created in the `.venv` directory of your local `AsymmeTree`
repository, and all dependencies of `AsymmeTree` will be installed in this virtual environment.
To activate this virtual environment, run the following command:

```bash
source .venv/bin/activate
```

### pre-commit

Please use [pre-commit](https://pre-commit.com) for automated code formatting and linting.
To install it and initialize it for your local `AsymmeTree` repository, follow these steps:

Install `uv` as described in the previous section.
Run the following command (after which you should be able to run `pre-commit` from anywhere):

```bash
uv tool install pre-commit --with pre-commit-uv
```

Navigate to the root directory of your local `AsymmeTree` repository and install `pre-commit`
as a git hook in the `AsymmeTree` repository:

```bash
cd <MY_PATH_TO>/AsymmeTree
pre-commit install
```

After this, `pre-commit` will automatically run the configured hooks (e.g., code formatting and
linting) before each commit, ensuring that your code adheres to the project's coding standards.

### MkDocs

The documentation of `AsymmeTree` is built using [MkDocs](https://www.mkdocs.org/).
To install it, run the following command:

```bash
cd <MY_PATH_TO>/AsymmeTree
uv sync --group docs

# or if you want to install all dependency groups, including the `docs` group:
uv sync --all-groups
```

After installing the dependencies, you can serve the documentation locally by running the following
command in the root directory of your local `AsymmeTree` repository:

```bash
uv run mkdocs serve
```

This will start a local development server, and you can access the documentation by navigating to
URL displayed in the terminal.

The documentation will automatically be deployed to GitHub Pages when pushing or merging to the
`main` or `dev` branch:

| Branch | URL |
| --- | --- |
| `main` | https://david-schaller.github.io/AsymmeTree/ |
| `dev` | https://david-schaller.github.io/AsymmeTree/dev/ |
