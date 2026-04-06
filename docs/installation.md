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

If you want to contribute to `asymmetree`, please use the package and project manager
[uv](https://docs.astral.sh/uv/).
See [this page](https://docs.astral.sh/uv/getting-started/installation/) for installation
instructions.

Moreover, please use [pre-commit](https://pre-commit.com) for automated code formatting and linting.

To install it and initialize it for your local `asymmetree` repository, follow these steps:

- Install `uv`
- Run the following command (after which you should be able to run `pre-commit` from anywhere)
    - `uv tool install pre-commit --with pre-commit-uv`
- Navigate to the root directory of your local `asymmetree` repository
    - `cd <MY_PATH_TO>/asymmetree`
- Install `pre-commit` as a git hook in the `asymmetree` repository
    - `pre-commit install`
