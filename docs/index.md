# Python Template

Welcome to the **python-template** documentation.

## Overview

This is a Python project template providing a consistent, feature-rich starting point for new Python packages. It includes:

- **[Hatch](https://hatch.pypa.io/)** for build and environment management
- **[Ruff](https://docs.astral.sh/ruff/)** for linting and formatting
- **[pytest](https://docs.pytest.org/)** with coverage reporting
- **[mypy](https://mypy.readthedocs.io/)** for static type checking
- **[MkDocs](https://www.mkdocs.org/)** with Material theme for documentation
- **[Docker](https://www.docker.com/)** multi-stage builds (dev, hatch, prod)
- **[pre-commit](https://pre-commit.com/)** hooks for consistent commits
- **GitHub Actions** CI/CD workflows
- **Dev Container** support for VS Code / GitHub Codespaces

## Quick Start

### Installation

```console
pip install python-template
```

### Development Setup

Clone the repository and install [Hatch](https://hatch.pypa.io/latest/install/):

```console
git clone https://github.com/LavLabInfrastructure/python-template.git
cd python-template
pip install hatch
```

### Common Commands

| Task | Command |
|------|---------|
| Run tests | `hatch run test:test` |
| Run tests with coverage | `hatch run test:cov` |
| Lint code | `hatch run lint:check` |
| Format code | `hatch run lint:format` |
| Auto-fix lint issues | `hatch run lint:fix` |
| Type check | `hatch run types:check` |
| Build docs | `hatch run docs:build-docs` |
| Serve docs locally | `hatch run docs:serve-docs` |
| Build wheel | `hatch build` |

### CLI Usage

The template includes a CLI entry point pattern:

```console
python-template --version
python-template example "hello world"
```

## Project Structure

```text
python-template/
├── src/
│   └── python_template/      # Package source code
│       ├── __init__.py        # Public API exports
│       ├── __about__.py       # Version info
│       ├── cli.py             # CLI entry point (wraps library functions)
│       ├── example.py         # Example library module
│       └── py.typed           # PEP 561 type marker
├── tests/                     # Test suite
│   ├── conftest.py            # Shared fixtures
│   └── test_example.py        # Example tests
├── docs/                      # Documentation source
├── requirements/              # Locked dependency files (auto-generated)
├── .devcontainer/             # Dev container configuration
├── .github/                   # GitHub Actions & Dependabot
├── pyproject.toml             # Project metadata & tool config
├── Dockerfile                 # Multi-stage Docker build
├── Makefile                   # Common task shortcuts
├── mkdocs.yml                 # Documentation config
└── .pre-commit-config.yaml    # Pre-commit hooks
```

## Design Philosophy

This template follows a **library-first** approach:

1. All business logic lives in library modules under `src/python_template/`
2. The CLI (`cli.py`) is a thin wrapper that parses arguments and delegates to library functions
3. This makes the code importable, testable, and reusable — whether invoked from the command line, another Python package, or a notebook

## API Reference

::: python_template