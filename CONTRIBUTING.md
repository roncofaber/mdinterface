# Contributing to mdinterface

Thank you for your interest in contributing to mdinterface! This document provides guidelines and instructions for contributing to the project.

## Development Setup

### Prerequisites

- Python 3.8 or higher
- Git
- PACKMOL (for molecular packing functionality)

### Setting up the development environment

1. Fork the repository on GitLab
2. Clone your fork locally:
   ```bash
   git clone https://gitlab.com/yourusername/mdinterface.git
   cd mdinterface
   ```

3. Create a virtual environment:
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

4. Install the package in development mode:
   ```bash
   make install-dev
   ```

   Or manually:
   ```bash
   pip install -e .[dev,test,docs]
   pre-commit install
   ```

## Development Workflow

### Making Changes

1. Create a new branch for your feature or bugfix:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. Make your changes, following the coding standards below

3. Add or update tests for your changes

4. Run the test suite:
   ```bash
   make test
   ```

5. Commit your changes with a descriptive commit message:
   ```bash
   git commit -m "Add feature: your feature description"
   ```

6. Push to your fork and create a merge request

### Coding Standards

#### Code Style

- Follow PEP 8 style guidelines
- Use Black for code formatting (line length: 88 characters)
- Use isort for import sorting
- Add type hints for all functions and methods
- Write comprehensive docstrings in NumPy style

#### Code Quality Tools

Before submitting your changes, run:

```bash
make format     # Format code with black and isort
make lint       # Run flake8 and bandit security checks
make type-check # Run mypy type checking
make test       # Run all tests
```

#### Pre-commit Hooks

The project uses pre-commit hooks to automatically check code quality. These are installed when you run `make install-dev`. The hooks will run automatically on `git commit`.

### Testing

#### Running Tests

```bash
# Run all tests
make test

# Run only fast tests (skip slow integration tests)
make test-fast

# Run specific test file
pytest tests/test_specie.py -v

# Run with coverage
pytest tests/ --cov=mdinterface --cov-report=html
```

#### Test Categories

- **Unit tests**: Fast tests that test individual functions/methods
- **Integration tests**: Tests that verify component interactions
- **Slow tests**: Tests that take longer to run (marked with `@pytest.mark.slow`)

#### Writing Tests

- Write tests for all new functionality
- Use pytest fixtures for common test data
- Mark slow tests with `@pytest.mark.slow`
- Mark integration tests with `@pytest.mark.integration`
- Aim for good test coverage (>80%)

### Documentation

#### Docstring Style

Use NumPy-style docstrings:

```python
def example_function(param1: int, param2: str = "default") -> bool:
    """Short description of the function.

    Longer description if needed, explaining the purpose
    and behavior of the function.

    Parameters
    ----------
    param1 : int
        Description of param1
    param2 : str, default="default"
        Description of param2

    Returns
    -------
    bool
        Description of return value

    Raises
    ------
    ValueError
        When param1 is negative

    Examples
    --------
    >>> example_function(5, "test")
    True
    """
    pass
```

#### Type Hints

All functions should have proper type hints:

```python
from typing import List, Optional, Union
from numpy.typing import NDArray

def process_atoms(
    atoms: ase.Atoms,
    charges: Optional[NDArray[np.float64]] = None
) -> List[str]:
    """Process atoms and return element symbols."""
    pass
```

## Merge Request Process

1. **Fork and branch**: Create a feature branch from `main`
2. **Implement**: Make your changes with tests and documentation
3. **Test**: Ensure all tests pass and coverage is maintained
4. **Format**: Run code formatting and linting tools
5. **Commit**: Use clear, descriptive commit messages
6. **Push**: Push to your fork and create a merge request
7. **Review**: Address any feedback from code review
8. **Merge**: Maintainer will merge once approved

### Merge Request Checklist

- [ ] Code follows project style guidelines
- [ ] Tests added for new functionality
- [ ] All tests pass
- [ ] Documentation updated if needed
- [ ] Type hints added for new code
- [ ] No new security issues (bandit passes)
- [ ] Commit messages are clear and descriptive

## Code of Conduct

- Be respectful and inclusive
- Focus on constructive feedback
- Help newcomers and be patient with questions
- Follow scientific best practices for reproducible research

## Getting Help

- Open an issue for bugs or feature requests
- Start a discussion for questions about development
- Contact maintainers for sensitive issues

## Release Process

Releases are handled by maintainers:

1. Update version in `mdinterface/__init__.py`
2. Update `CHANGELOG.md`
3. Create a git tag: `git tag v1.x.x`
4. Push tag: `git push origin v1.x.x`
5. GitHub Actions will automatically build and publish to PyPI

## Development Tools

### Useful Commands

```bash
# Quick development cycle
make format && make lint && make type-check && make test-fast

# Full check before committing
make format && make lint && make type-check && make test

# Clean build artifacts
make clean

# Build distribution packages
make build
```

### IDE Setup

#### VS Code

Recommended extensions:
- Python
- Pylance
- Black Formatter
- isort
- GitLens

Settings (`.vscode/settings.json`):
```json
{
    "python.formatting.provider": "black",
    "python.linting.enabled": true,
    "python.linting.flake8Enabled": true,
    "python.linting.mypyEnabled": true
}
```

#### PyCharm

- Enable Black as the code formatter
- Configure flake8 and mypy as external tools
- Set line length to 88 characters

Thank you for contributing to mdinterface!