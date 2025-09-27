.PHONY: help install install-dev test test-fast lint format type-check clean build docs

help:
	@echo "Available commands:"
	@echo "  install      Install package in production mode"
	@echo "  install-dev  Install package in development mode with all extras"
	@echo "  test         Run full test suite"
	@echo "  test-fast    Run fast tests only (skip slow/integration tests)"
	@echo "  lint         Run linting checks"
	@echo "  format       Format code with black and isort"
	@echo "  type-check   Run type checking with mypy"
	@echo "  clean        Clean build artifacts"
	@echo "  build        Build distribution packages"
	@echo "  docs         Build documentation"

install:
	pip install -e .

install-dev:
	pip install -e .[dev,test,docs]
	pre-commit install

test:
	pytest tests/ -v --cov=mdinterface --cov-report=term-missing --cov-report=html

test-fast:
	pytest tests/ -v -m "not slow and not integration"

lint:
	flake8 mdinterface tests
	bandit -r mdinterface -f json -o bandit-report.json

format:
	black mdinterface tests
	isort mdinterface tests

type-check:
	mypy mdinterface --ignore-missing-imports

clean:
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	rm -rf .pytest_cache/
	rm -rf .mypy_cache/
	rm -rf htmlcov/
	rm -rf .coverage
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete

build: clean
	python -m build

docs:
	@echo "Documentation build not yet configured"
	@echo "Would run: sphinx-build -b html docs docs/_build/html"