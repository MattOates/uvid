.PHONY: help setup dev build test test-rust test-python lint fmt fmt-check clippy ruff ruff-check mypy clean docs docs-serve

help: ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2}'

setup: ## Install Rust tooling and sync Python dependencies
	rustup component add clippy rustfmt
	uv sync --all-extras

dev: ## Build the extension in-place (debug mode)
	maturin develop

build: ## Build the extension in-place (release mode)
	maturin develop --release

test: test-rust test-python ## Run all tests

test-rust: ## Run Rust tests
	cargo test --lib

test-python: build ## Run Python tests (builds release extension first)
	python -m pytest tests/python -v

lint: fmt-check clippy ruff-check mypy ## Run all linters

fmt: ## Format all code (Rust + Python)
	cargo fmt
	ruff format python tests
	ruff check --fix python tests

fmt-check: ## Check Rust formatting
	cargo fmt --check

clippy: ## Run Clippy
	cargo clippy -- -D warnings

ruff: ## Format and fix Python code
	ruff format python tests
	ruff check --fix python tests

ruff-check: ## Check Python formatting and lint rules
	ruff format --check python tests
	ruff check python tests

mypy: ## Type-check Python code
	mypy

docs: build ## Build the documentation site
	cargo doc --no-deps
	mkdir -p docs/_rust_api
	cp -r target/doc/* docs/_rust_api/
	python -m mkdocs build --strict

docs-serve: build ## Serve the documentation site locally
	python -m mkdocs serve

clean: ## Remove build artifacts
	cargo clean
	rm -rf site/ docs/_rust_api/ target/
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name '*.egg-info' -exec rm -rf {} + 2>/dev/null || true
