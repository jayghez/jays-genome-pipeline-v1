.PHONY: setup test lint sample check-tools clean

setup:
	python -m venv .venv
	. .venv/bin/activate && python -m pip install --upgrade pip
	. .venv/bin/activate && python -m pip install -e ".[dev]"

test:
	python -m pytest

check-tools:
	./scripts/check_tools.sh

sample:
	genome-pipeline run --vcf data/sample/HG001_chr20_subset.vcf.gz --objective hereditary_cancer

clean:
	rm -rf .pytest_cache .coverage htmlcov
