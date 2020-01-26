test:
	pytest
cov:
	pytest --cov-report html --cov=python_code python_tests/

