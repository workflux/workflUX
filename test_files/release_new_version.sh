# update version in setup and init

# convert Readme.md to Readme.rst
pandoc -f markdown -t rst -o README.rst README.md

# make dist:
rm -rf dist/*
python setup.py sdist

# check dist
twine check dist/*

# upload to test pypi
twine upload --repository-url https://test.pypi.org/legacy/ dist/*

# upload to production pypi
twine upload dist/*