#Clean and reinstall requirements
pip uninstall -r requirements.txt --yes
pip uninstall PPP --yes
pip install .

#fast test
python -m PPP.main -h
cd tests
rm 1l63_proc_processed.pdb
pytest
cd ..

#Clean and build
rm -r dist PPP.egg*
python setup.py sdist
twine upload  --repository-url https://test.pypi.org/legacy/ dist/*
#twine upload  dist/*

