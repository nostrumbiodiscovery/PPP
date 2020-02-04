#Clean and reinstall requirements
pip uninstall -r requirements.txt --yes
pip uninstall PPP --yes
pip install .

#Clean and build
rm -r dist PPP.egg*
python setup.py sdist
twine upload  dist/*

