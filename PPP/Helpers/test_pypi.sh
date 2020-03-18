#Clean after build
pip uninstall -r requirements.txt --yes
pip uninstall PPPPele --yes
pip install PPPele
cd tests
pytest
cd ..


