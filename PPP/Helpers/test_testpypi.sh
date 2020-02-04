#Clean after build
pip uninstall -r requirements.txt --yes
pip uninstall PPP --yes
pip install -r requirements.txt
pip install --index-url https://test.pypi.org/simple/ PPPELE
cd tests
pytest


