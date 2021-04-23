#using pyCRAC on katana

#pyCRAC is installed locally using a virual environment. Itâ€™s installed in the dir

my_python_dir/

#The installation required the I remove python from .bashrc setup so that 
#python/2.7.10 can be loaded - the virualenv command has changed in python3
# so after restarting the katana session

mkdir my_python_dir
cd my_python_dir
module load python/2.7.15
virtualenv --system-site-packages my_python_env
. my_python_env/bin/activate           #The environment needs to be activated
wget https://bitbucket.org/sgrann/pycrac/get/cc0dcde5ea3f.zip
unzip cc0dcde5ea3f.zip
cd sgrann-pycrac-cc0dcde5ea3f/
python setup.py install

#have had to install some additional packages: 

wget http://www.parallelpython.com/downloads/pp/pp-1.6.4.zip
unzip pp-1.6.4.zip
cd pp-1.6.4
python setup.py install

#pip seems to be a better way to install 

pip install greenlet
pip install curtsies
pip install requests
pip install pygments
pip install ruffus --upgrade

#py ReadCoutners seems to work now
#to exit the session

deactivate