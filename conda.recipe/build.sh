
git clone https://github.com/fabm-model/fabm.git
git clone https://github.com/gotm-model/code.git gotm
cd gotm && git submodule update --init --recursive && cd ..

mkdir build_pyfabm && cd build_pyfabm
cmake ../fabm/src/drivers/python -DFABM_ERSEM_BASE=.. -DCMAKE_INSTALL_PREFIX:PATH=$PREFIX 
# This makes sure pyfabm is built in the correct enviroment
# Will be changed once pyfabm install is updated.
sed -i -e 's/--user//g' cmake_install.cmake
make -j${CPU_COUNT}
make install DESTDIR=$PREFIX
cd ..

mkdir build_0d && cd build_0d
cmake -DCMAKE_INSTALL_PREFIX:PATH=$PREFIX ../fabm/src/drivers/0d -DGOTM_BASE=../gotm -DFABM_ERSEM_BASE=..
make -j${CPU_COUNT}
make install
cd ..

mkdir build_gotm && cd build_gotm
cmake -DCMAKE_INSTALL_PREFIX:PATH=$PREFIX ../gotm -DFABM_BASE=../fabm -DFABM_ERSEM_BASE=..
make -j${CPU_COUNT}
make install
cd ..

