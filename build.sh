#!/bin/bash



# copy git directory structure for stargraph to conda package location
cd ${SRC_DIR}
mkdir -p $PREFIX/main $PREFIX/aux $PREFIX/db $PREFIX/lib
cp -r bin/* $PREFIX/bin/
cp -r db/* $PREFIX/db/


# compile CNEFinder
git clone --branch main https://github.com/egluckthaler/starfish.git starfish
cd starfish/CNEFinder/
./pre-install.sh
make -f Makefile
mv cnef ../bin/
cd ..

# copy git directory structure for starfish to conda package location
cp -r aux/* $PREFIX/aux/
cp -r main/* $PREFIX/main/
cp -r bin/* $PREFIX/bin/
cp -r db/* $PREFIX/db/
cp -r lib/* $PREFIX/lib/
