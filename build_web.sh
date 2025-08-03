mkdir -p build/web
cd build/web
emcmake cmake ../..
emmake make
mkdir -p ../../web/public
cp ../../dist/* ../../web/public
