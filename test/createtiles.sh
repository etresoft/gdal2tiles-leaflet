#!/bin/bash

rm -rf tiles

export GDAL_ALLOW_LARGE_LIBJPEG_MEM_ALLOC=1

case $1 in
  mpz)
    ../gdal2tiles-multiprocess.py -l -p raster -z 0-5 -w none karta.jpg tiles
    ;;
  mp)
    ../gdal2tiles-multiprocess.py -l -p raster -w none karta.jpg tiles
    ;;
  z)
    ../gdal2tiles.py karta.jpg -z 0-5 tiles
    ;;
  *)
    ../gdal2tiles.py karta.jpg tiles
    ;;
esac
