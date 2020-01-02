#!/bin/bash

rm -rf tiles

export GDAL_LIBJPEG_LARGEST_MEM_ALLOC=10691937

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
