#!/bin/bash

# convert pngs to tifs
convert figure2.png -compress zip -define tiff:alpha=associated figure2.tif
convert figure3.png -compress zip -define tiff:alpha=associated figure3.tif

