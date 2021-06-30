# Aperture Obstruction Calculator (MOCCA)

## About

MOCCA (which stands for **M**ick's aperture **O**bstruction **C**al**C**ulat**O**r) allows you to compute the percentage obstruction of the aperture of a telescope aperture, on an equatorial mount, by a hemispherical dome using basic ray tracing techniques. Particularly designed for the [Gratama telescope](https://www.rug.nl/research/kapteyn/sterrenwacht/gratama?lang=en) and dome of the [Blaauw Observatory](https://www.rug.nl/research/kapteyn/sterrenwacht/). However, by changing the properties in `config.ini` it could be adapted to any telescope on an equatorial mount.

## Usage

With the script `mocca.py`, one can determine the percentage obstruction of the telescope aperture by the dome. Add the `-h` argument to display additional options, such as what aperture (in the case of the Gratama telescope: telescope/guider/finder) and more.