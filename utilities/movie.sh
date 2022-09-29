#!/bin/bash

x=1;
for i in $(ls -t -r ../plots/*png); do
 counter=$(printf %03d $x);
 ln -s "$i" ./img"$counter".png
 x=$(($x+1))
done

avconv -f image2 -r 20 -vb 4096 -i ./img%03d.png ./movie.mp4

rm *.png
