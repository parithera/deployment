#!/bin/bash

for img in $(docker compose config --images); do
  images="$images $img"
done

docker save -o services.img $images