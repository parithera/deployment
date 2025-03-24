#!/bin/bash

for img in $(docker compose config --images); do
  if [[ "$img" == *codeclarity* || "$img" == *parithera* ]]; then
    echo "$img"
  fi
  images="$images $img"
done

docker save -o services.img $images