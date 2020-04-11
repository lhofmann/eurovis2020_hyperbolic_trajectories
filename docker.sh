#!/bin/bash
set -e

readonly image=eurovis2020-dht
readonly container=eurovis2020-dht
readonly cwd="$(dirname "$(readlink -f "$0")")"

# build base image
if [[ "$(docker images -q ${image} 2> /dev/null)" == "" ]]; then
  docker build --rm --network=host -t "${image}" -f "${cwd}/Dockerfile" "${cwd}"
fi

# start container
docker top ${container} >/dev/null 2>&1 || \
docker start ${container} >/dev/null 2>&1 || \
docker run -itd \
    --name ${container} \
    --volume="${cwd}:/mnt/source:ro" \
    ${image}

# build
docker exec ${container} cmake \
	-DVTK_PYTHONPATH=/usr/lib/python2.7/dist-packages \
    -B/tmp/build \
    -H/mnt/source
docker exec ${container} cmake --build /tmp/build --target install

# run example script
docker exec --workdir /tmp ${container} python /mnt/source/example.py

# copy result data to host
docker cp ${container}:/tmp/vector_field.vti .
docker cp ${container}:/tmp/candidate_line.vtp .
docker cp ${container}:/tmp/dht.vtp .
