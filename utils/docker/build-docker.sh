#!/bin/bash

# Utility script to start the Docker build process.

set -x
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

GIT_DESCRIBE=$(git describe --tags | cut -d - -f 1)
GIT_TAG=${GIT_TAG-$GIT_DESCRIBE}
DOCKER_VERSION=$(echo $GIT_TAG | sed -e 's/^v//')

ORG=bihealth
REPO=mehari

docker build . \
    --file utils/docker/Dockerfile \
    --pull \
    -t ghcr.io/$ORG/$REPO:$DOCKER_VERSION
