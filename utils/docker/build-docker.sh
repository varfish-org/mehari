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

GIT_DEPTH=$(($(git rev-list HEAD ^$(git describe --abbrev=0 --tags) --count) + 1))
GIT_URL=https://github.com/$ORG/$REPO.git

docker build . \
    --file utils/docker/Dockerfile \
    --build-arg git_treeish=$GIT_TAG \
    --build-arg git_url=$GIT_URL \
    --pull \
    -t ghcr.io/$ORG/$REPO:$DOCKER_VERSION
