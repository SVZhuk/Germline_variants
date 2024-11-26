#!/bin/bash

function log() {
  local message="$1"
  echo "$(date +'%Y-%m-%d %H:%M:%S') : $message" | tee -a "$LOG_FILE"
}

function log_and_run() {
  local cmd="$*"
  log "Running command: $cmd"
  eval "$cmd"
}

function activate_conda() {
  log "Activating nf-core conda environment."
  # shellcheck disable=SC1090
  source "$CONDA_ENV_PATH"
  if ! conda activate "$CONDA_ENV_NAME"; then
    echo "ERROR: Failed to activate conda environment $CONDA_ENV_NAME"
    exit 1
  fi
}

function cleanup_docker() {
  log "Cleaning up docker images."
  docker system prune --all --force
  docker image prune --all --force
}
