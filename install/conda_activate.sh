# activate an openvisus env

#usage:
# source ~/bin/conda_activate.sh <conda_envname>

#echo "PS1: $PS1, \$1: $1"
if [ -z "$PS1" ] || [ -z "$1" ]; then
  echo "ERROR: must call this script using \"source ~/bin/conda_activate.sh <conda_envname>\")"
else
  echo "activating $1..."
  conda activate $1
  export VISUS_HOME=${CONDA_PREFIX}/OpenVisus
fi

