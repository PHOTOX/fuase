ENV_BASE_DIR=$(dirname $(readlink -f ${BASH_SOURCE[0]}))

export PYTHONPATH=$ENV_BASE_DIR:$PYTHONPATH
export PATH=$ENV_BASE_DIR/tools:$PATH

unset ENV_BASE_DIR
