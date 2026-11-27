source /opt/foam/foam-extend-4.0/etc/bashrc >/dev/null 2>&1
export FOAM_USER_LIBBIN=$HOME/tukovic-fsi40/platform/lib
export FOAM_USER_APPBIN=$HOME/tukovic-fsi40/platform/bin
export PATH=$FOAM_USER_APPBIN:$PATH
export LD_LIBRARY_PATH=$FOAM_USER_LIBBIN:$LD_LIBRARY_PATH
mkdir -p $FOAM_USER_LIBBIN $FOAM_USER_APPBIN
