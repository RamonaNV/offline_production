#!/bin/sh

# -------------------------------------------------------
# GridFTP wrapper script, based on an example by Heath
# Skarlupa: https://wiki.icecube.wisc.edu/index.php/GZK_Cluster
# -------------------------------------------------------

uid=51509

function print_usage () {
  echo "  options:"
  echo "    -t TARBALL  IceTray tarball to unpack"
  echo
  echo "    -d DEST     Destination directory: files will be copied to "
  echo "                gsiftp://data.icecube.wisc.edu/DEST/BASENAME"
  echo
  echo "    -i          Switch to input mode. All following arguments will be interpreted"
  echo "                as files to copy from gsiftp://data.icecube.wisc.edu"
  echo
  echo "    -o          Switch to output mode. All following arguments will be interpreted"
  echo "                as files to copy back to gsiftp://data.icecube.wisc.edu"
  echo
  echo "    --          End argument parsing. All following arguments will be"
  echo "                run as a script."
  exit 1
}

function link_jvm () {
  # Try to figure out where the heck JAVA_HOME is
  sdk=jdk1.5.0_21_x86_64
  for dir in /home/icecube/tools /data2/icecube/tools /usr/java; do
    if [[ -d $dir/$sdk ]]; then
      export JAVA_HOME=$dir/$sdk
      break
    fi
  done

  if [[ -z $JAVA_HOME ]]; then
    echo "Can't find JAVA_HOME!"
    exit 1
  fi

  for path in `find $JAVA_HOME -name libjvm.so`; do
    export JVM=`dirname $path`;
    break;
  done
  if [[ -z "$JVM" ]]; then
    echo 'could not find a jvm'
    exit 1
  fi
  
  JAVA_LIB=`dirname $JVM`/lib
  
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$JVM:$JAVA_LIB

  # ln -s $localjvm $I3_BUILD/lib/tools/
  
  echo "local jvm: $JVM"
  echo LD_LIBRARY_PATH: $LD_LIBRARY_PATH
  # ls -lh $I3_BUILD/lib/tools/$localjvm
}


if [[ "$#" == "0" ]]; then
  print_usage
fi

next="input"
while (( "$#" )); do
  # break out if we've seen "--"
  if [[ $next == "args" ]]; then
    break
  fi
  case $1 in
    "-d")
      shift;
      dest=$1
      shift;;
    "-t")
      shift;
      tarball=$1
      shift;;
    "-i")
      next="input";
      shift;;
    "-o")
      next="output";
      shift;;
    "--")
      next="args";
      shift;;
    *)
      case $next in 
        "input")
          inputs="$inputs $1";
          shift;;
        "output")
          outputs="$outputs $1";
          shift;;
        "args")
          break;;
      esac
  esac
done

if [[ -z $dest ]]; then
  echo "You must specify a destination!"
  print_usage
fi

if [[ -z $tarball ]]; then
  echo "You must specify a tarball!"
  print_usage
fi

# Unpack and set up GridFTP executables
tar xzf globus.tar.gz
# Location of CA Certs Directory
export X509_CERT_DIR=${_CONDOR_SCRATCH_DIR}/globus/certificates
export PATH=$PATH:${_CONDOR_SCRATCH_DIR}/globus/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${_CONDOR_SCRATCH_DIR}/globus/lib
# Location of my user proxy
export X509_USER_PROXY=${_CONDOR_SCRATCH_DIR}/x509up_u${uid}

# Transfer inputs
for infile in $inputs $tarball; do
  globus-url-copy -p 8 -vb gsiftp://data.icecube.wisc.edu/$infile file://${_CONDOR_SCRATCH_DIR}/`basename $infile`
done

# Unpack the tarball
tar xzf `basename $tarball`

export I3_BUILD=${_CONDOR_SCRATCH_DIR}/`basename $tarball .tar.gz`

env

# Run the script in an env-shell
$I3_BUILD/env-shell.sh $@

# Copy outputs back
for outfile in $outputs; do
  globus-url-copy -p 8 -vb file://${_CONDOR_SCRATCH_DIR}/$outfile gsiftp://data.icecube.wisc.edu/$dest/$outfile
done
