#!/bin/bash

# ----------------
# Copyright 2016 Jaroslaw Glowacki
# jarek (dot) glowacki (at) gmail (dot) com
#
# This file is part of cgDNArecon.
#
# cgDNArecon is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# cgDNArecon is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with cgDNArecon.  If not, see <http://www.gnu.org/licenses/>.
# ----------------
#
# This script provides a reference of how the cgDNAmc project can be built
# together with the two dependencies: a,gebra3d and cgDNArecon

# A helper macro for building a project
# (see below)
build_project() {
	proj_name=$1
	echo
	echo "*******************************************************************************"
	echo "Building $proj_name"
	cd $proj_name
	make
	if [ $? != 0 ]; then
		echo "Building $proj_name failed!" >&2
		echo "(see $proj_name/Makefile.config)" >&2
		exit 2
	fi
	cd -
}

echo "*******************************************************************************"
echo "Unpacking source files"

# Try unpacking all source files
unzip algebra3d.zip \
&& unzip cgDNArecon.zip \
&& unzip cgDNAmc.zip

if [ $? != 0 ]; then
	echo "Could not unpack source files!" >&2
	echo "Please place:" >&2
	echo "* algebra3d.zip" >&2
	echo "* cgDNArecon.zip" >&2
	echo "* cgDNAmc.zip" >&2
	echo "in the current directory" >&2
	exit 1
fi

echo
echo "*******************************************************************************"
echo "Setting up Makefile.config s"
# Make changes in Makefile.config of:
# * cgDNArecon to set path to algebra3d
sed -i.orig 's/\(ALG3D_DIR[ \t]*=\).*/\1 \..\/algebra3d/' cgDNArecon/Makefile.config
# * cgDNAmc to set paths to algebra3d and cgDNArecon
sed -i.orig 's/\(ALG3D_DIR[ \t]*=\).*/\1 \..\/algebra3d/' cgDNAmc/Makefile.config
sed -i.orig 's/\(CGDNA_RECON_DIR[ \t]*=\).*/\1 \..\/cgDNArecon/' cgDNAmc/Makefile.config

# Everything should be ready for building
# If the build phase fails some changes migh be necessary for the target
# platform. Adjustments should be made in Makefile.config files in the failing
# project. Most certainly building algebra3d and cgDNArecon causes no trouble.
#
# cdDNAmc, however, requires BLAS and LAPACK which may be installed in
# in a non-standard way (e.g. partucular implementations like OpenBLAS)
# 
# It is important to build the projects in the provided order as the subsequent
# depend on the preceeding ones
build_project algebra3d
build_project cgDNArecon
build_project cgDNAmc

