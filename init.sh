#!/bin/bash

root=$(cd $(dirname $0); pwd)

mkdir -p "$root/../data"
mkdir -p "$root/../data/MALDI-TOF"
mkdir -p "$root/../data/HPLC-MS"
mkdir -p "$root/../summary"
mkdir -p "$root/../summary/MALDI-TOF"
mkdir -p "$root/../summary/HPLC-MS"
mkdir -p "$root/../plot"
echo
read -p "Finished initializing directories."
