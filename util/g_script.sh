#!/bin/bash

cat ../gt.t
git add .
git commit -m "$1"
git push
