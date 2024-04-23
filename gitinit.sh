#!/bin/bash

#cat > temp << EOF

git config --global color.branch auto
git config --global color.diff   auto
git config --global color.status auto

git config --global user.name "Jeremy-Morere"
git config --global user.email jeremymorere57370@gmail.com

git config --global push.default current

git init

touch .gitignore

#EOF
#cat temp
#rm temp
