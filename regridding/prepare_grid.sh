#! /bin/bash
INFILE=${1}

echo $INFILE
sed -i 's^'"\["'^'""'^g' $INFILE
sed -i 's^'"\]"'^'""'^g' $INFILE
sed -i 's^'"', '"'^'""'^g' $INFILE
sed -i 's^'"'"'^'""'^g' $INFILE

