# Usage: cat wn.xml | ./wn2txt.txt > input.txt
sed -rn 's/.*turl="([^"]*)".*/\1/p' | sort
