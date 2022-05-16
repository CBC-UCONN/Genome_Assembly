#/bin/bash

echo "Long Read Coverage (taking genome size 300Mb) "
length=$(awk '{s+=$7}END{print s}' centrifuge_4514407.out)

coverage=$((${length}/300000000))
echo "Before centrifuge coverage = " ${coverage}"X"
echo -e "\n"


l2=$(awk '$2=="unclassified" {s+=$7}END{print s}' centrifuge_4514407.out)
c2=$((${l2}/300000000))

echo "Post centrifuge coverage = " ${c2}"X"




