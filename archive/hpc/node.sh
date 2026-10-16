#!/bin/sh




read -p "Enter node no. you want to submit job in chandra!: " no
var="bose"$no

#grep -n "bose" node_target.sh
awk -v cuv1="bose" -v cuv2="$var" 'NR%2==0{gsub(cuv1,cuv2)}1' node_target.sh > node_target.sh.tmp #&& mv node_target.sh.tmp node_target.sh


#decision=true
read -p "Editing done!. You want to submit? " decision

if [ "$decision" = yes  ]
then
   chmod +x node_target.sh.tmp submit.sh
   ./node_target.sh.tmp && ./submit.sh
   echo "Job Submitted to node:   "$var 
else
   echo "ok! Never mind!"
fi
