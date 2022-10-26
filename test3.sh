#!/bin/bash
echo "-----start testing----"
PATH=/usr/bin:$PATH
startTime=`date +%Y%m%d-%H:%M:%S`
startTime_s=`date +%s`
for i in {1..100}
do
    echo "Test"${i}": "
    results=`./type3_2 | grep -E "Till Final IOcost:|TEST"`
    echo "$results"
done
endTime=`date +%Y%m%d-%H:%M:%S`
endTime_s=`date +%s`

sumTime=$[ $endTime_s - $startTime_s ]

echo "$startTime ---> $endTime" "Total:$sumTime seconds"
echo "-----finish testing----"
