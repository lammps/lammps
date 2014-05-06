

for ((i=0; i<60; i++)); do echo "$((i+1)) " `echo "$i*0.05" | bc` 0 0; done

echo 61 3.0 0 -5

for ((i=61; i<=4000; i++)); do echo "$((i+1)) " `echo "$i*0.05" | bc` `echo "($i-60)*0.5"|bc` -10; done

