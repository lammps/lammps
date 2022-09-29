for (( ; ; ))
do 
    python test_en.py
    # terminate loop if seg fault
    if [[ $? -eq 139 ]]; then 
        break 
    fi
done
