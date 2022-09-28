for (( ; ; ))
do 
    python test_en.py
    if [[ $? -eq 139 ]]; then 
        break 
    fi
done
