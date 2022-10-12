# transparently fetch external files for a given package

fetch_potentials() {
    pdir="$1"
    shift

    type curl > /dev/null 2>&1 && have_curl=1 || have_curl=0
    type wget > /dev/null 2>&1 && have_wget=1 || have_wget=0
    if [ $have_curl -ne 1 ] && [ $have_wget -ne 1 ]
    then \
        echo "Need 'curl' or 'wget' to fetch external potential files"
        return
    fi

    while [ $# -gt 1 ]
    do \
        file=$1; sum=$2
        shift; shift

        echo ${sum} ${pdir}/${file} | md5sum -c - > /dev/null 2>&1 \
            && need_fetch=0 || need_fetch=1
        if [ ${need_fetch} -eq 1 ]
        then \
            url="https://download.lammps.org/potentials/${file}.${sum}"
            echo "Fetching external potential file ${file} from ${url}"
            if [ ${have_curl} ]
            then \
                curl -L -o ${pdir}/${file} ${url}
            elif [ ${have_wget} ]
            then \
                wget -O ${pdir}/${file} ${url}
            fi
        fi
   done
}

if [ -f "$1/potentials.txt" ]
then
    fetch_potentials "$1/../../potentials" `sed -e 's/#.*$//' "$1/potentials.txt"`
fi
