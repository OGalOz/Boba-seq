
if [ -z "$1" ]
then 
    echo "MISSING Arg!!"
else
    cp -r $1/* ~/transfer/NewPlots
fi

