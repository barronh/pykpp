for input in ~/Development/kpp/models/*.def
do
    python -m pykpp ${input} >& ${input}.log
done