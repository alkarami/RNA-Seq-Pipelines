#!/bin/bash
pip install --user umi_tools
export PY_USER_BIN=$(python -c 'import site; print(site.USER_BASE + "/bin")')
export PATH=$PY_USER_BIN:$PATH
for f1 in *_R1.fastq.gz
do
        umi_tools whitelist --stdin $f1 --extract-method=regex --bc-pattern="(?P<cell_1>.{8,12})(?P<discard_1>GAGTGATTGCTTGTGACGCCTT)(?P<cell_2>.{8})(?P<umi_1>.{6})T{3}.*" --set-cell-number=100 --log2stderr > "$f1_whitelist.txt" 
done