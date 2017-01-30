#!/bin/bash
aspell --lang=en create master ./Scripts/IKthesis.rws < ./Scripts/IKthesis-dict.txt
if [ x$1 == x"" ]; then
for a in $(ls -1 IKthesis-*.tex); do
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
echo $a
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
aspell  --add-extra-dicts=./Scripts/IKthesis.rws  -t --dont-tex-check-comments -d en_GB -c $a
done
exit
fi
echo "1111111111111"
aspell  --add-extra-dicts=./Scripts/IKthesis.rws  -t --dont-tex-check-comments -d en_GB -c $1
