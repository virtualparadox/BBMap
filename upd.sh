#!/bin/sh

# extract BBMap tgz to ./t/ first.
tar -xzv -C ./t --strip-components=1 -f $1

rm -fr current && mv t/current .
rm -fr docs && mv t/docs .
rm -fr resources && mv t/resources .
rm -fr sh && mkdir sh && mv t/*.sh sh/
mv t/build.xml . && mv t/license.txt .
mv t/README.md tt/

sed -i.bak 's/^# *//g' tt/README.md
cat tt/README.md
ln -s ../current sh/
ls -la t tt sh/current

git add sh docs resources current

echo git commit . -m "'Extract Version 3?.?? from `basename $1`'"
# git tag -a v31.56 -F tt/README.md
# git push && git push --tags
